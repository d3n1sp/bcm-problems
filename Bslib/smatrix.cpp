//------------------------------------------------------------------------------------------------
// File: smatrix.cpp
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

#include "smatrix.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrix: Memory allocation zero data constructor
//========================================================================================
CSMatrix::CSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CSMatrix_00";

	m = 0;
	n = 0;
	nsupr = 0;
	nsupc = 0;
	nlist = 0;
	nlist2 = 0;
	nlist3 = 0;
	nzja = 0;
	nzja2 = 0;
	nzja3 = 0;
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) {
		list[i] = i;
	};
	list2 = new int [nlist2];
	if (!list2) MemoryFail (funcname);
	for (i=0;i<nlist2;i++) {
		list2[i] = i;
	};
	list3 = new int [nlist3];
	if (!list3) MemoryFail (funcname);
	for (i=0;i<nlist3;i++) {
		list3[i] = i;
	};
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	for (i=0;i<=nlist;i++) {
		ia[i] = 0;
	};
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		ja[i] = 0;
	};
	ja2 = new int [nzja2];
	if (!ja2) MemoryFail (funcname);
	for (i=0;i<nzja2;i++) {
		ja2[i] = 0;
	};
	ja3 = new int [nzja3];
	if (!ja3) MemoryFail (funcname);
	for (i=0;i<nzja3;i++) {
		ja3[i] = 0;
	};
};

// Author: Kharchenko S.A.
// CSMatrix: Memory allocation zero data constructor
//========================================================================================
CSMatrix::CSMatrix (int _nlist, int _nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrix_01";

	m = _nlist;
	n = _nlist;
	nsupr = _nlist;
	nsupc = _nlist;
	nlist = _nlist;
	nlist2 = 0;
	nlist3 = 0;
	nzja = _nzja;
	nzja2 = 0;
	nzja3 = 0;
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) {
		list[i] = i;
	};
	list2 = new int [nlist2];
	if (!list2) MemoryFail (funcname);
	for (i=0;i<nlist2;i++) {
		list2[i] = i;
	};
	list3 = new int [nlist3];
	if (!list3) MemoryFail (funcname);
	for (i=0;i<nlist3;i++) {
		list3[i] = i;
	};
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	for (i=0;i<=nlist;i++) {
		ia[i] = 0;
	};
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		ja[i] = 0;
	};
	ja2 = new int [nzja2];
	if (!ja2) MemoryFail (funcname);
	for (i=0;i<nzja2;i++) {
		ja2[i] = 0;
	};
	ja3 = new int [nzja3];
	if (!ja3) MemoryFail (funcname);
	for (i=0;i<nzja3;i++) {
		ja3[i] = 0;
	};
};

// Author: Kharchenko S.A.
// CSMatrix: Memory allocation zero data constructor
//========================================================================================
CSMatrix::CSMatrix (int _nlist, int _nlist2, int _nzja, int _nzja2) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrix_02";

	m = _nlist;
	n = _nlist;
	nsupr = _nlist;
	nsupc = _nlist;
	nlist = _nlist;
	nlist2 = _nlist2;
	nlist3 = 0;
	nzja = _nzja;
	nzja2 = _nzja2;
	nzja3 = 0;
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) {
		list[i] = i;
	};
	list2 = new int [nlist2];
	if (!list2) MemoryFail (funcname);
	for (i=0;i<nlist2;i++) {
		list2[i] = i;
	};
	list3 = new int [nlist3];
	if (!list3) MemoryFail (funcname);
	for (i=0;i<nlist3;i++) {
		list3[i] = i;
	};
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	for (i=0;i<=nlist;i++) {
		ia[i] = 0;
	};
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		ja[i] = 0;
	};
	ja2 = new int [nzja2];
	if (!ja2) MemoryFail (funcname);
	for (i=0;i<nzja2;i++) {
		ja2[i] = 0;
	};
	ja3 = new int [nzja3];
	if (!ja3) MemoryFail (funcname);
	for (i=0;i<nzja3;i++) {
		ja3[i] = 0;
	};
};

// Author: Kharchenko S.A.
// CSMatrix: Destructor
//========================================================================================
CSMatrix::~CSMatrix () { // Destructor
//	std::cout << " On entry to CSMatrix destructor " << std::endl;
	m     = 0;
	n     = 0;
	nsupr = 0;
	nsupc = 0;
	nlist = 0;
	nlist2 = 0;
	nlist3 = 0;
	nzja  = 0;
	nzja2  = 0;
	nzja3  = 0;
	delete [] list;
	delete [] list2;
	delete [] list3;
	delete [] ia;
	delete [] ja;
	delete [] ja2;
	delete [] ja3;
//	std::cout << " On return from CSMatrix destructor " << std::endl;
};

// Author: Kharchenko S.A.
// CSMatrix: Copy constructor
//========================================================================================
CSMatrix::CSMatrix (const CSMatrix &_aa) { // Copy constructor

	const char *funcname = "CSMatrix_copy";

	m = _aa.m;
	n = _aa.n;
	nsupr = _aa.nsupr;
	nsupc = _aa.nsupc;
	nlist = _aa.nlist;
	nlist2 = _aa.nlist2;
	nlist3 = _aa.nlist3;
	nzja = _aa.nzja;
	nzja2 = _aa.nzja2;
	nzja3 = _aa.nzja3;
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) {
		list[i] = _aa.list[i];
	};
	list2 = new int [nlist2];
	if (!list2) MemoryFail (funcname);
	for (i=0;i<nlist2;i++) {
		list2[i] = _aa.list2[i];
	};
	list3 = new int [nlist3];
	if (!list3) MemoryFail (funcname);
	for (i=0;i<nlist3;i++) {
		list3[i] = _aa.list3[i];
	};
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	for (i=0;i<=nlist;i++) {
		ia[i] = _aa.ia[i];
	};
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		ja[i] = _aa.ja[i];
	};
	ja2 = new int [nzja2];
	if (!ja2) MemoryFail (funcname);
	for (i=0;i<nzja2;i++) {
		ja2[i] = _aa.ja2[i];
	};
	ja3 = new int [nzja3];
	if (!ja3) MemoryFail (funcname);
	for (i=0;i<nzja3;i++) {
		ja3[i] = _aa.ja3[i];
	};
};

// Author: Kharchenko S.A.
// CSMatrix: Equality operator
//========================================================================================
CSMatrix &CSMatrix::operator= (const CSMatrix &_aa) { // Equality operator

	const char *funcname = "CSMatrix_=";

	delete [] list;
	delete [] list2;
	delete [] list3;
	delete [] ia;
	delete [] ja;
	delete [] ja2;
	delete [] ja3;
	m = _aa.m;
	n = _aa.n;
	nsupr = _aa.nsupr;
	nsupc = _aa.nsupc;
	nlist = _aa.nlist;
	nlist2 = _aa.nlist2;
	nlist3 = _aa.nlist3;
	nzja = _aa.nzja;
	nzja2 = _aa.nzja2;
	nzja3 = _aa.nzja3;
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) {
		list[i] = _aa.list[i];
	};
	list2 = new int [nlist2];
	if (!list2) MemoryFail (funcname);
	for (i=0;i<nlist2;i++) {
		list2[i] = _aa.list2[i];
	};
	list3 = new int [nlist3];
	if (!list3) MemoryFail (funcname);
	for (i=0;i<nlist3;i++) {
		list3[i] = _aa.list3[i];
	};
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	for (i=0;i<=nlist;i++) {
		ia[i] = _aa.ia[i];
	};
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	for (i=0;i<nzja;i++) {
		ja[i] = _aa.ja[i];
	};
	ja2 = new int [nzja2];
	if (!ja2) MemoryFail (funcname);
	for (i=0;i<nzja2;i++) {
		ja2[i] = _aa.ja2[i];
	};
	ja3 = new int [nzja3];
	if (!ja3) MemoryFail (funcname);
	for (i=0;i<nzja3;i++) {
		ja3[i] = _aa.ja3[i];
	};
	return *this;
};

// Author: Kharchenko S.A.
// CSMatrix: Compute sparsity of the 3D matrix
//========================================================================================
CSMatrix CSMatrix::Sparsity3D (int _nx, int _ny, int _nz) { // Compute sparsity of the 3D matrix

	const char *funcname = "Sparsity3D";

// Allocate the memory

	int ntot = _nx*_ny*_nz;
	int nxy = _nx*_ny;

	CSMatrix as (ntot,ntot*7);

	as.ia[0] = 0;

	int i;

	for (i=1;i<=ntot;i++) as.ia[i] = 1;

	int ix, iy, iz, ind;

	for (iz=0;iz<_nz;iz++) {
		for (iy=0;iy<_ny;iy++) {
			for (ix=0;ix<_nx;ix++) {
				ind = iz*nxy+iy*_nx+ix;
				if (ix>0) as.ia[ind+1]++;
				if (ix<_nx-1) as.ia[ind+1]++;
				if (iy>0) as.ia[ind+1]++;
				if (iy<_ny-1) as.ia[ind+1]++;
				if (iz>0) as.ia[ind+1]++;
				if (iz<_nz-1) as.ia[ind+1]++;
			};
		};
	};

	for (i=0;i<ntot;i++) as.ia[i+1] = as.ia[i]+as.ia[i+1];

	int *iptr;

	iptr = new int [ntot];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<ntot;i++) iptr[i] = as.ia[i];

	int k;

	for (i=0;i<ntot;i++) {
		k = iptr[i];
		as.ja[k] = i;
		iptr[i]++;
	};

	int ind1;

	for (iz=0;iz<_nz;iz++) {
		for (iy=0;iy<_ny;iy++) {
			for (ix=0;ix<_nx;ix++) {
				ind = iz*nxy+iy*_nx+ix;
				if (ix>0) {
					k = iptr[ind];
					ind1 = iz*nxy+iy*_nx+(ix-1);
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (ix<_nx-1) {
					k = iptr[ind];
					ind1 = iz*nxy+iy*_nx+(ix+1);
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iy>0) {
					k = iptr[ind];
					ind1 = iz*nxy+(iy-1)*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iy<_ny-1) {
					k = iptr[ind];
					ind1 = iz*nxy+(iy+1)*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iz>0) {
					k = iptr[ind];
					ind1 = (iz-1)*nxy+iy*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iz<_nz-1) {
					k = iptr[ind];
					ind1 = (iz+1)*nxy+iy*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				qsort (as.ja+as.ia[ind],as.ia[ind+1]-as.ia[ind],sizeof(int),compint);
			};
		};
	};

	for (i=0;i<ntot;i++) as.list[i] = i;

	as.m = ntot;
	as.n = ntot;
	as.nsupr = ntot;
	as.nsupc = ntot;
	as.nlist = ntot;
	as.nzja = as.ia[ntot];

// Free work arrays

	delete [] iptr;

// Return the result

	return as;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute sparsity of the 3D matrix
//========================================================================================
CSMatrix CSMatrix::Sparsity3D (int _nx, int _ny, int _nz, int _ibegz, int _iendz) { // Compute sparsity of the 3D matrix

	const char *funcname = "Sparsity3D";

// Allocate the memory

	int nzloc = _iendz-_ibegz+1;

	int ntot = _nx*_ny*_nz;
	int nloc = _nx*_ny*nzloc;
	int nxy = _nx*_ny;

	CSMatrix as (nloc,nloc*7);

	as.ia[0] = 0;

	int i;

	for (i=1;i<=nloc;i++) as.ia[i] = 1;

	int ix, iy, iz, ind;

	for (iz=_ibegz;iz<=_iendz;iz++) {
		for (iy=0;iy<_ny;iy++) {
			for (ix=0;ix<_nx;ix++) {
				ind = (iz-_ibegz)*nxy+iy*_nx+ix;
				if (ix>0) as.ia[ind+1]++;
				if (ix<_nx-1) as.ia[ind+1]++;
				if (iy>0) as.ia[ind+1]++;
				if (iy<_ny-1) as.ia[ind+1]++;
				if (iz>0) as.ia[ind+1]++;
				if (iz<_nz-1) as.ia[ind+1]++;
			};
		};
	};

	for (i=0;i<nloc;i++) as.ia[i+1] = as.ia[i]+as.ia[i+1];

	int *iptr;

	iptr = new int [nloc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nloc;i++) iptr[i] = as.ia[i];

	int k, jj;

	for (iz=_ibegz;iz<=_iendz;iz++) {
		for (iy=0;iy<_ny;iy++) {
			for (ix=0;ix<_nx;ix++) {
				ind = (iz-_ibegz)*nxy+iy*_nx+ix;
				jj = iz*nxy+iy*_nx+ix;
				k = iptr[ind];
				as.ja[k] = jj;
				iptr[ind]++;
			};
		};
	};

	int ind1;

	for (iz=_ibegz;iz<=_iendz;iz++) {
		for (iy=0;iy<_ny;iy++) {
			for (ix=0;ix<_nx;ix++) {
				ind = (iz-_ibegz)*nxy+iy*_nx+ix;
				jj = iz*nxy+iy*_nx+ix;
				if (ix>0) {
					k = iptr[ind];
					ind1 = iz*nxy+iy*_nx+(ix-1);
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (ix<_nx-1) {
					k = iptr[ind];
					ind1 = iz*nxy+iy*_nx+(ix+1);
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iy>0) {
					k = iptr[ind];
					ind1 = iz*nxy+(iy-1)*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iy<_ny-1) {
					k = iptr[ind];
					ind1 = iz*nxy+(iy+1)*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iz>0) {
					k = iptr[ind];
					ind1 = (iz-1)*nxy+iy*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				if (iz<_nz-1) {
					k = iptr[ind];
					ind1 = (iz+1)*nxy+iy*_nx+ix;
					as.ja[k] = ind1;
					iptr[ind]++;
				};
				qsort (as.ja+as.ia[ind],as.ia[ind+1]-as.ia[ind],sizeof(int),compint);
			};
		};
	};

	for (iz=_ibegz;iz<=_iendz;iz++) {
		for (iy=0;iy<_ny;iy++) {
			for (ix=0;ix<_nx;ix++) {
				ind = (iz-_ibegz)*nxy+iy*_nx+ix;
				jj = iz*nxy+iy*_nx+ix;
				as.list[ind] = jj;
			};
		};
	};

	as.m = ntot;
	as.n = ntot;
	as.nsupr = ntot;
	as.nsupc = ntot;
	as.nlist = nloc;
	as.nzja = as.ia[nloc];

// Free work arrays

	delete [] iptr;

// Return the result

	return as;

};

// Author: Kharchenko S.A.
// CSMatrix: Super sparsify the matrix
//========================================================================================
void CSMatrix::SuperSparsify () { // Super sparsify the matrix

	const char *funcname = "SuperSparsify";

	int *listloc, *ialoc;

	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);

	int i;
	int nlistnew = 0;

	ialoc[0] = 0;

	for (i=0;i<nlist;i++) {
		if (ia[i+1]>ia[i]) {
			listloc[nlistnew] = list[i];
			ialoc[nlistnew+1] = ia[i+1];
			nlistnew++;
		};
	};

	for (i=0;i<nlistnew;i++) {
		list[i] = listloc[i];
	};

	for (i=0;i<=nlistnew;i++) {
		ia[i] = ialoc[i];
	};

	nlist = nlistnew;

	delete [] listloc;
	delete [] ialoc;

};

// Author: Kharchenko S.A.
// CSMatrix: Super sparsify the matrix
//========================================================================================
void CSMatrix::SuperSparsify2Index () { // Super sparsify the matrix

	const char *funcname = "SuperSparsify2Index";

	int *listloc, *list2loc, *ialoc;

	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlist];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);

	int i;
	int nlistnew = 0;

	ialoc[0] = 0;

	for (i=0;i<nlist;i++) {
		if (ia[i+1]>ia[i]) {
			listloc[nlistnew] = list[i];
			list2loc[nlistnew] = list2[i];
			ialoc[nlistnew+1] = ia[i+1];
			nlistnew++;
		};
	};

	for (i=0;i<nlistnew;i++) {
		list[i] = listloc[i];
		list2[i] = list2loc[i];
	};

	for (i=0;i<=nlistnew;i++) {
		ia[i] = ialoc[i];
	};

	nlist = nlistnew;

	delete [] listloc;
	delete [] list2loc;
	delete [] ialoc;

};

// Author: Kharchenko S.A.
// CSMatrix: Shift the sparsity
//========================================================================================
void CSMatrix::ShiftSparsity (bool _is_ja_shift, int _ishift) { // Shift the sparsity

	int i;

	for (i=0;i<nlist;i++) list[i] += _ishift;

	if (_is_ja_shift) {
		for (i=0;i<nzja;i++) ja[i] += _ishift;
	};

};

// Author: Kharchenko S.A.
// CSMatrix: Add two matrices in super sparse format
//========================================================================================
CSMatrix CSMatrix::operator+ (const CSMatrix &_mtr2) const { // Add two matrices in super sparse format

	const char *funcname = "CSMatrix_+";

// Create the list of rows

	int nlistmax = nlist+_mtr2.nlist;

	int *listloc, *ialoc;

	listloc = new int [nlistmax];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlistmax+1];
	if (!ialoc) MemoryFail (funcname);

	int isup1, isup2;

	int nlistloc = 0;

	int iend1 = nlist-1;
	int iend2 = _mtr2.nlist-1;

	int ip1 = 0;
	int ip2 = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			isup1 = list[ip1];
			isup2 = _mtr2.list[ip2];
			if (isup1 == isup2) {
				listloc[nlistloc] = isup1;
				ip1++;
				ip2++;
			} else if (isup1 < isup2) {
				listloc[nlistloc] = isup1;
				ip1++;
			} else {
				listloc[nlistloc] = isup2;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			isup1 = list[ip1];
			listloc[nlistloc] = isup1;
			ip1++;
		} else {
			isup2 = _mtr2.list[ip2];
			listloc[nlistloc] = isup2;
			ip2++;
		};
		nlistloc++;
	};

// Fill ja array

	int nzjamax = nzja + _mtr2.nzja;

	int *jaloc;

	jaloc = new int [nzjamax];
	if (!jaloc) MemoryFail (funcname);

	int ilist, irow, jj1, jj2, jp1, jp2;
	int irow1, irow2, jend1, jend2, nz, iblk1, iblk2;

	ialoc[0] = 0;
	nz = 0;

	ip1 = 0;
	ip2 = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		irow = listloc[ilist];
		iblk1 = 0;
		iblk2 = 0;
		if (ip1 <= iend1) {
			irow1 = list[ip1];
			if (irow1 != irow) iblk1 = -1;
		} else {
			iblk1 = -1;
		};
		if (ip2 <= iend2) {
			irow2 = _mtr2.list[ip2];
			if (irow2 != irow) iblk2 = -1;
		} else {
			iblk2 = -1;
		};
		if (iblk1 >= 0 && iblk2 >= 0) {
			jend1 = ia[ip1+1]-1;
			jend2 = _mtr2.ia[ip2+1]-1;
			jp1 = ia[ip1];
			jp2 = _mtr2.ia[ip2];
			while (jp1 <= jend1 || jp2 <= jend2) {
				if (jp1 <= jend1 && jp2 <= jend2) {
					jj1 = ja[jp1];
					jj2 = _mtr2.ja[jp2];
					if (jj1 == jj2) {
						jaloc[nz] = jj1;
						nz++;
						jp1++;
						jp2++;
					} else if (jj1 < jj2) {
						jaloc[nz] = jj1;
						nz++;
						jp1++;
					} else if (jj1 > jj2) {
						jaloc[nz] = jj2;
						nz++;
						jp2++;
					};
				} else if (jp1 <= jend1) {
					jaloc[nz] = ja[jp1];
					nz++;
					jp1++;
				} else if (jp2 <= jend2) {
					jaloc[nz] = _mtr2.ja[jp2];
					nz++;
					jp2++;
				};
			};
			ip1++;
			ip2++;
		} else if (iblk1 >= 0) {
			for (jp1=ia[ip1];jp1<ia[ip1+1];jp1++) {
				jaloc[nz] = ja[jp1];
				nz++;
			};
			ip1++;
		} else {
			for (jp2=_mtr2.ia[ip2];jp2<_mtr2.ia[ip2+1];jp2++) {
				jaloc[nz] = _mtr2.ja[jp2];
				nz++;
			};
			ip2++;
		};
		ialoc[ilist+1] = nz;
	};

// Store computed data

	CSMatrix mtr;

	delete [] mtr.list;
	delete [] mtr.ia;
	delete [] mtr.ja;

	mtr.list  = listloc;
	mtr.ia    = ialoc;
	mtr.ja    = jaloc;

	int ntot = _mtr2.m;

	mtr.m      = ntot;
	mtr.n      = ntot;
	mtr.nsupr  = ntot;
	mtr.nsupc  = ntot;
	mtr.nlist  = nlistloc;
	mtr.nzja   = nz;

	return mtr;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute diagonal matrix 
//========================================================================================
CSMatrix CSMatrix::DiagonalMtr (int _n) { // Compute diagonal matrix 

//	const char *funcname = "DiagonalMtr";

	CSMatrix asp (_n,_n);

	int i;

	asp.ia[0] = 0;

	int nz = 0;

	for (i=0;i<_n;i++) {
		asp.list[i] = i;
		asp.ja[nz] = i;
		nz++;
		asp.ia[i+1] = nz;
	};

	asp.nzja = nz;

	return asp;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute tridiagonal matrix 
//========================================================================================
CSMatrix CSMatrix::TridiagonalMtr (int _n) { // Compute tridiagonal matrix 

//	const char *funcname = "TridiagonalMtr";

	CSMatrix asp (_n,3*_n);

	int i;

	asp.ia[0] = 0;

	int nz = 0;

	for (i=0;i<_n;i++) {
		asp.list[i] = i;
		if (i>0) {
			asp.ja[nz] = i-1;
			nz++;
		};
		asp.ja[nz] = i;
		nz++;
		if (i<_n-1) {
			asp.ja[nz] = i+1;
			nz++;
		};
		asp.ia[i+1] = nz;
	};

	asp.nzja = nz;

	return asp;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed matrix
//========================================================================================
CSMatrix CSMatrix::TranspMtr () const { // Compute transposed matrix 

	const char *funcname = "TranspMtr";

	CSMatrix at(nsupr,nzja);

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

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	int k;
	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			pat[jj]++;
			k = pat[jj];
			at.ja[k-1] = i;
		};
	};

	for (i=0;i<nsupr;i++) at.list[i] = i;

	delete [] pat;

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupc = nsupc;
	at.nsupr = nsupr;
	at.nlist = nlist;
	at.nzja = nzja;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed to the rectangular matrix stored by columns
//========================================================================================
CSMatrix CSMatrix::TranspMtrRect () const { // Compute transposed to the rectangular matrix stored by columns

	const char *funcname = "TranspMtrRect";

	CSMatrix at(nsupr,nzja);

	int i, j, jj;

	for (i=0;i<=nsupr;i++) at.ia[i] = 0;

	for (i=0;i<nsupc;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			at.ia[jj+1]++;
		};
	};

	for (i=0;i<nsupr;i++) at.ia[i+1] += at.ia[i];

	int *pat;

	pat = new int [nsupr];

	if (!pat) MemoryFail (funcname);

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];
	
	int k;
	for (i=0;i<nsupc;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			pat[jj]++;
			k = pat[jj];
			at.ja[k-1] = i;
		};
	};

	for (i=0;i<nsupr;i++) at.list[i] = i;

	delete [] pat;

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupc = nsupc;
	at.nsupr = nsupr;
	at.nlist = nsupr;
	at.nzja = nzja;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed to the rectangular matrix stored by rows
//========================================================================================
CSMatrix CSMatrix::TranspMtrRectRows () const { // Compute transposed to the rectangular matrix stored by rows

	const char *funcname = "TranspMtrRectRows";

	int nsupmx = nsupc;
	if (nsupr > nsupmx) nsupmx = nsupr;

	CSMatrix at(nsupmx,nzja);

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

	for (i=0;i<nsupc;i++) pat[i] = at.ia[i];
	
	int k;
	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			pat[jj]++;
			k = pat[jj];
			at.ja[k-1] = i;
		};
	};

	for (i=0;i<nsupc;i++) at.list[i] = i;

	delete [] pat;

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupc = nsupc;
	at.nsupr = nsupr;
	at.nlist = nsupc;
	at.nzja = nzja;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed to the rectangular matrix stored by columns
//========================================================================================
CSMatrix CSMatrix::TranspMtrRect (int *_indarr) const { // Compute transposed to the rectangular matrix stored by columns

	const char *funcname = "TranspMtrRect";

	CSMatrix at(nsupr,nzja);

	int i, j, jj;

	for (i=0;i<=nsupr;i++) at.ia[i] = 0;

	for (i=0;i<nsupc;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			at.ia[jj+1]++;
		};
	};

	for (i=0;i<nsupr;i++) at.ia[i+1] += at.ia[i];

	int *pat;

	pat = new int [nsupr];

	if (!pat) MemoryFail (funcname);

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];
	
	int k;
	for (i=0;i<nsupc;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			_indarr[pat[jj]] = j;
			pat[jj]++;
			k = pat[jj];
			at.ja[k-1] = i;
		};
	};

	for (i=0;i<nsupr;i++) at.list[i] = i;

	delete [] pat;

// Assign control data

	at.m = n;
	at.n = m;
	at.nsupc = nsupr;
	at.nsupr = nsupc;
	at.nlist = nsupr;
	at.nzja = nzja;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed matrix via the mask (list is incomplete)
//========================================================================================
CSMatrix CSMatrix::TranspMtrMask () const { // Compute transposed matrix via the mask (list is incomplete)

	const char *funcname = "TranspMtrMask";

// Init the mask data

	int *imask, *iptr, *listloc, *ialoc;

	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	iptr = new int [nsupr];
	if (!iptr) MemoryFail (funcname);
	listloc = new int [nsupr];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nsupr+1];
	if (!ialoc) MemoryFail (funcname);

	int i, j, jj;

	for (i=0;i<nsupr;i++) imask[i] = 0;

// Count the sparsity of at

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			imask[jj]++;
		};
	};

	int nlistloc = 0;

	int nzloc = 0;

	ialoc[0] = 0;

	for (i=0;i<nsupr;i++) {
		iptr[i] = nzloc;
		if (imask[i] > 0) {
			listloc[nlistloc] = i;
			nzloc += imask[i];
			ialoc[nlistloc+1] = nzloc;
			nlistloc++;
		};
	};

// Init transposed matrix

	CSMatrix at(nlistloc,nzja);

	for (i=0;i<nlistloc;i++) at.list[i] = listloc[i];
	for (i=0;i<=nlistloc;i++) at.ia[i] = ialoc[i];

	int k, irow;
	for (i=0;i<nlist;i++) {
		irow = list[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			iptr[jj]++;
			k = iptr[jj];
			at.ja[k-1] = irow;
		};
	};

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupc = nsupc;
	at.nsupr = nsupr;
	at.nlist = nlistloc;
	at.nzja = nzja;

// Free the memory

	delete [] imask;
	delete [] iptr;
	delete [] listloc;
	delete [] ialoc;

// Return the result

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed matrix via the mask (list is incomplete)
//========================================================================================
CSMatrix CSMatrix::TranspMtrMask (int &_icycle, int *_imask) const { // Compute transposed matrix via the mask (list is incomplete)

	const char *funcname = "TranspMtrMask";

// Init the mask data

	int *imask, *iptr, *listloc, *ialoc;

	imask = new int [nsupr];
	if (!imask) MemoryFail (funcname);
	iptr = new int [nsupr];
	if (!iptr) MemoryFail (funcname);
	listloc = new int [nsupr];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nsupr+1];
	if (!ialoc) MemoryFail (funcname);

	int i, j, jj;

// Count the sparsity of at

	int nlistloc = 0;

	_icycle++;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (_imask[jj] != _icycle) {
				listloc[nlistloc] = jj;
				nlistloc++;
				_imask[jj] = _icycle;
				imask[jj] = 0;
			};
			imask[jj]++;
		};
	};

	sort (listloc, listloc + nlistloc);

//	int nzloc = 0;

	ialoc[0] = 0;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		ialoc[i+1] = ialoc[i]+imask[jj];
		iptr[i] = ialoc[i];
		imask[jj] = i;
	};

// Init transposed matrix

	CSMatrix at(nlistloc,nzja);

	for (i=0;i<nlistloc;i++) at.list[i] = listloc[i];
	for (i=0;i<=nlistloc;i++) at.ia[i] = ialoc[i];

	int k, irow, ind;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			ind = imask[jj];
			iptr[ind]++;
			k = iptr[ind];
			at.ja[k-1] = irow;
		};
	};

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupc = nsupc;
	at.nsupr = nsupr;
	at.nlist = nlistloc;
	at.nzja = nzja;

// Free the memory

	delete [] imask;
	delete [] iptr;
	delete [] listloc;
	delete [] ialoc;

// Return the result

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute transposed 2 index matrix via the mask (list is incomplete)
//========================================================================================
CSMatrix CSMatrix::TranspMtrMask2Index () const { // Compute transposed 2 index matrix via the mask (list is incomplete)

	const char *funcname = "TranspMtrMask2Index";

// Find largest index in matrix and in the list

	int nblks, *blks;

	ComputeBlocks2Index (nblks, blks);

	int *imaskblk;
	int *listblk;

	imaskblk = new int [nblks];
	if (!imaskblk) MemoryFail (funcname);
	listblk = new int [nblks];
	if (!listblk) MemoryFail (funcname);

	int i, j, jj, jjblk;

	for (i=0;i<nblks;i++) imaskblk[i] = -1;

	int nlistblk = 0;

	for (i=0;i<nlist;i++) {
		jjblk = list2[i];
		if (imaskblk[jjblk] == -1) {
			imaskblk[jjblk] = 1;
			listblk[nlistblk] = jjblk;
			nlistblk++;
		};
		for (j=ia[i];j<ia[i+1];j++) {
			jjblk = ja2[j];
			if (imaskblk[jjblk] == -1) {
				imaskblk[jjblk] = 1;
				listblk[nlistblk] = jjblk;
				nlistblk++;
			};
		};
	};

	qsort (listblk, nlistblk, sizeof(int), compint);

// Create block mask array

	int *ibsblks;

	ibsblks = new int [nblks];
	if (!ibsblks) MemoryFail (funcname);

	for (i=0;i<nblks;i++) ibsblks[i] = -1;

	int nz = 0;

	for (i=0;i<nlistblk;i++) {
		jjblk = listblk[i];
		ibsblks[jjblk] = nz;
		nz += blks[jjblk+1]-blks[jjblk];
	};

// Create imask, list and list2 arrays

	int *imask, *listloc, *list2loc, *iptr, *ialoc;

	imask = new int [nz];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nz];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nz];
	if (!list2loc) MemoryFail (funcname);
	iptr = new int [nz];
	if (!iptr) MemoryFail (funcname);
	ialoc = new int [nz+1];
	if (!ialoc) MemoryFail (funcname);

	for (i=0;i<nz;i++) imask[i] = -1;

	int nlistloc = 0;

	int ibs;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jjblk = ja2[j];
			ibs = ibsblks[jjblk]+jj;
			if (imask[ibs] < 0) {
				imask[ibs] = 1;
				listloc[nlistloc] = jj;
				list2loc[nlistloc] = jjblk;
				nlistloc++;
			};
		};
	};

	CInd2Int *i2arr;

	i2arr = new CInd2Int [nz];
	if (!i2arr) MemoryFail (funcname);

	for (i=0;i<nlistloc;i++) {
		i2arr[i].indx = list2loc[i];
		i2arr[i].indy = listloc[i];
	};

	qsort (i2arr, nlistloc, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

	for (i=0;i<nlistloc;i++) {
		list2loc[i] = i2arr[i].indx;
		listloc[i] = i2arr[i].indy;
	};

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		jjblk = list2loc[i];
		ibs = ibsblks[jjblk]+jj;
		imask[ibs] = i;
	};

	delete [] i2arr;

// Count the sparsity of at

	for (i=0;i<=nlistloc;i++) ialoc[i] = 0;

	int ind;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jjblk = ja2[j];
			ibs = ibsblks[jjblk]+jj;
			ind = imask[ibs];
			ialoc[ind+1]++;
		};
	};

	for (i=0;i<nlistloc;i++) ialoc[i+1] += ialoc[i];

	for (i=0;i<nlistloc;i++) iptr[i] = ialoc[i];

	ialoc[0] = 0;

// Init transposed matrix

	CSMatrix at(nlistloc,nlistloc,nzja,nzja);

	for (i=0;i<nlistloc;i++) at.list[i] = listloc[i];
	for (i=0;i<nlistloc;i++) at.list2[i] = list2loc[i];
	for (i=0;i<=nlistloc;i++) at.ia[i] = ialoc[i];

	int k, irow, irow2;
	for (i=0;i<nlist;i++) {
		irow = list[i];
		irow2 = list2[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jjblk = ja2[j];
			ibs = ibsblks[jjblk]+jj;
			ind = imask[ibs];
			iptr[ind]++;
			k = iptr[ind];
			at.ja[k-1] = irow;
			at.ja2[k-1] = irow2;
		};
	};

// Assign control data

	at.m = m;
	at.n = n;
	at.nsupc = nsupc;
	at.nsupr = nsupr;
	at.nlist = nlistloc;
	at.nlist2 = nlistloc;
	at.nzja = nzja;
	at.nzja2 = nzja;

// Free the memory

	delete [] blks;
	delete [] imaskblk;
	delete [] listblk;
	delete [] ibsblks;
	delete [] imask;
	delete [] listloc;
	delete [] list2loc;
	delete [] iptr;
	delete [] ialoc;

// Return the result

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute symmetrized matrix
//========================================================================================
CSMatrix CSMatrix::SymmMtr () const { // Compute symmetrized matrix 

//	const char *funcname = "SymmMtr";

// Compute transposed matrix

	CSMatrix at;

	at = TranspMtr();

// Compute the symmetric matrix

	int nzjas;
	nzjas = 2*nzja+nsupr;

	CSMatrix as(nsupr,nzjas);

	int i, nz;
	int ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;

	as.ia[0] = 0;
	nz = 0;

	for (i=0;i<nsupr;i++) {
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
					as.ja[nz] = ja[jind1];
					nz++;
					jind1++;
				} else if (jj1 > jj2) {
					as.ja[nz] = at.ja[jind2];
					nz++;
					jind2++;
				} else if (jj1 == jj2) {
					as.ja[nz] = ja[jind1];
					nz++;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				as.ja[nz] = ja[jind1];
				nz++;
				jind1++;
			} else if (jind2 <= iend2) {
				as.ja[nz] = at.ja[jind2];
				nz++;
				jind2++;
			};
		};
		as.ia[i+1] = nz;
	};

	for (i=0;i<nsupr;i++) as.list[i] = i;

	as.nzja = as.ia[nsupr];

	return as;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute general symmetrized rectangular matrix
//========================================================================================
CSMatrix CSMatrix::SymmMtrRect () const { // Compute general symmetrized rectangular matrix

	const char *funcname = "SymmMtrRect";

// Check that matrix on entry satisfies assumptions

	int i, j, jj;

	for (i=0;i<nlist-1;i++) {
		if (list[i] > list[i+1]) throw " CSMatrix::SymmMtrRect: nonmonotone list ";
	};

// Find largest index in matrix and in the list

	int nmax = 0;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj > nmax) nmax = jj;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj > nmax) nmax = jj;
		};
	};

	nmax++;

// Create imask and list arrays

	int *imask, *listloc;

	imask = new int [nmax];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nmax];
	if (!listloc) MemoryFail (funcname);

	for (i=0;i<nmax;i++) imask[i] = -1;

	int nlistloc = 0;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (imask[jj] < 0) {
			imask[jj] = 1;
			listloc[nlistloc] = jj;
			nlistloc++;
		};
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (imask[jj] < 0) {
				imask[jj] = 1;
				listloc[nlistloc] = jj;
				nlistloc++;
			};
		};
	};

	qsort (listloc, nlistloc, sizeof(int), compint);

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		imask[jj] = i;
	};

// Compute transposed rectangular matrix

	CSMatrix at = TranspMtrMask ();

// Create ia and ja arrays

	int *ialoc;

	ialoc = new int [nlistloc+1];
	if (!ialoc) MemoryFail (funcname);

	for (i=0;i<=nlistloc;i++) ialoc[i] = 0;

	int ind1 = 0;
	int ind2 = 0;

	int jj1, jj2;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		if (ind1 < nlist && ind2 < at.nlist) {
			jj1 = list[ind1];
			jj2 = at.list[ind2];
			if (jj1 == jj && jj2 == jj) {
				ialoc[i+1] = (ia[ind1+1]-ia[ind1]) + (at.ia[ind2+1]-at.ia[ind2]);
				ind1++;
				ind2++;
			} else if (jj1 == jj) {
				ialoc[i+1] = (ia[ind1+1]-ia[ind1]);
				ind1++;
			} else if (jj2 == jj) {
				ialoc[i+1] = (at.ia[ind2+1]-at.ia[ind2]);
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
			};
		} else if (ind1 < nlist) {
			jj1 = list[ind1];
			if (jj1 == jj) {
				ialoc[i+1] = (ia[ind1+1]-ia[ind1]);
				ind1++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
			};
		} else if (ind2 < at.nlist) {
			jj2 = at.list[ind2];
			if (jj2 == jj) {
				ialoc[i+1] = (at.ia[ind2+1]-at.ia[ind2]);
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
			};
		};
	};

	for (i=0;i<nlistloc;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	int nzjaloc = ialoc[nlistloc];

	int *jaloc;

	jaloc = new int [nzjaloc];
	if (!jaloc) MemoryFail (funcname);

	ind1 = 0;
	ind2 = 0;

	nzjaloc = 0;

	int nz, k, kk;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		if (ind1 < nlist && ind2 < at.nlist) {
			jj1 = list[ind1];
			jj2 = at.list[ind2];
			if (jj1 == jj && jj2 == jj) {
				nz = 0;
				AddLists (ia[ind1+1]-ia[ind1], ja+ia[ind1],
							at.ia[ind2+1]-at.ia[ind2], at.ja+at.ia[ind2],
							nz, jaloc+nzjaloc);
				nzjaloc += nz;
				ind1++;
				ind2++;
			} else if (jj1 == jj) {
				for (k=ia[ind1];k<ia[ind1+1];k++) {
					kk = ja[k];
					jaloc[nzjaloc] = kk;
					nzjaloc++;
				};
				ind1++;
			} else if (jj2 == jj) {
				for (k=at.ia[ind2];k<at.ia[ind2+1];k++) {
					kk = at.ja[k];
					jaloc[nzjaloc] = kk;
					nzjaloc++;
				};
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
			};
		} else if (ind1 < nlist) {
			jj1 = list[ind1];
			if (jj1 == jj) {
				for (k=ia[ind1];k<ia[ind1+1];k++) {
					kk = ja[k];
					jaloc[nzjaloc] = kk;
					nzjaloc++;
				};
				ind1++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
			};
		} else if (ind2 < at.nlist) {
			jj2 = at.list[ind2];
			if (jj2 == jj) {
				for (k=at.ia[ind2];k<at.ia[ind2+1];k++) {
					kk = at.ja[k];
					jaloc[nzjaloc] = kk;
					nzjaloc++;
				};
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
			};
		} else {
			assert (false); throw " CSMatrix::SymmMtrRect: index not found ";
		};
		ialoc[i+1] = nzjaloc;
	};

// Store the result

// Init transposed matrix

	CSMatrix as(nlistloc,nzjaloc);

	for (i=0;i<nlistloc;i++) as.list[i] = listloc[i];
	for (i=0;i<=nlistloc;i++) as.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) as.ja[i] = jaloc[i];

	as.m = m;
	as.n = n;
	as.nsupc = nsupc;
	as.nsupr = nsupr;
	as.nlist = nlistloc;
	as.nzja = nzjaloc;

// Free the memory

	delete [] imask;
	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;

// Return the result

	return as;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute general symmetrized rectangular matrix
//========================================================================================
CSMatrix CSMatrix::SymmMtrRect2Index () const { // Compute general symmetrized rectangular matrix

	const char *funcname = "SymmMtrRect2Index";

// Check that matrix on entry satisfies assumptions

	int i, j, jj, jjblk;

	for (i=0;i<nlist-1;i++) {
		if (list2[i] > list2[i+1]) {
			throw " CSMatrix::SymmMtrRect2Index: nonmonotone list ";
		};
		if (list2[i] == list2[i+1]) {
			if (list[i] > list[i+1]) {
				throw " CSMatrix::SymmMtrRect2Index: nonmonotone list ";
			};
		};
	};

// Find largest index in matrix and in the list

	int nblks, *blks;

	ComputeBlocks2Index (nblks, blks);

	int *imaskblk;
	int *listblk;

	imaskblk = new int [nblks];
	if (!imaskblk) MemoryFail (funcname);
	listblk = new int [nblks];
	if (!listblk) MemoryFail (funcname);

	for (i=0;i<nblks;i++) imaskblk[i] = -1;

	int nlistblk = 0;

	for (i=0;i<nlist;i++) {
		jjblk = list2[i];
		if (imaskblk[jjblk] == -1) {
			imaskblk[jjblk] = 1;
			listblk[nlistblk] = jjblk;
			nlistblk++;
		};
		for (j=ia[i];j<ia[i+1];j++) {
			jjblk = ja2[j];
			if (imaskblk[jjblk] == -1) {
				imaskblk[jjblk] = 1;
				listblk[nlistblk] = jjblk;
				nlistblk++;
			};
		};
	};

	qsort (listblk, nlistblk, sizeof(int), compint);

// Create block mask array

	int *ibsblks;

	ibsblks = new int [nblks];
	if (!ibsblks) MemoryFail (funcname);

	for (i=0;i<nblks;i++) ibsblks[i] = -1;

	int nz = 0;

	for (i=0;i<nlistblk;i++) {
		jjblk = listblk[i];
		ibsblks[jjblk] = nz;
		nz += blks[jjblk+1]-blks[jjblk];
	};

// Create imask, list and list2 arrays

	int *imask, *listloc, *list2loc;

	imask = new int [nz];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nz];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nz];
	if (!list2loc) MemoryFail (funcname);

	for (i=0;i<nz;i++) imask[i] = -1;

	int nlistloc = 0;

	int ibs;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		jjblk = list2[i];
		ibs = ibsblks[jjblk]+jj;
		if (imask[ibs] < 0) {
			imask[ibs] = 1;
			listloc[nlistloc] = jj;
			list2loc[nlistloc] = jjblk;
			nlistloc++;
		};
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jjblk = ja2[j];
			ibs = ibsblks[jjblk]+jj;
			if (imask[ibs] < 0) {
				imask[ibs] = 1;
				listloc[nlistloc] = jj;
				list2loc[nlistloc] = jjblk;
				nlistloc++;
			};
		};
	};

	CInd2Int *i2arr;

	i2arr = new CInd2Int [nz];
	if (!i2arr) MemoryFail (funcname);

	for (i=0;i<nlistloc;i++) {
		i2arr[i].indx = list2loc[i];
		i2arr[i].indy = listloc[i];
	};

	qsort (i2arr, nlistloc, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

	for (i=0;i<nlistloc;i++) {
		list2loc[i] = i2arr[i].indx;
		listloc[i] = i2arr[i].indy;
	};

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		jjblk = list2loc[i];
		ibs = ibsblks[jjblk]+jj;
		imask[ibs] = i;
	};

	delete [] i2arr;

// Compute transposed rectangular matrix

	CSMatrix at = TranspMtrMask2Index ();

// Create ia and ja arrays

	int *ialoc;

	ialoc = new int [nlistloc+1];
	if (!ialoc) MemoryFail (funcname);

	for (i=0;i<=nlistloc;i++) ialoc[i] = 0;

	int ind1 = 0;
	int ind2 = 0;

	int jj1, jj2, jjblk1, jjblk2;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		jjblk = list2loc[i];
		if (ind1 < nlist && ind2 < at.nlist) {
			jj1 = list[ind1];
			jjblk1 = list2[ind1];
			jj2 = at.list[ind2];
			jjblk2 = at.list2[ind2];
			if (jjblk1 == jjblk && jjblk2 == jjblk && jj1 == jj && jj2 == jj) {
				ialoc[i+1] = (ia[ind1+1]-ia[ind1]) + (at.ia[ind2+1]-at.ia[ind2]);
				ind1++;
				ind2++;
			} else if (jjblk1 == jjblk && jj1 == jj) {
				ialoc[i+1] = (ia[ind1+1]-ia[ind1]);
				ind1++;
			} else if (jjblk2 == jjblk && jj2 == jj) {
				ialoc[i+1] = (at.ia[ind2+1]-at.ia[ind2]);
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
			};
		} else if (ind1 < nlist) {
			jj1 = list[ind1];
			jjblk1 = list2[ind1];
			if (jjblk1 == jjblk && jj1 == jj) {
				ialoc[i+1] = (ia[ind1+1]-ia[ind1]);
				ind1++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
			};
		} else if (ind2 < at.nlist) {
			jj2 = at.list[ind2];
			jjblk2 = at.list2[ind2];
			if (jjblk2 == jjblk && jj2 == jj) {
				ialoc[i+1] = (at.ia[ind2+1]-at.ia[ind2]);
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
			};
		};
	};

	for (i=0;i<nlistloc;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	int nzjaloc = ialoc[nlistloc];

	int *jaloc, *ja2loc;

	jaloc = new int [nzjaloc];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjaloc];
	if (!ja2loc) MemoryFail (funcname);

	ind1 = 0;
	ind2 = 0;

	nzjaloc = 0;

	int k;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		jjblk = list2loc[i];
		if (ind1 < nlist && ind2 < at.nlist) {
			jj1 = list[ind1];
			jjblk1 = list2[ind1];
			jj2 = at.list[ind2];
			jjblk2 = at.list2[ind2];
			if (jjblk1 == jjblk && jjblk2 == jjblk && jj1 == jj && jj2 == jj) {
				nz = 0;
				AddLists (ia[ind1+1]-ia[ind1], ja+ia[ind1], ja2+ia[ind1],
							at.ia[ind2+1]-at.ia[ind2], at.ja+at.ia[ind2], at.ja2+at.ia[ind2],
							nz, jaloc+nzjaloc, ja2loc+nzjaloc);
				nzjaloc += nz;
				ind1++;
				ind2++;
			} else if (jjblk1 == jjblk && jj1 == jj) {
				for (k=ia[ind1];k<ia[ind1+1];k++) {
					jaloc[nzjaloc] = ja[k];
					ja2loc[nzjaloc] = ja2[k];
					nzjaloc++;
				};
				ind1++;
			} else if (jjblk2 == jjblk && jj2 == jj) {
				for (k=at.ia[ind2];k<at.ia[ind2+1];k++) {
					jaloc[nzjaloc] = at.ja[k];
					ja2loc[nzjaloc] = at.ja2[k];
					nzjaloc++;
				};
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
			};
		} else if (ind1 < nlist) {
			jj1 = list[ind1];
			jjblk1 = list2[ind1];
			if (jjblk1 == jjblk && jj1 == jj) {
				for (k=ia[ind1];k<ia[ind1+1];k++) {
					jaloc[nzjaloc] = ja[k];
					ja2loc[nzjaloc] = ja2[k];
					nzjaloc++;
				};
				ind1++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
			};
		} else if (ind2 < at.nlist) {
			jj2 = at.list[ind2];
			jjblk2 = at.list2[ind2];
			if (jjblk2 == jjblk && jj2 == jj) {
				for (k=at.ia[ind2];k<at.ia[ind2+1];k++) {
					jaloc[nzjaloc] = at.ja[k];
					ja2loc[nzjaloc] = at.ja2[k];
					nzjaloc++;
				};
				ind2++;
			} else {
				assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
			};
		} else {
			assert (false); throw " CSMatrix::SymmMtrRect2Index: index not found ";
		};
		ialoc[i+1] = nzjaloc;
	};

// Store the result

// Init transposed matrix

	CSMatrix as(nlistloc,nlistloc,nzjaloc,nzjaloc);

	for (i=0;i<nlistloc;i++) as.list[i] = listloc[i];
	for (i=0;i<nlistloc;i++) as.list2[i] = list2loc[i];
	for (i=0;i<=nlistloc;i++) as.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) as.ja[i] = jaloc[i];
	for (i=0;i<nzjaloc;i++) as.ja2[i] = ja2loc[i];

	as.m = m;
	as.n = n;
	as.nsupc = nsupc;
	as.nsupr = nsupr;
	as.nlist = nlistloc;
	as.nlist2 = nlistloc;
	as.nzja = nzjaloc;
	as.nzja2 = nzjaloc;

// Free the memory

	delete [] blks;
	delete [] imaskblk;
	delete [] listblk;
	delete [] ibsblks;
	delete [] imask;
	delete [] listloc;
	delete [] list2loc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] ja2loc;

// Return the result

	return as;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute minimal blocks partitioning for 2 index matrix
//========================================================================================
void CSMatrix::ComputeBlocks2Index (int &_nblks, int *&_blks) const { // Compute minimal blocks partitioning for 2 index matrix

	const char *funcname = "ComputeBlocks2Index";

	int nblkmax = 0;

	int i, jj, jjblk, j;

	for (i=0;i<nlist;i++) {
		jjblk = list2[i];
		if (jjblk > nblkmax) nblkmax = jjblk;
		for (j=ia[i];j<ia[i+1];j++) {
			jjblk = ja2[j];
			if (jjblk > nblkmax) nblkmax = jjblk;
		};
	};

	nblkmax++;

	_nblks = nblkmax;

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);

	for (i=0;i<=_nblks;i++) _blks[i] = 0;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		jjblk = list2[i];
		if (_blks[jjblk+1] < jj) _blks[jjblk+1] = jj;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jjblk = ja2[j];
			if (_blks[jjblk+1] < jj) _blks[jjblk+1] = jj;
		};
	};

	for (i=0;i<_nblks;i++) _blks[i+1]++;

	for (i=0;i<_nblks;i++) _blks[i+1] += _blks[i];

};

// Author: Kharchenko S.A.
// CSMatrix: Transform 2Index matrix sparsity into one index
//========================================================================================
void CSMatrix::Transform2IndexIntoContinuous (int *_blks) { // Transform 2Index matrix sparsity into one index

	int i, j, jj, jjblk, jjnew;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		jjblk = list2[i];
		jjnew = _blks[jjblk]+jj;
		list[i] = jjnew;
		list2[i] = -1;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jjblk = ja2[j];
			jjnew = _blks[jjblk]+jj;
			ja[j] = jjnew;
			ja2[j] = -1;
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrix: Transform 2Index matrix sparsity into one index
//========================================================================================
void CSMatrix::Transform2IndexIntoContinuous (CMPIComm &_comm) { // Transform 2Index matrix sparsity into one index

	const char *funcname = "Transform2IndexIntoContinuous";

	int myid = _comm.GetMyid ();
	int nproc = _comm.GetNproc ();

// Determine the maximal sp index number

	int indspmax = 0;

	int i, isup;

	for (i=0;i<nlist;i++) {
		isup = list2[i];
		if (isup > indspmax) indspmax = isup;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, MAXIMUM, 1, &indspmax, &indspmax);

	int nsuptot = indspmax+1;

// Compute local partitioning

	int nsuploc = 0;
	int isupprev = -1;

	if (nlist > 0) {
		nsuploc = 1;
		isupprev = list2[0];
	};

	for (i=1;i<nlist;i++) {
		isup = list2[i];
		if (isup != isupprev) {
			isupprev = isup;
			nsuploc++;
		};
	};

	int *listsp;
	int *sprndsloc;

	listsp = new int [nsuploc];
	if (!listsp) MemoryFail (funcname);
	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	nsuploc = 0;
	isupprev = -1;

	sprndsloc[0] = 0;

	if (nlist > 0) {
		isupprev = list2[0];
		listsp[nsuploc] = isupprev;
		nsuploc = 1;
		sprndsloc[nsuploc] = 1;
	};

	for (i=1;i<nlist;i++) {
		sprndsloc[nsuploc] = i+1;
		isup = list2[i];
		if (isup != isupprev) {
			sprndsloc[nsuploc] = i;
			isupprev = isup;
			listsp[nsuploc] = isupprev;
			nsuploc++;
		};
	};

	sprndsloc[nsuploc] = nlist;

// Check that data inside supernodes are also contiguous

	int j;

	for (i=0;i<nsuploc;i++) {
		for (j=sprndsloc[i];j<sprndsloc[i+1];j++) {
			if (list2[j] != listsp[i] || list[j] != j-sprndsloc[i]) {
				throw " CSMatrix::Transform2IndexIntoContinuous: error in input list data ";
			};
		};
	};

// Compute local contiguous blocks

	int *blks2sploc;

	blks2sploc = new int [nsuploc+1];
	if (!blks2sploc) MemoryFail (funcname);

	int nblksloc = 0;
	isupprev = -1;

	blks2sploc[0] = 0;

	if (nsuploc > 0) {
		isupprev = listsp[0];
		nblksloc = 1;
		blks2sploc[nblksloc] = 1;
	};

	for (i=1;i<nsuploc;i++) {
		blks2sploc[nblksloc] = i+1;
		isup = listsp[i];
		if (isup != isupprev+1) {
			blks2sploc[nblksloc] = i;
			nblksloc++;
		};
		isupprev = isup;
	};

	blks2sploc[nblksloc] = nsuploc;

	int *blksloc;
	int *ibegloc, *iendloc;

	blksloc = new int [nblksloc+1];
	if (!blksloc) MemoryFail (funcname);
	ibegloc = new int [nblksloc];
	if (!ibegloc) MemoryFail (funcname);
	iendloc = new int [nblksloc];
	if (!iendloc) MemoryFail (funcname);

	blksloc[0] = 0;

	int ispbeg, ispend, niloc;

	for (i=0;i<nblksloc;i++) {
		ispbeg = blks2sploc[i];
		ispend = blks2sploc[i+1]-1;
		niloc = sprndsloc[ispend+1]-sprndsloc[ispbeg];
		blksloc[i+1] = blksloc[i] + niloc;
		ibegloc[i] = listsp[ispbeg];
		iendloc[i] = listsp[ispend];
	};

// Count the global number of blocks

	int *iabl2cpu;

	iabl2cpu = new int [nproc+1];
	if (!iabl2cpu) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iabl2cpu[i] = 0;

	iabl2cpu[myid+1] = nblksloc;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, iabl2cpu, iabl2cpu);

	for (i=0;i<nproc;i++) iabl2cpu[i+1] = iabl2cpu[i]+iabl2cpu[i+1];

	int nblkstot = iabl2cpu[nproc];

// Prepare global arrays

	int *ibegtot, *iendtot, *niblktot, *icputot;

	ibegtot = new int [nblkstot];
	if (!ibegtot) MemoryFail (funcname);
	iendtot = new int [nblkstot];
	if (!iendtot) MemoryFail (funcname);
	niblktot = new int [nblkstot];
	if (!niblktot) MemoryFail (funcname);
	icputot = new int [nblkstot];
	if (!icputot) MemoryFail (funcname);

	for (i=0;i<nblkstot;i++) ibegtot[i] = 0;
	for (i=0;i<nblkstot;i++) iendtot[i] = 0;
	for (i=0;i<nblkstot;i++) niblktot[i] = 0;

	for (j=iabl2cpu[myid];j<iabl2cpu[myid+1];j++) ibegtot[j] = ibegloc[j-iabl2cpu[myid]];
	for (j=iabl2cpu[myid];j<iabl2cpu[myid+1];j++) iendtot[j] = iendloc[j-iabl2cpu[myid]];
	for (j=iabl2cpu[myid];j<iabl2cpu[myid+1];j++) niblktot[j] = blksloc[j-iabl2cpu[myid]+1]-blksloc[j-iabl2cpu[myid]];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nblkstot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nblkstot, iendtot, iendtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nblkstot, niblktot, niblktot);

	for (i=0;i<nproc;i++) {
		for (j=iabl2cpu[i];j<iabl2cpu[i+1];j++) {
			icputot[j] = i;
		};
	};

// Sort global data

	CIntInt *iiarr;

	iiarr = new CIntInt [nblkstot];
	if (!iiarr) MemoryFail (funcname);

	for (i=0;i<nblkstot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nblkstot, sizeof(CIntInt), CIntInt::CompareIntInt);

	int *order;

	order = new int [nblkstot];
	if (!order) MemoryFail (funcname);

	for (i=0;i<nblkstot;i++) {
		j = iiarr[i].int2value;
		order[j] = i;
	};

	int *imask;

	imask = new int [nblkstot];
	if (!imask) MemoryFail (funcname);

	for (i=0;i<nblkstot;i++) {
		j = order[i];
		imask[j] = ibegtot[i];
	};
	for (i=0;i<nblkstot;i++) ibegtot[i] = imask[i];

	for (i=0;i<nblkstot;i++) {
		j = order[i];
		imask[j] = iendtot[i];
	};
	for (i=0;i<nblkstot;i++) iendtot[i] = imask[i];

	for (i=0;i<nblkstot;i++) {
		j = order[i];
		imask[j] = niblktot[i];
	};
	for (i=0;i<nblkstot;i++) niblktot[i] = imask[i];

	for (i=0;i<nblkstot;i++) {
		j = order[i];
		imask[j] = icputot[i];
	};
	for (i=0;i<nblkstot;i++) icputot[i] = imask[i];

// Create global blocks partitioning

	int *blks2sptot;
	int *blkstot;

	blks2sptot = new int [nblkstot+1];
	if (!blks2sptot) MemoryFail (funcname);
	blkstot = new int [nblkstot+1];
	if (!blkstot) MemoryFail (funcname);

	blks2sptot[0] = 0;
	for (i=0;i<nblkstot;i++) blks2sptot[i+1] = blks2sptot[i] + (iendtot[i]-ibegtot[i]+1);

	if (blks2sptot[nblkstot] != nsuptot) throw " CSMatrix::Transform2IndexIntoContinuous: error in input data ";

	blkstot[0] = 0;
	for (i=0;i<nblkstot;i++) blkstot[i+1] = blkstot[i] + niblktot[i];

	int ntot = blkstot[nblkstot];

// Find global to local mapping

	int *blkl2blkg;
	int *blkg2blkl;

	blkl2blkg = new int [nblksloc];
	if (!blkl2blkg) MemoryFail (funcname);
	blkg2blkl = new int [nblkstot];
	if (!blkg2blkl) MemoryFail (funcname);

	int iblkloc = 0;

	int iblk;

	for (iblk=0;iblk<nblkstot;iblk++) {
		blkg2blkl[iblk] = -1;
		if (icputot[iblk] == myid) {
			blkl2blkg[iblkloc] = iblk;
			blkg2blkl[iblk] = iblkloc;
			iblkloc++;
		};
	};

// Compute shifts for all local sups

	int *ibssploc;

	ibssploc = new int [nsuploc];
	if (!ibssploc) MemoryFail (funcname);

	int isuploc = 0;

	int ibeg;

	for (iblkloc=0;iblkloc<nblksloc;iblkloc++) {
		iblk = blkl2blkg[iblkloc];
		ibeg = sprndsloc[isuploc];
		for (i=blks2sptot[iblk];i<blks2sptot[iblk+1];i++) {
			ibssploc[isuploc] = blkstot[iblk] + (sprndsloc[isuploc]-ibeg);
			isuploc++;
		};
	};

// Reassign the list values

	int ibs;

	for (i=0;i<nsuploc;i++) {
		ibs = ibssploc[i];
		for (j=sprndsloc[i];j<sprndsloc[i+1];j++) {
			list[j] = ibs+j-sprndsloc[i];
		};
	};

	delete [] list2;

	list2 = new int [0];
	if (!list2) MemoryFail (funcname);

	nlist2 = 0;

// Compute block numbers for the column indices

	int *ja2blk;

	ja2blk = new int [nzja];
	if (!ja2blk) MemoryFail (funcname);

	int iblkprev = 0;

	int jjsup, jjblk, iblkbeg, iblkend;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jjsup = ja2[j];
			if (jjsup >= blks2sptot[iblkprev] && jjsup < blks2sptot[iblkprev+1]) {
				jjblk = iblkprev;
			} else {
				if (j==ia[i]) {
					iblkbeg = 0;
					iblkend = nblkstot-1;
				} else {
					iblkbeg = iblkprev+1;
					iblkend = nblkstot-1;
				};
				while (iblkbeg != iblkend) {
					iblk = (iblkbeg+iblkend)/2;
					if (iblk < iblkbeg) iblk = iblkbeg;
					if (iblk > iblkend) iblk = iblkend;
					if (jjsup >= blks2sptot[iblk] && jjsup < blks2sptot[iblk+1]) {
						iblkbeg = iblk;
						iblkend = iblk;
						break;
					} else {
						if (jjsup < blks2sptot[iblk]) iblkend = iblk-1;
						if (jjsup >= blks2sptot[iblk+1]) iblkbeg = iblk+1;
					};
				};
				jjblk = iblkbeg;
			};
			ja2blk[j] = jjblk;
			iblkprev = jjblk;
		};
	};

// Compute cpu partitioning of the column indices

	int *iacpu;
	int *iptrcpu;
	int *jaindcpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);
	jaindcpu = new int [nzja];
	if (!jaindcpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int iproc;

	for (j=0;j<nzja;j++) {
		jjblk = ja2blk[j];
		iproc = icputot[jjblk];
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	int k;

	for (j=0;j<nzja;j++) {
		jjblk = ja2blk[j];
		iproc = icputot[jjblk];
		k = iptrcpu[iproc];
		jaindcpu[k] = j;
		iptrcpu[iproc]++;
	};

// Create lists of send/receive data

	int **plists;

	plists = new int * [nproc];
	if (!plists) MemoryFail (funcname);

	int ni, ind;

	for (i=0;i<nproc;i++) {
		ni = iacpu[i+1]-iacpu[i];
		plists[i] = new int [2*ni];
		if (!plists) MemoryFail (funcname);
		for (j=iacpu[i];j<iacpu[i+1];j++) {
			ind = jaindcpu[j];
			iblk = ja2blk[ind];
			isup = ja2[ind];
			plists[i][(j-iacpu[i])*2] = iblk;
			plists[i][(j-iacpu[i])*2+1] = isup;
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

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 2*(iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) plists[i];
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
	if(info) throw " CSMatrix::Transform2IndexIntoContinuous: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Renumber the indices

	int *piarr;

	int isup0, iblkl, isuploc0;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			iblk = piarr[j*2];
			isup = piarr[j*2+1];
			isup0 = blks2sptot[iblk];
			iblkl = blkg2blkl[iblk];
			if (iblkl < 0) throw " CSMatrix::Transform2IndexIntoContinuous: Block not found";
			isuploc0 = blks2sploc[iblkl];
			isuploc = isuploc0 + (isup-isup0);
			ibs = ibssploc[isuploc];
			piarr[j*2] = ibs;
		};
	};

// Exchange the data once again

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CDecomp::TransformSparsity: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Use received data

	int jloc, jj;

	for (i=0;i<NObjSend;i++) {
		iproc = CpuIDSend[i];
		ni = ObjSizeSend[i] / (2*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ibs = piarr[j*2];
			jloc = iacpu[iproc]+j;
			ind = jaindcpu[jloc];
			jj = ja[ind];
			ja[ind] = ibs + jj;
		};
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Replace ja2 data

	delete [] ja2;

	ja2 = new int [0];
	if(!ja2) MemoryFail(funcname);

	nzja2 = 0;

	m = ntot;
	n = ntot;
	nsupc = nlist;
	nsupr = nlist;

// Free work arrays

	delete [] listsp;
	delete [] sprndsloc;
	delete [] blks2sploc;
	delete [] blksloc;
	delete [] ibegloc;
	delete [] iendloc;
	delete [] iabl2cpu;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] niblktot;
	delete [] icputot;
	delete [] iiarr;
	delete [] order;
	delete [] imask;
	delete [] blks2sptot;
	delete [] blkstot;
	delete [] blkl2blkg;
	delete [] blkg2blkl;
	delete [] ibssploc;
	delete [] ja2blk;
	delete [] iacpu;
	delete [] iptrcpu;
	delete [] jaindcpu;

	for (i=0;i<nproc;i++) delete [] plists[i];

	delete [] plists;

};

// Author: Kharchenko S.A.
// CSMatrix: Combine sparsity structures of the two matrices
//========================================================================================
CSMatrix CSMatrix::CombineStructures (const CSMatrix &_a1, const CSMatrix &_a2) { // Combine sparsity structures of the two matrices

//	const char *funcname = "CombineStructures";

// Allocate combined matrix

	int nlistloc = _a1.nlist + _a2.nlist;
	int nzjaloc = _a1.nzja + _a2.nzja;

	CSMatrix aloc (nlistloc,nzjaloc);

// Compute combined matrix

	int ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;
	int ilist1, ilist2, isup1, isup2, j;

	aloc.ia[0] = 0;

	nlistloc = 0;
	nzjaloc = 0;

	ilist1 = 0;
	ilist2 = 0;

	while (ilist1 < _a1.nlist || ilist2 < _a2.nlist) {
		if (ilist1 < _a1.nlist && ilist2 < _a2.nlist) {
			isup1 = _a1.list[ilist1];
			isup2 = _a2.list[ilist2];
			if (isup1 == isup2) {
				ibeg1 = _a1.ia[ilist1]; 
				iend1 = _a1.ia[ilist1+1]-1; 
				ibeg2 = _a2.ia[ilist2]; 
				iend2 = _a2.ia[ilist2+1]-1; 
				jind1 = ibeg1;
				jind2 = ibeg2;
				while (jind1 <= iend1 || jind2 <= iend2) {
					if (jind1 <= iend1 && jind2 <= iend2) {
						jj1 = _a1.ja[jind1];
						jj2 = _a2.ja[jind2];
						if (jj1 < jj2) {
							aloc.ja[nzjaloc] = _a1.ja[jind1];
							nzjaloc++;
							jind1++;
						} else if (jj1 > jj2) {
							aloc.ja[nzjaloc] = _a2.ja[jind2];
							nzjaloc++;
							jind2++;
						} else if (jj1 == jj2) {
							aloc.ja[nzjaloc] = _a1.ja[jind1];
							nzjaloc++;
							jind1++;
							jind2++;
						};
					} else if (jind1 <= iend1) {
						aloc.ja[nzjaloc] = _a1.ja[jind1];
						nzjaloc++;
						jind1++;
					} else if (jind2 <= iend2) {
						aloc.ja[nzjaloc] = _a2.ja[jind2];
						nzjaloc++;
						jind2++;
					};
				};
				aloc.list[nlistloc] = _a1.list[ilist1];
				ilist1++;
				ilist2++;
			} else if (isup1 < isup2) {
				for (j=_a1.ia[ilist1];j<_a1.ia[ilist1+1];j++) {
					aloc.ja[nzjaloc] = _a1.ja[j];
					nzjaloc++;
				};
				aloc.list[nlistloc] = _a1.list[ilist1];
				ilist1++;
			} else if (isup1 > isup2) {
				for (j=_a2.ia[ilist2];j<_a2.ia[ilist2+1];j++) {
					aloc.ja[nzjaloc] = _a2.ja[j];
					nzjaloc++;
				};
				aloc.list[nlistloc] = _a2.list[ilist2];
				ilist2++;
			};
		} else if (ilist1 < _a1.nlist) {
			for (j=_a1.ia[ilist1];j<_a1.ia[ilist1+1];j++) {
				aloc.ja[nzjaloc] = _a1.ja[j];
				nzjaloc++;
			};
			aloc.list[nlistloc] = _a1.list[ilist1];
			ilist1++;
		} else if (ilist2 < _a2.nlist) {
			for (j=_a2.ia[ilist2];j<_a2.ia[ilist2+1];j++) {
				aloc.ja[nzjaloc] = _a2.ja[j];
				nzjaloc++;
			};
			aloc.list[nlistloc] = _a2.list[ilist2];
			ilist2++;
		};
		aloc.ia[nlistloc+1] = nzjaloc;
		nlistloc++;
	};

// Register the head

	aloc.m = _a1.m;
	aloc.n = _a1.n;
	aloc.nsupr = _a1.nsupr;
	aloc.nsupc = _a1.nsupc;
	aloc.nlist = nlistloc;
	aloc.nzja = nzjaloc;

// Return the result

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute block sparsity of a sparse matrix (incomplete lists version)
//========================================================================================
CSMatrix CSMatrix::BlockSparsity (int _nblks, int *_sp2blk) const { // Compute block sparsity of a sparse matrix (incomplete lists version)

	const char *funcname = "BlockSparsity";

// Allocate the data

	int *imask, *iptr;
	int *ialist, *jalist, *listloc;

	imask = new int [_nblks];
	if (!imask) MemoryFail (funcname);
	iptr = new int [_nblks];
	if (!iptr) MemoryFail (funcname);
	ialist = new int [_nblks+1];
	if (!ialist) MemoryFail (funcname);
	jalist = new int [nsupr];
	if (!jalist) MemoryFail (funcname);
	listloc = new int [_nblks];
	if (!listloc) MemoryFail (funcname);

	int i, j, jj, k, iblk, ilist, jblk;

	int icycle = -1;

	for (i=0;i<_nblks;i++) imask[i] = icycle;
	for (i=0;i<=_nblks;i++) ialist[i] = 0;

// Create the lists of supernodes to be scanned per each block

	icycle++;

	ialist[0] = 0;

	for (i=0;i<nlist;i++) {
		j = list[i];
		jj = _sp2blk[j];
		ialist[jj+1]++;
	};

	for (i=0;i<_nblks;i++) ialist[i+1] = ialist[i] + ialist[i+1];

	for (i=0;i<_nblks;i++) iptr[i] = ialist[i];

	for (i=0;i<nlist;i++) {
		j = list[i];
		jj = _sp2blk[j];
		k = iptr[jj];
		jalist[k] = i;
		iptr[jj]++;
	};

// Count the total number of nonzero blocks

	int nlistblk = 0;
	int nzblk = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (ialist[iblk+1] > ialist[iblk]) nlistblk++;
		icycle++;
		for (i=ialist[iblk];i<ialist[iblk+1];i++) {
			ilist = jalist[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jblk = _sp2blk[jj];
				if (imask[jblk] != icycle) {
					nzblk++;
					imask[jblk] = icycle;
				};
			};
		};
	};

// Allocate and init the block matrix

	CSMatrix aloc (nlistblk,nzblk);

	nlistblk = 0;
	nzblk = 0;

	aloc.ia[0] = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		icycle++;
		int nlistloc = 0;
		for (i=ialist[iblk];i<ialist[iblk+1];i++) {
			ilist = jalist[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jblk = _sp2blk[jj];
				if (imask[jblk] != icycle) {
					listloc[nlistloc] = jblk;
					nlistloc++;
					imask[jblk] = icycle;
				};
			};
		};
		if (ialist[iblk+1] > ialist[iblk]) {
			if (nlistloc != 0) qsort (listloc, nlistloc, sizeof(int), compint);
			for (j=0;j<nlistloc;j++) aloc.ja[nzblk+j] = listloc[j];
			nzblk += nlistloc;
			aloc.list[nlistblk] = iblk;
			aloc.ia[nlistblk+1] = nzblk;
			nlistblk++;
		};
	};

// Register the head

	aloc.m = nsupr;
	aloc.n = nsupr;
	aloc.nsupr = _nblks;
	aloc.nsupc = _nblks;
	aloc.nlist = nlistblk;
	aloc.nzja = nzblk;

// Free the memory

	delete [] imask;
	delete [] iptr;
	delete [] ialist;
	delete [] jalist;
	delete [] listloc;

// Return the result

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute block sparsity of a sparse matrix (incomplete lists version) with added diagonal
//========================================================================================
CSMatrix CSMatrix::BlockSparsityDiag (int _nblks, int *_sp2blk) const { // Compute block sparsity of a sparse matrix (incomplete lists version) with added diagonal

	const char *funcname = "BlockSparsityDiag";

// Allocate the data

	int *imask, *iptr;
	int *ialist, *jalist, *listloc;

	imask = new int [_nblks];
	if (!imask) MemoryFail (funcname);
	iptr = new int [_nblks];
	if (!iptr) MemoryFail (funcname);
	ialist = new int [_nblks+1];
	if (!ialist) MemoryFail (funcname);
	jalist = new int [nsupr];
	if (!jalist) MemoryFail (funcname);
	listloc = new int [_nblks];
	if (!listloc) MemoryFail (funcname);

	int i, j, jj, k, iblk, ilist, jblk;

	int icycle = -1;

	for (i=0;i<_nblks;i++) imask[i] = icycle;
	for (i=0;i<=_nblks;i++) ialist[i] = 0;

// Create the lists of supernodes to be scanned per each block

	icycle++;

	ialist[0] = 0;

	for (i=0;i<nlist;i++) {
		j = list[i];
		jj = _sp2blk[j];
		ialist[jj+1]++;
	};

	for (i=0;i<_nblks;i++) ialist[i+1] = ialist[i] + ialist[i+1];

	for (i=0;i<_nblks;i++) iptr[i] = ialist[i];

	for (i=0;i<nlist;i++) {
		j = list[i];
		jj = _sp2blk[j];
		k = iptr[jj];
		jalist[k] = i;
		iptr[jj]++;
	};

// Count the total number of nonzero blocks

	int nlistblk = 0;
	int nzblk = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
//		if (ialist[iblk+1] > ialist[iblk]) nlistblk++;
		nlistblk++;
		icycle++;
		for (i=ialist[iblk];i<ialist[iblk+1];i++) {
			ilist = jalist[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jblk = _sp2blk[jj];
				if (imask[jblk] != icycle) {
					nzblk++;
					imask[jblk] = icycle;
				};
			};
		};
		if (ialist[iblk+1] == ialist[iblk]) nzblk++;
	};

// Allocate and init the block matrix

	CSMatrix aloc (nlistblk,nzblk);

	nlistblk = 0;
	nzblk = 0;

	aloc.ia[0] = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		icycle++;
		int nlistloc = 0;
		for (i=ialist[iblk];i<ialist[iblk+1];i++) {
			ilist = jalist[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jblk = _sp2blk[jj];
				if (imask[jblk] != icycle) {
					listloc[nlistloc] = jblk;
					nlistloc++;
					imask[jblk] = icycle;
				};
			};
		};
		if (ialist[iblk+1] > ialist[iblk]) {
			if (nlistloc != 0) qsort (listloc, nlistloc, sizeof(int), compint);
			for (j=0;j<nlistloc;j++) aloc.ja[nzblk+j] = listloc[j];
			nzblk += nlistloc;
			aloc.list[nlistblk] = iblk;
			aloc.ia[nlistblk+1] = nzblk;
			nlistblk++;
		} else {
			aloc.list[nlistblk] = iblk;
			aloc.ja[nzblk] = iblk;
			nzblk++;
			aloc.ia[nlistblk+1] = nzblk;
			nlistblk++;
		};
	};

// Register the head

	aloc.m = nsupr;
	aloc.n = nsupr;
	aloc.nsupr = _nblks;
	aloc.nsupc = _nblks;
	aloc.nlist = nlistblk;
	aloc.nzja = nzblk;

// Free the memory

	delete [] imask;
	delete [] iptr;
	delete [] ialist;
	delete [] jalist;
	delete [] listloc;

// Return the result

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute optimal boundary
//========================================================================================
void CSMatrix::OptimalBoundary (int *_imask, // Compute optimal boundary
											int &_nlistbnd, int *_listbnd) const {

	const char *funcname = "OptimalBoundary";

// Determine the two-sided boundaries

	int *imasknd;
	int *listloc;

	imasknd = new int [n];
	if(!imasknd) MemoryFail(funcname);
	listloc = new int [n];
	if(!listloc) MemoryFail(funcname);

	int i, ipart;

	for (i=0;i<n;i++) {
		ipart = _imask[i];
		imasknd[i] = ipart;
	};

	int nlistloc = 0;

	int j, jj, jpart;

	for (i=0;i<n;i++) {
		ipart = _imask[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jpart = _imask[jj];
			if (jpart != ipart) {
				if (imasknd[i] >= 0) {
					imasknd[i] = -3;
					listloc[nlistloc] = i;
					nlistloc++;
				};
			};
		};
	};

// Split the list into two separate lists

	int *listloc1;
	int *listloc2;

	listloc1 = new int [n];
	if(!listloc1) MemoryFail(funcname);
	listloc2 = new int [n];
	if(!listloc2) MemoryFail(funcname);

	int nlistloc1 = 0;
	int nlistloc2 = 0;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		if (_imask[jj] == 0) {
			listloc1[nlistloc1] = jj;
			nlistloc1++;
		} else {
			listloc2[nlistloc2] = jj;
			nlistloc2++;
		};
	};

// Main cycle

	int *listval1;
	int *listval2;
	int *listfailed;

	listval1 = new int [n];
	if(!listval1) MemoryFail(funcname);
	listval2 = new int [n];
	if(!listval2) MemoryFail(funcname);
	listfailed = new int [n];
	if(!listfailed) MemoryFail(funcname);

	int kk, ivalue, ival1max, ival2max, ntry1, nadd1, ntry2, nadd2;
	bool bcheck;

	int ithresh, icycle, icycle1;

	int ncyclemax = 3;
//	int ncyclemax = 1;

	int icase = 0;

	for (icycle=0;icycle<ncyclemax;icycle++) {

		ithresh = -1;

		while (true) {

// For each active list value compute its characteristics

			for (i=0;i<nlistloc1;i++) listval1[i] = 0;
			for (i=0;i<nlistloc2;i++) listval2[i] = 0;

			for (i=0;i<nlistloc1;i++) {
				kk = listloc1[i];
				if (imasknd[kk] == -3) {
					ivalue = 0;
					for (j=ia[kk];j<ia[kk+1];j++) {
						jj = ja[j];
						if (_imask[jj] == 1 && imasknd[jj] != -2) {
							ivalue++;
						};
					};
					listval1[i] = ivalue;
				};
			};

			for (i=0;i<nlistloc2;i++) {
				kk = listloc2[i];
				if (imasknd[kk] == -3) {
					ivalue = 0;
					for (j=ia[kk];j<ia[kk+1];j++) {
						jj = ja[j];
						if (_imask[jj] == 0 && imasknd[jj] != -1) {
							ivalue++;
						};
					};
					listval2[i] = ivalue;
				};
			};

// Find maximal characteristics

			ival1max = 0;

			for (i=0;i<nlistloc1;i++) {
				if (listval1[i] > ival1max) ival1max = listval1[i];
			};

			ival2max = 0;

			for (i=0;i<nlistloc2;i++) {
				if (listval2[i] > ival2max) ival2max = listval2[i];
			};

			if (ival1max == 0 && ival2max == 0) break;

// Modify threshold value

			if (ithresh == -1) {
				ithresh = ival1max;
				if (ival2max > ithresh) ithresh = ival2max;
				ithresh++;
			};

			ithresh--;

// Try to add nodes according to the threshold

			ntry1 = 0;
			ntry2 = 0;
			nadd1 = 0;
			nadd2 = 0;

			for (icycle1=0;icycle1<2;icycle1++) {

				icase = (icase+1) % 2;

				switch (icase) {

					case 0:

						for (i=0;i<nlistloc1;i++) {
							kk = listloc1[i];
							ivalue = listval1[i];
							if (ivalue >= ithresh && imasknd[kk] == -3) {
								ntry1++;
								imasknd[kk] = -1;
								bcheck = CheckNeibours (_imask, imasknd, 1, &kk);
								if (!bcheck) {
									imasknd[kk] = -3;
								} else {
									nadd1++;
								};
							};
						};

						break;

					case 1:

						for (i=0;i<nlistloc2;i++) {
							kk = listloc2[i];
							ivalue = listval2[i];
							if (ivalue >= ithresh && imasknd[kk] == -3) {
								ntry2++;
								imasknd[kk] = -2;
								bcheck = CheckNeibours (_imask, imasknd, 1, &kk);
								if (!bcheck) {
									imasknd[kk] = -3;
								} else {
									nadd2++;
								};
							};
						};

						break;

				};

			};

//			cout << "     Ithresh = " << ithresh << " Ntry1 = " << ntry1 << " Nadd1 = " << nadd1;
//			cout << " Ntry2 = " << ntry2 << " Nadd2 = " << nadd2 << endl;

// Check if iterations should be stopped

			if (ithresh == 0) break;

		};

// Check incorrectly inserted boundaries

		int nfailed = 0;

		CheckList (_imask, imasknd,
										nlistloc, listloc, nfailed, listfailed);

//		cout << "    The number of boundaries having no connections = " << nfailed << endl;
//		OutArr (cout," Listfailed_bnd = ",nfailed,listfailed);

// Check that all boundary values has been processed

		CheckWholes (_imask, imasknd, nlistloc, listloc, nfailed, listfailed);

//		cout << "    The number of failed found boundaries = " << nfailed << endl;
//		OutArr (cout," Listfailed_wholes = ",nfailed,listfailed);

		if (nfailed == 0 || icycle == ncyclemax-1) break;

// According to the failed nodes modify status of some boundaries

		for (i=0;i<nfailed;i++) {
			kk = listfailed[i];
			for (j=ia[kk];j<ia[kk+1];j++) {
				jj = ja[j];
				if (imasknd[jj] == -1 || imasknd[jj] == -2) imasknd[jj] = -3;
			};
		};

		icase = (icase+1) % 2;

	};

// Perform some search of the remaining data

	int *imaskloc;
	int *listsmall;
	int *listext;
	int *imaskchk;
	int *listchk;
	int *listsave;

	imaskloc = new int [n];
	if(!imaskloc) MemoryFail(funcname);
	listsmall = new int [n];
	if(!listsmall) MemoryFail(funcname);
	listext = new int [n];
	if(!listext) MemoryFail(funcname);
	imaskchk = new int [n];
	if(!imaskchk) MemoryFail(funcname);
	listchk = new int [n];
	if(!listchk) MemoryFail(funcname);
	listsave = new int [n];
	if(!listsave) MemoryFail(funcname);

	int nlistsmall = 0;

	CheckWholes (_imask, imasknd, nlistloc, listloc, nlistsmall, listsmall);

	for (i=0;i<n;i++) imaskloc[i] = -1;
	for (i=0;i<n;i++) imaskchk[i] = -1;

	int nlistext = 0;

	for (i=0;i<nlistsmall;i++) {
		jj = listsmall[i];
		listsave[jj] = imasknd[jj];
		imaskloc[jj] = nlistext;
		listext[nlistext] = jj;
		nlistext++;
	};

	for (i=0;i<nlistsmall;i++) {
		kk = listsmall[i];
		for (j=ia[kk];j<ia[kk+1];j++) {
			jj = ja[j];
			if (imasknd[jj] == -1 || imasknd[jj] == -2) {
				listsave[jj] = imasknd[jj];
				imasknd[jj] = -3;
			};
			if (imaskloc[jj] == -1 && imasknd[jj] == -3) {
				imaskloc[jj] = nlistext;
				listext[nlistext] = jj;
				nlistext++;
			};
		};
	};

	int *imaskext;
	int *bitmask;

	imaskext = new int [n];
	if(!imaskext) MemoryFail(funcname);
	bitmask = new int [n];
	if(!bitmask) MemoryFail(funcname);

	int icycleext = -1;

	int nlistsave, ind, nlistchk, ntrymax, itry, icheck;
	int nfailed;

	for (i=0;i<nlistext;i++) imaskext[i] = icycleext;

	while (true) {

		icycleext++;

// Find not marked nodes from extended list

		nlistsmall = 0;

		for (i=0;i<nlistext;i++) {
			if (imaskext[i] == -1) {
				imaskext[i] = icycleext;
				listsmall[nlistsmall] = listext[i];
				nlistsmall++;
				break;
			};
		};

		if (nlistsmall == 0) break;

// Find all other connected nodes if any

		while (true) {
			nlistsave = nlistsmall;
			for (i=0;i<nlistsave;i++) {
				kk = listsmall[i];
				for (j=ia[kk];j<ia[kk+1];j++) {
					jj = ja[j];
					ind = imaskloc[jj];
					if (ind > 0) {
						if (imaskext[ind] == -1) {
							imaskext[ind] = icycleext;
							listsmall[nlistsmall] = listext[ind];
							nlistsmall++;
						};
					};
				};
			};
			if (nlistsave == nlistsmall) break;
		};

// Create extended list according to the current list

		ntrymax = 4;

		for (itry=0;itry<ntrymax;itry++) {

			nlistchk = 0;

			for (i=0;i<nlistsmall;i++) {
				jj = listsmall[i];
				imaskchk[jj] = icycleext;
				listchk[nlistchk] = jj;
				nlistchk++;
			};

			for (i=0;i<nlistsmall;i++) {
				kk = listsmall[i];
				for (j=ia[kk];j<ia[kk+1];j++) {
					jj = ja[j];
					if (imasknd[jj] < 0 && imaskchk[jj] != icycleext) {
						imaskchk[jj] = icycleext;
						listchk[nlistchk] = jj;
						nlistchk++;
					};
				};
			};

// Perform binary search

			cout << "    Before binary search: Nlist = " << nlistsmall << " Nlistchk = " << nlistchk << endl;

//			OutArr (cout," Listsmall = ",nlistsmall,listsmall);
//			OutArr (cout," Listchk = ",nlistchk,listchk);

//			bcheck = BinarySearch (_imask, imasknd,
//											nlistsmall, listsmall, nlistchk, listchk,
//											bitmask, listfailed);
			icheck = BinarySearchByPairs (_imask, imasknd,
											nlistsmall, listsmall, nlistchk, listchk,
											bitmask, listfailed);

// Restore previous data if necessary 

			if (icheck != 0 && itry == ntrymax-1) {
				icheck = -1;
			};

			if (icheck < 0) {
				for (i=0;i<nlistsmall;i++) {
					jj = listsmall[i];
					imasknd[jj] = listsave[jj];
				};
				break;
			};

			if (icheck == 0) {
//				for (i=0;i<nlistsmall;i++) {
//					jj = listsmall[i];
//					cout << " Ind = " << jj << " Part = " << _imask[jj] << " Bnd = " << imasknd[jj] << endl;
//				};
			};
			if (icheck == 0) break;

			nlistchk = nlistsmall;

			for (i=0;i<nlistsmall;i++) {
				kk = listsmall[i];
				for (j=ia[kk];j<ia[kk+1];j++) {
					jj = ja[j];
					if (imasknd[jj] == -1 || imasknd[jj] == -2) {
						listsave[jj] = imasknd[jj];
						imasknd[jj] = -3;
					};
					if (imaskloc[jj] == -1 && imasknd[jj] == -3) {
						imaskloc[jj] = nlistext;
						listext[nlistext] = jj;
						nlistext++;
						listsmall[nlistchk] = jj;
						nlistchk++;
					};
				};
			};

			nlistsmall = nlistchk;

			icycleext++;

		};

	};

// Perform final check of the result

	CheckList (_imask, imasknd,
									nlistloc, listloc, nfailed, listfailed);

//	cout << "    Final check: Nfailed bnd = " << nfailed << endl;
//	OutArr (cout," Listfailed_bnd = ",nfailed,listfailed);

	CheckWholes (_imask, imasknd,
						nlistloc, listloc, nfailed, listfailed);
//	cout << "    Final check: Nfailed wholes = " << nfailed << endl;
	for (i=0;i<nfailed;i++) {
		jj = listfailed[i];
		if (imasknd[jj] == -3) {
			ipart = _imask[jj];
			if (ipart == 0) imasknd[jj] = -1;
			if (ipart == 1) imasknd[jj] = -2;
		};
	};
	CheckWholes (_imask, imasknd,
						nlistloc, listloc, nfailed, listfailed);
//	cout << "    Again final check: Nfailed wholes = " << nfailed << endl;

	if (nfailed != 0) {
		cout << "    Error: Fail on check ";
		throw " CSMatrix::OptimalBoundary: failed on check ";
	};

// Create condensed boundary list

	_nlistbnd = 0;

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		if (imasknd[jj] != -3) {
			_listbnd[_nlistbnd] = jj;
			_nlistbnd++;
		};
	};

// Free work arrays

	delete [] imasknd;
	delete [] listloc;
	delete [] listloc1;
	delete [] listloc2;
	delete [] listval1;
	delete [] listval2;
	delete [] listfailed;
	delete [] imaskloc;
	delete [] listsmall;
	delete [] listext;
	delete [] imaskchk;
	delete [] listchk;
	delete [] listsave;
	delete [] imaskext;
	delete [] bitmask;

};

// Author: Kharchenko S.A.
// CSMatrix: For the list of nodes perform binary search of the correct combination
//========================================================================================
bool CSMatrix::BinarySearch (int *_imask, int *_imasknd, // For the list of nodes perform binary search of the correct combination
										int _nlist, int *_list, int _nlistchk, int *_listchk,
										int *_bitmask, int *_listfailed) const {

//	const char *funcname = "BinarySearch";

	bool bcheck = false;

	int i;
	int icurr = _nlist-1;

	for (i=0;i<_nlist;i++) _bitmask[i] = 0;

// Perform binary search

	bool bcheckloc = true;
	int jj, ivalue, nfailed = 0;

	while (true) {

// Transform current bitmask into the data for check

		for (i=0;i<_nlist;i++) {
			jj = _list[i];
			ivalue = _bitmask[i];
			if (ivalue == 0) {
				_imasknd[jj] = -3;
			} else {
				if (_imask[jj] == 0) _imasknd[jj] = -1;
				if (_imask[jj] == 1) _imasknd[jj] = -2;
			};
		};

// Check current combination

		bcheckloc = CheckList (_imask, _imasknd,
										_nlistchk, _listchk);
		if (bcheckloc) {
			CheckWholes (_imask, _imasknd,
				_nlistchk, _listchk, nfailed, _listfailed);
		} else {
			nfailed = 1;
		};

		if (bcheckloc && nfailed == 0) {
			bcheck = true;
			break;
		};

// Modify the bitmask data

		if (_bitmask[icurr] == 0) {
			_bitmask[icurr] = 1;
			for (i=icurr+1;i<_nlist;i++) _bitmask[i] = 0;
			icurr = _nlist-1;
		} else {
			while (icurr > 0 && _bitmask[icurr] == 1) icurr--;
			if (icurr == 0 && _bitmask[icurr] == 1) break;
			_bitmask[icurr] = 1;
			for (i=icurr+1;i<_nlist;i++) _bitmask[i] = 0;
			icurr = _nlist-1;
		};
	};

// Restore original situation if failed

	if (!bcheck) {
		for (i=0;i<_nlist;i++) {
			jj = _list[i];
			_imasknd[jj] = -3;
		};
	};

	return bcheck;

};

// Author: Kharchenko S.A.
// CSMatrix: For the list of nodes perform binary search of the correct combination
//========================================================================================
int CSMatrix::BinarySearchByPairs (int *_imask, int *_imasknd, // For the list of nodes perform binary search of the correct combination
												int _nlist, int *_list, int _nlistchk, int *_listchk,
												int *_bitmask, int *_listfailed) const {

	const char *funcname = "BinarySearchByPairs";

// Allocate work memory

	int *iapairs;
	int *japairs;
	int *ind2list;
	int *imasklist;
	int *bitmask2;
	int *listloc;

	iapairs = new int [_nlist+1];
	if(!iapairs) MemoryFail(funcname);
	japairs = new int [_nlist];
	if(!japairs) MemoryFail(funcname);
	ind2list = new int [n];
	if(!ind2list) MemoryFail(funcname);
	imasklist = new int [_nlist];
	if(!imasklist) MemoryFail(funcname);
	bitmask2 = new int [_nlist];
	if(!bitmask2) MemoryFail(funcname);
	listloc = new int [_nlist];
	if(!listloc) MemoryFail(funcname);

	int i, jj;

	for (i=0;i<n;i++) ind2list[i] = -1;

	for (i=0;i<_nlist;i++) {
		jj = _list[i];
		ind2list[jj] = i;
		imasklist[i] = -1;
	};

// Split current list into the pairs

	int ilist, ipart1, ilist1;
	int ipart, j, kk;

	int npairs = 0;
	int nzpairs = 0;

	iapairs[0] = 0;

	while (true) {

// Find first node if any

		ilist = -1;
		jj = -1;
		for (i=0;i<_nlist;i++) {
			if (imasklist[i] == -1) {
				ilist = i;
				jj = _list[i];
				break;
			};
		};
		if (ilist == -1) break;

// Register current node

		imasklist[ilist] = npairs;
		japairs[nzpairs] = jj;
		nzpairs++;

// Try to find a pair

		ipart = _imask[jj];

		for (j=ia[jj];j<ia[jj+1];j++) {
			kk = ja[j];
			ipart1 = _imask[kk];
			if (ipart != ipart1) {
				ilist1 = ind2list[kk];
				if (ilist1 >= 0) {
					if (imasklist[ilist1] == -1) {
						imasklist[ilist1] = npairs;
						japairs[nzpairs] = kk;
						nzpairs++;
						break;
					};
				};
			};
		};

		iapairs[npairs+1] = nzpairs;
		npairs++;

	};

	cout << "         Npairs = " << npairs << endl;

//	OutArr (cout," IaPairs = ",npairs+1,iapairs);
//	OutArr (cout," JaPairs = ",_nlist,japairs);

// Reproduce the modified list

	int npairsloc = 0;
	int nlistloc = 0;

	int ni, ibeg;

	for (i=0;i<npairs;i++) {
		ni = iapairs[i+1]-iapairs[i];
		ibeg = iapairs[i];
		if (ni == 2) {
			listloc[nlistloc] = japairs[ibeg];
			nlistloc++;
			listloc[nlistloc] = japairs[ibeg+1];
			nlistloc++;
			npairsloc++;
		};
	};

	for (i=0;i<npairs;i++) {
		ni = iapairs[i+1]-iapairs[i];
		ibeg = iapairs[i];
		if (ni == 1) {
			listloc[nlistloc] = japairs[ibeg];
			nlistloc++;
		};
	};

// Perform pairs binary search

	int iend_pairs = npairsloc*2;

	int icheck = 1;

	int icurr = npairs-1;

	for (i=0;i<npairs;i++) bitmask2[i] = 0;

	int npairsmax = 20;

	if (npairs > npairsmax) icheck = -1;

// Perform binary search

	bool bcheckloc = true;
	int ivalue, nfailed = 0;

	while (true && icheck > 0) {

// Create bitmask array

		for (i=0;i<npairsloc;i++) {
			if (bitmask2[i] == 0) {
				_bitmask[i*2] = 0;
				_bitmask[i*2+1] = 1;
			} else if (bitmask2[i] == 1) {
				_bitmask[i*2] = 1;
				_bitmask[i*2+1] = 0;
			} else if (bitmask2[i] == 2) {
				_bitmask[i*2] = 1;
				_bitmask[i*2+1] = 1;
			};
		};

		for (i=iend_pairs;i<_nlist;i++) _bitmask[i] = bitmask2[i+npairsloc-iend_pairs];

// Transform current bitmask into the data for check

		for (i=0;i<_nlist;i++) {
			jj = listloc[i];
			ivalue = _bitmask[i];
			if (ivalue == 0) {
				_imasknd[jj] = -3;
			} else {
				if (_imask[jj] == 0) _imasknd[jj] = -1;
				if (_imask[jj] == 1) _imasknd[jj] = -2;
			};
		};

// Check current combination

		bcheckloc = CheckList (_imask, _imasknd,
										_nlistchk, _listchk);
		if (bcheckloc) {
			CheckWholes (_imask, _imasknd,
				_nlistchk, _listchk, nfailed, _listfailed);
		} else {
			nfailed = 1;
		};

		if (bcheckloc && nfailed == 0) {
			icheck = 0;
			break;
		};

// Modify the bitmask2 data

		if (icurr >= npairsloc) {
			if (bitmask2[icurr] == 0) {
				bitmask2[icurr] = 1;
				for (i=icurr+1;i<npairs;i++) bitmask2[i] = 0;
				icurr = npairs-1;
			} else {
				while (icurr > 0 && bitmask2[icurr] == 1 && icurr >= npairsloc) icurr--;
				if (icurr < npairsloc) {
					while (icurr > 0 && bitmask2[icurr] == 2) icurr--;
					if (icurr == 0 && bitmask2[icurr] == 2) break;
					bitmask2[icurr]++;
					for (i=icurr+1;i<npairs;i++) bitmask2[i] = 0;
					icurr = npairs-1;
				} else {
					if (icurr == 0 && bitmask2[icurr] == 1) break;
					bitmask2[icurr] = 1;
					for (i=icurr+1;i<npairs;i++) bitmask2[i] = 0;
					icurr = npairs-1;
				};
			};
		} else {
			if (bitmask2[icurr] == 0) {
				bitmask2[icurr] = 1;
				for (i=icurr+1;i<npairs;i++) bitmask2[i] = 0;
				icurr = npairs-1;
			} else if (bitmask2[icurr] == 1) {
				bitmask2[icurr] = 2;
				for (i=icurr+1;i<npairs;i++) bitmask2[i] = 0;
				icurr = npairs-1;
			} else if (bitmask2[icurr] == 2) {
				while (icurr > 0 && bitmask2[icurr] == 2) icurr--;
				if (icurr == 0 && bitmask2[icurr] == 2) break;
				bitmask2[icurr]++;
				for (i=icurr+1;i<npairs;i++) bitmask2[i] = 0;
				icurr = npairs-1;
			};
		};
	};

// Restore original situation if failed

	if (icheck != 0) {
		for (i=0;i<_nlist;i++) {
			jj = _list[i];
			_imasknd[jj] = -3;
		};
	};

// Free work data

	delete [] iapairs;
	delete [] japairs;
	delete [] ind2list;
	delete [] imasklist;
	delete [] bitmask2;
	delete [] listloc;

	return icheck;

};

// Author: Kharchenko S.A.
// CSMatrix: For the list of nodes check their neibours for correctness
//========================================================================================
bool CSMatrix::CheckNeibours (int *_imask, int *_imasknd, // For the list of nodes check their neibours for correctness
											int _nlist, int *_list) const {

//	const char *funcname = "CheckNeibours";

	int ilist, k, kk, kkk, j, jj, ipart, ipart1, icheck1loc, icheck2loc;

	bool bcheck = true;

	for (ilist=0;ilist<_nlist;ilist++) {

		kk = _list[ilist];

		for (j=ia[kk];j<ia[kk+1];j++) {
			jj = ja[j];
			ipart = _imask[jj];
			ipart1 = (ipart+1)%2;
			if (_imasknd[jj] == -1 || _imasknd[jj] == -2) {
				icheck1loc = 0;
				icheck2loc = 0;
				for (k=ia[jj];k<ia[jj+1];k++) {
					kkk = ja[k];
					if ((_imask[kkk] == ipart && _imasknd[kkk] == -3) || _imasknd[kkk] == ipart) icheck1loc = 1;
					if ((_imask[kkk] == ipart1 && _imasknd[kkk] == -3) || _imasknd[kkk] == ipart1) icheck2loc = 1;
				};
				if (icheck1loc != 1 || icheck2loc != 1) {
					bcheck = false;
					return bcheck;
				};
			};
		};
	};

	return bcheck;

};

// Author: Kharchenko S.A.
// CSMatrix: For the list of nodes check their correctness
//========================================================================================
bool CSMatrix::CheckList (int *_imask, int *_imasknd, // For the list of nodes check their correctness
											int _nlist, int *_list) const {

//	const char *funcname = "CheckList";

	int ilist, k, kkk, jj, ipart, ipart1, icheck1loc, icheck2loc;

	bool bcheck = true;

	for (ilist=0;ilist<_nlist;ilist++) {

		jj = _list[ilist];

		ipart = _imask[jj];
		ipart1 = (ipart+1)%2;
		if (_imasknd[jj] == -1 || _imasknd[jj] == -2) {
			icheck1loc = 0;
			icheck2loc = 0;
			for (k=ia[jj];k<ia[jj+1];k++) {
				kkk = ja[k];
				if ((_imask[kkk] == ipart && _imasknd[kkk] == -3) || _imasknd[kkk] == ipart) icheck1loc = 1;
				if ((_imask[kkk] == ipart1 && _imasknd[kkk] == -3) || _imasknd[kkk] == ipart1) icheck2loc = 1;
			};
			if (icheck1loc != 1 || icheck2loc != 1) {
				bcheck = false;
				return bcheck;
			};
		};

	};

	return bcheck;

};

// Author: Kharchenko S.A.
// CSMatrix: For the list of nodes check their correctness
//========================================================================================
void CSMatrix::CheckList (int *_imask, int *_imasknd, // For the list of nodes check their correctness
											int _nlist, int *_list, int &_nfailed, int *_listfailed) const {

//	const char *funcname = "CheckList";

	int ilist, k, kkk, jj, ipart, ipart1, icheck1loc, icheck2loc;

	_nfailed = 0;

	for (ilist=0;ilist<_nlist;ilist++) {

		jj = _list[ilist];

		ipart = _imask[jj];
		ipart1 = (ipart+1)%2;
		if (_imasknd[jj] == -1 || _imasknd[jj] == -2) {
			icheck1loc = 0;
			icheck2loc = 0;
			for (k=ia[jj];k<ia[jj+1];k++) {
				kkk = ja[k];
				if ((_imask[kkk] == ipart && _imasknd[kkk] == -3) || _imasknd[kkk] == ipart) icheck1loc = 1;
				if ((_imask[kkk] == ipart1 && _imasknd[kkk] == -3) || _imasknd[kkk] == ipart1) icheck2loc = 1;
			};
			if (icheck1loc != 1 || icheck2loc != 1) {
				_listfailed[_nfailed] = jj;
				_nfailed++;
			};
		};

	};

};

// Author: Kharchenko S.A.
// CSMatrix: For the list of nodes check that there are no wholes
//========================================================================================
void CSMatrix::CheckWholes (int *_imask, int *_imasknd, // For the list of nodes check that there are no wholes
										int _nlist, int *_list, int &_nfailed, int *_listfailed) const {

//	const char *funcname = "CheckWholes";

	int i, kk, ipart, ivalue, j, jj;

	_nfailed = 0;

	for (i=0;i<_nlist;i++) {
		kk = _list[i];
		ipart = _imask[kk];
		if (_imasknd[kk] == -3) {
			ivalue = 0;
			for (j=ia[kk];j<ia[kk+1];j++) {
				jj = ja[j];
				if (ipart == 0) {
					if (_imask[jj] == 1 && _imasknd[jj] != -2) {
						ivalue++;
					};
				} else {
					if (_imask[jj] == 0 && _imasknd[jj] != -1) {
						ivalue++;
					};
				};
			};
			if (ivalue != 0) {
				_listfailed[_nfailed] = kk;
				_nfailed++;
//					throw " CSMatrix::OptimalBoundary: not all boundaries were processed ";
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrix: Compute the upper triangular part of the sparsity of the AtA matrix
//========================================================================================
CSMatrix CSMatrix::AtAMatrSpars () const { // Compute the upper triangular part of the sparsity of the AtA matrix

	const char *funcname = "AtAMatrSpars";

// Compute transposed to the matrix

	CSMatrix at;

	at = TranspMtrRect ();

// Count the number of elements in AtA matrix

	int *imask;
	int *lstloc;

	int nsupmx = nsupr;

	if (nsupc < nsupmx) nsupmx = nsupc;

	imask = new int [nsupmx];
	if (!imask) MemoryFail (funcname);
	lstloc = new int [nsupmx];
	if (!lstloc) MemoryFail (funcname);

	int icycle = -1;
	int i;
	for (i=0;i<nsupmx;i++) {
		imask[i] = icycle;
	};

	int nz = 0;

	int isup;
	for (isup=0;isup<nsupc;isup++) {
		icycle++;
		for (int j=ia[isup];j<ia[isup+1];j++) {
			int jj = ja[j];
			for (int k=at.ia[jj];k<at.ia[jj+1];k++) {
				int kk = at.ja[k];
				if (imask[kk] != icycle) {
					nz++;
					imask[kk] = icycle;
				};
			};
		};
	};

// Allocate local ia and ja arrays

	int *ialoc;
	int *jaloc;

	ialoc = new int [nsupc+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nz];
	if (!jaloc) MemoryFail (funcname);

// Fill arrays

	ialoc[0] = 0;
	nz = 0;

	for (isup=0;isup<nsupc;isup++) {
		icycle++;
		int nzloc = 0;
		int j;
		for (j=ia[isup];j<ia[isup+1];j++) {
			int jj = ja[j];
			for (int k=at.ia[jj];k<at.ia[jj+1];k++) {
				int kk = at.ja[k];
				if (imask[kk] != icycle) {
					imask[kk] = icycle;
					lstloc[nzloc] = kk;
					nzloc++;
				};
			};
		};
		if (nzloc != 0) qsort (lstloc, nzloc, sizeof(int), compint);
		for (j=0;j<nzloc;j++) jaloc[nz+j] = lstloc[j];
		nz += nzloc;
		ialoc[isup+1] = nz;
	};

// Filter lower triangular part

	int *ianew;

	ianew = new int [nsupc+1];
	if (!ianew) MemoryFail (funcname);

	ianew[0] = 0;
	nz = 0;

	for (isup=0;isup<nsupc;isup++) {
		for (int j=ialoc[isup];j<ialoc[isup+1];j++) {
			int jj = jaloc[j];
			if (jj >= isup) {
				jaloc[nz] = jaloc[j];
				nz++;
			};
		};
		ianew[isup+1] = nz;
	};

// Create ata structure

	CSMatrix ata (nsupc,nz);

	for (i=0;i<nsupc;i++) ata.list[i] = i;
	for (i=0;i<=nsupc;i++) ata.ia[i] = ianew[i];
	for (i=0;i<nz;i++) ata.ja[i] = jaloc[i];

	ata.m = nsupc;
	ata.n = nsupc;
	ata.nsupc = nsupc;
	ata.nsupr = nsupc;
	ata.nlist = nsupc;
	ata.nzja = nz;

// Free work arrays

	delete [] imask;
	delete [] lstloc;
	delete [] ialoc;
	delete [] ianew;
	delete [] jaloc;

	return ata;

};

// Author: Kharchenko S.A.
// CSMatrix: Find separators of the symmetric square matrix if any
//========================================================================================
void CSMatrix::FindSeparator (int &_ind1, int &_ind2) const { // Find separators of the symmetric square matrix if any

//	const char *funcname = "FindSeparator";

	_ind1 = -1;
	_ind2 = -1;

// Initial guesses for the separators

	int ind1beg, ind1end;
	int ind2beg, ind2end;

	ind1beg = (1*n) / 4;
	ind1end = (3*n) / 4;

	ind2beg = ind1end+1;
	ind2end = n;

	if (ind1beg >= ind1end) {
		return;
	};
	if (ind2beg >= ind2end) {
		return;
	};
	if (ind1end >= ind2beg) {
		return;
	};

// Check that initial guess is a valid one

	int i, j, jj;

	for (i=ind1end;i<ind2beg;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj<ind1beg) return;
		};
	};

// Modify ind1end

	int ind1endnew = ind1end;

//	for (i=ind2beg-1;i>=ind1end;i--) {
	for (i=ind2beg-1;i>=ind1endnew;i--) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj < ind1endnew-1) ind1endnew = jj+1;
		};
	};

// Check that modified guess for ind1end is a valid one

	ind1end = ind1endnew;

	for (i=ind1end;i<ind2beg;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj<ind1beg) return;
		};
	};

// Improve ind1beg

	int ind1begnew = ind1beg;

//	for (i=0;i<ind1beg;i++) {
	for (i=0;i<ind1begnew;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj<ind2beg) {
//				if (jj > ind1begnew-1) ind1begnew = jj+1;
				if (jj > ind1begnew) ind1begnew = jj;
			};
		};
	};

// Check that modified guess for ind1beg is a valid one

	ind1beg = ind1begnew;

	for (i=ind1end;i<ind2beg;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj<ind1beg) return;
		};
	};

// Now perform not more than log2n iterations to find exact block splitting

	int ind1check;

	while (ind1beg != ind1end) {
		if (ind1end == ind1beg+1) {
			ind1check = ind1beg;
		} else {
			ind1check = (ind1beg+ind1end) / 2;
		};
		int icheck = 1;
		for (i=ind1check;i<ind2beg;i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				if (jj<ind1check) {
					icheck = 0;
				};
			};
		};
		if (icheck == 1) {
			if (ind1end != ind1beg+1) {
				ind1beg = ind1check;
			} else {
				ind1end = ind1check;
			};
		} else {
			ind1end = ind1check;
		};
	};

	_ind1 = ind1beg;

// Compute ind2

	_ind2 = n;

	for (i=0;i<_ind1;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj>_ind1+1) {
				if (jj < _ind2) _ind2 = jj;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrix: Matrix ordering routine
//========================================================================================
CSMatrix CSMatrix::OrdMtr (const int *_order) const { // Matrix ordering routine

	const char *funcname = "OrdMtr";

// Compute inverse order

	int *iord;
	iord = new int [n];
	if (!iord) MemoryFail (funcname);

	int i;
	for (i=0;i<n;i++) iord[_order[i]] = i;

// Reorder the matrix

	int nzjas = nzja;

	CSMatrix an(n,nzjas);

	for (i=0;i<=n;i++) an.ia[i] = 0;

	int j;

	for (i=0;i<n;i++) {
		j = _order[i];
		an.ia[j+1] = ia[i+1]-ia[i];
	};

	for (i=0;i<n;i++) an.ia[i+1] += an.ia[i];

	int nz = 0;

	int inew, nzloc, ibs, jold, jnew;
	int *pja;

	for (inew=0;inew<n;inew++) {
		i = iord[inew];
		nzloc = ia[i+1]-ia[i];
		ibs = nz;
		for (j=ia[i];j<ia[i+1];j++) {
			jold = ja[j];
			jnew = _order[jold];
			an.ja[nz] = jnew;
			nz++;
		};

// Sort elements

		pja = an.ja;

		if (nzloc != 0) qsort (pja+ibs, nzloc, sizeof(int), compint);

	};

// Delete ordering array

	delete [] iord;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrix: Sort matrix data by ordering the list (not indices)
//========================================================================================
CSMatrix CSMatrix::SortDataViaList (const int *_orderr) const { // CSMatrix: Sort matrix data by ordering the list (not indices)

//	const char *funcname = "SortDataViaList";

// Reorder the matrix

	int i, j;

	CSMatrix an(nlist,nzja);

	for (i=0;i<=nlist;i++) an.ia[i] = 0;

	for (i=0;i<nlist;i++) {
		j = _orderr[i];
		an.ia[j+1] = ia[i+1]-ia[i];
		an.list[j] = list[i];
	};

	for (i=0;i<nlist;i++) an.ia[i+1] += an.ia[i];

	int inew, ibs;

	for (i=0;i<nlist;i++) {
		inew = _orderr[i];
		ibs = an.ia[inew];
		for (j=ia[i];j<ia[i+1];j++) {
			an.ja[ibs+j-ia[i]] = ja[j];
		};
	};

	an.n = n;
	an.m = m;
	an.nsupr = nsupr;
	an.nsupc = nsupc;
	an.nlist = nlist;
	an.nzja = nzja;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrix: Sort matrix data by ordering the list (not indices)
//========================================================================================
CSMatrix CSMatrix::Sort2IndexDataViaList (const int *_orderr) const { // Sort matrix data by ordering the list (not indices)

//	const char *funcname = "Sort2IndexDataViaList";

// Reorder the matrix

	int i, j;

	CSMatrix an(nlist,nlist,nzja,nzja);

	for (i=0;i<=nlist;i++) an.ia[i] = 0;

	for (i=0;i<nlist;i++) {
		j = _orderr[i];
		an.ia[j+1] = ia[i+1]-ia[i];
		an.list[j] = list[i];
		an.list2[j] = list2[i];
	};

	for (i=0;i<nlist;i++) an.ia[i+1] += an.ia[i];

	int inew, ibs;

	for (i=0;i<nlist;i++) {
		inew = _orderr[i];
		ibs = an.ia[inew];
		for (j=ia[i];j<ia[i+1];j++) {
			an.ja[ibs+j-ia[i]] = ja[j];
			an.ja2[ibs+j-ia[i]] = ja2[j];
		};
	};

	an.n = n;
	an.m = m;
	an.nsupr = nsupr;
	an.nsupc = nsupc;
	an.nlist = nlist;
	an.nlist2 = nlist2;
	an.nzja = nzja;
	an.nzja2 = nzja2;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrix: Symmetric matrix ordering routine
//========================================================================================
CSMatrix CSMatrix::OrdMtrSymm (const int *_order) const { // Symmetric matrix ordering routine

	const char *funcname = "OrdMtrSymm";

// Compute inverse order

	int *iord;
	iord = new int [n];
	if (!iord) MemoryFail (funcname);

	int i;
	for (i=0;i<n;i++) iord[_order[i]] = i;

// Compute the symmetric matrix

	CSMatrix as;

	as = SymmMtr ();

// Reorder the matrix

	int nzjas = as.nzja;

	CSMatrix an(n,nzjas);

	for (i=0;i<=n;i++) an.ia[i] = 0;

	int j;

	for (i=0;i<n;i++) {
		j = _order[i];
		an.ia[j+1] = as.ia[i+1]-as.ia[i];
	};

	for (i=0;i<n;i++) an.ia[i+1] += an.ia[i];

	int nz = 0;

	int inew, nzloc, ibs, jold, jnew;
	int *pja;

	for (inew=0;inew<n;inew++) {
		i = iord[inew];
		nzloc = as.ia[i+1]-as.ia[i];
		ibs = nz;
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jold = as.ja[j];
			jnew = _order[jold];
			an.ja[nz] = jnew;
			nz++;
		};

// Sort elements

		pja = an.ja;

		if (nzloc != 0) qsort (pja+ibs, nzloc, sizeof(int), compint);

	};

// Filter unnecessary elements

	CSMatrix at(n,nzja);

	at.ia[0] = 0;
	nz = 0;

	int k;

	for (i=0;i<n;i++) {
		for (j=an.ia[i];j<an.ia[i+1];j++) {
			k = an.ja[j];
			if (k >= i) {
				at.ja[nz] = an.ja[j];
				nz++;
			};
		};
		at.ia[i+1] = nz;
	};

// Delete ordering array

	delete [] iord;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrix: Multiply two matrices in supersparse format stored by rows
//========================================================================================
CSMatrix CSMatrix::SuperSparseMultiplyByRows (int &_icycle, int *_imask, // Multiply two matrices in supersparse format stored by rows
																int &_icycle1, int *_imask1, 
																CSMatrix &_amatr2) const {

	const char *funcname = "SuperSparseMultiplyByRows";

// Main cycle over the rows

	int nlist2ini = _amatr2.GetNlist ();
	int *plist2 = _amatr2.GetList ();
	int *pia2 = _amatr2.GetIa ();
	int *pja2 = _amatr2.GetJa ();

	int ntot = n;

	int *listloc;
	int *list1loc;
	int *list2loc;

	listloc = new int [ntot];
	if (!listloc) MemoryFail (funcname);
	list1loc = new int [ntot];
	if (!list1loc) MemoryFail (funcname);
	list2loc = new int [ntot];
	if (!list2loc) MemoryFail (funcname);

	_icycle1++;

	int i, jj;

	for (i=0;i<nlist2ini;i++) {
		jj = plist2[i];
		list1loc[jj] = i;
		_imask1[jj] = _icycle1;
	};

	int *ialoc;

	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);

	int nzjaloc = 0;
	ialoc[0] = 0;
	int nlist2loc = 0;

	int j, ind1, k, kk, ibeg, nzjalocsave, nlistloc;

	for (i=0;i<nlist;i++) {
		_icycle++;
		nzjalocsave = nzjaloc;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (_imask1[jj] == _icycle1) {
				ind1 = list1loc[jj];
				for (k=pia2[ind1];k<pia2[ind1+1];k++) {
					kk = pja2[k];
					if (_imask[kk] != _icycle) {
						nzjaloc++;
						_imask[kk] = _icycle;
					};
				};
			};
		};
		if (nzjalocsave < nzjaloc) {
			ialoc[nlist2loc+1] = nzjaloc;
			list2loc[nlist2loc] = list[i];
			nlist2loc++;
		};
	};

	int *jaloc;

	jaloc = new int [nzjaloc];
	if (!jaloc) MemoryFail (funcname);

	nlist2loc = 0;

	for (i=0;i<nlist;i++) {
		_icycle++;
		nlistloc = 0;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (_imask1[jj] == _icycle1) {
				ind1 = list1loc[jj];
				for (k=pia2[ind1];k<pia2[ind1+1];k++) {
					kk = pja2[k];
					if (_imask[kk] != _icycle) {
						listloc[nlistloc] = kk;
						nlistloc++;
						_imask[kk] = _icycle;
					};
				};
			};
		};
		if (nlistloc > 0) {
			sort (listloc, listloc + nlistloc);
			ibeg = ialoc[nlist2loc];
			for (j=0;j<nlistloc;j++) jaloc[ibeg+j] = listloc[j];
			nlist2loc++;
		};
	};

// Create new matrix

	CSMatrix temp (nlist2loc,nzjaloc);

// Store new data

	for (i=0;i<nlist2loc;i++) temp.list[i] = list2loc[i];
	for (i=0;i<=nlist2loc;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) temp.ja[i] = jaloc[i];

	temp.m      = m;
	temp.n      = n;
	temp.nsupc  = nsupc;
	temp.nsupr  = nsupr;
	temp.nlist  = nlist2loc;
	temp.nzja   = nzjaloc;

// Free work memory

	delete [] listloc;
	delete [] list1loc;
	delete [] list2loc;
	delete [] ialoc;
	delete [] jaloc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrix: Multiply two matrices in supersparse format stored by rows with block list information
//========================================================================================
void CSMatrix::SuperSparseMultiplyByRowsWithBlockList (int &_icycle, int *_imask, // Multiply two matrices in supersparse format stored by rows with block list information
																			int &_icycle1, int *_imask1, 
																			CSMatrix &_amatr2,
																			int _iblkstr, CSMatrix &_amatr2str,
																			CSMatrix &_mvm, CSMatrix &_mvmstr) const {

	const char *funcname = "SuperSparseMultiplyByRowsWithBlockList";

// Perform multiplication

	_mvm = SuperSparseMultiplyByRows (_icycle, _imask, _icycle1, _imask1, 
													_amatr2);

// Store index references in list arrays

	int ntot = n;

	int *list2ind;
	int *list2ind2;

	list2ind = new int [ntot];
	if (!list2ind) MemoryFail (funcname);
	list2ind2 = new int [ntot];
	if (!list2ind2) MemoryFail (funcname);

	int i, j, irow, icol;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		list2ind[irow] = i;
	};

	int nlist2loc = _amatr2.GetNlist ();
//	int nzja2loc = _amatr2.GetNzja ();
	int *plist2 = _amatr2.GetList ();
	int *pia2 = _amatr2.GetIa ();
	int *pja2 = _amatr2.GetJa ();
	int nzja2strloc = _amatr2str.GetNzja ();
	int *pia2str = _amatr2str.GetIa ();
	int *pja2str = _amatr2str.GetJa ();

	for (i=0;i<nlist2loc;i++) {
		irow = plist2[i];
		list2ind2[irow] = i;
	};

	int iblkmax = 0;

	int iblk;

	for (i=0;i<nzja2strloc;i++) {
		iblk = pja2str[i];
		if (iblk > iblkmax) iblkmax = iblk;
	};

	iblkmax += 2;

	int *listblk;

	listblk = new int [iblkmax];
	if (!listblk) MemoryFail (funcname);

// Find reference to the transposed structure

	CSMatrix amatr2t;

	amatr2t = _amatr2.TranspMtrMask (_icycle, _imask);

	int nzja2t = amatr2t.GetNzja ();
	int nlist2tloc = amatr2t.GetNlist ();
	int *plist2t = amatr2t.GetList ();
	int *pia2t = amatr2t.GetIa ();
	int *pja2t = amatr2t.GetJa ();

	int *ja2tind;

	ja2tind = new int [nzja2t];
	if (!ja2tind) MemoryFail (funcname);

	int ind, indnew, k, kk;

	for (i=0;i<nlist2tloc;i++) {
		icol = plist2t[i];
		for (j=pia2t[i];j<pia2t[i+1];j++) {
			irow = pja2t[j];
			ind = list2ind2[irow];
			indnew = -1;
			for (k=pia2[ind];k<pia2[ind+1];k++) {
				kk = pja2[k];
				if (kk == icol) indnew = k;
			};
			if (indnew < 0) throw " CSMatrix::SuperSparseMultiplyByRowsWithBlockList: index not found ";
			ja2tind[j] = indnew;
		};
	};

	for (i=0;i<nlist2tloc;i++) {
		irow = plist2t[i];
		list2ind2[irow] = i;
	};

// Allocate the data for blocks computation

	int nzjamvm = _mvm.GetNzja ();
	int nlistmvmloc = _mvm.GetNlist ();
	int *plistmvm = _mvm.GetList ();
	int *piamvm = _mvm.GetIa ();
	int *pjamvm = _mvm.GetJa ();

	int *iablkstr;

	iablkstr = new int [nzjamvm+1];
	if (!iablkstr) MemoryFail (funcname);

	for (i=0;i<=nzjamvm;i++) iablkstr[i] = 0;

	int indrow, ind2col, nlistblk, ibeg1, iend1, ibeg2, iend2, ip1, ip2;
	int kk1, kk2, kblk, ibeg;

	for (i=0;i<nlistmvmloc;i++) {
		irow = plistmvm[i];
		indrow = list2ind[irow];
		for (j=piamvm[i];j<piamvm[i+1];j++) {
			icol = pjamvm[j];
			ind2col = list2ind2[icol];

			_icycle++;

			_imask[_iblkstr] = _icycle;
			listblk[0] = _iblkstr;
			nlistblk = 1;

			ibeg1 = ia[indrow];
			iend1 = ia[indrow+1]-1;
			ibeg2 = pia2t[ind2col];
			iend2 = pia2t[ind2col+1]-1;
			ip1 = ibeg1;
			ip2 = ibeg2;
			while (ip1 <= iend1 && ip2 <= iend2) {
				kk1 = ja[ip1];
				kk2 = pja2t[ip2];
				if (kk1 == kk2) {
					ind = ja2tind[ip2];
					for (k=pia2str[ind];k<pia2str[ind+1];k++) {
						kblk = pja2str[k];
						if (_imask[kblk] != _icycle) {
							listblk[nlistblk] = kblk;
							nlistblk++;
							_imask[kblk] = _icycle;
						};
					};
					ip1++;
					ip2++;
				} else if (kk1 < kk2) {
					ip1++;
				} else if (kk1 > kk2) {
					ip2++;
				};
			};
			iablkstr[j+1] = nlistblk;
		};
	};

	for (i=0;i<nzjamvm;i++) iablkstr[i+1] = iablkstr[i]+iablkstr[i+1];

	int nzjablkstr = iablkstr[nzjamvm];

	int *jablkstr;

	jablkstr = new int [nzjablkstr];
	if (!jablkstr) MemoryFail (funcname);

	for (i=0;i<nlistmvmloc;i++) {
		irow = plistmvm[i];
		indrow = list2ind[irow];
		for (j=piamvm[i];j<piamvm[i+1];j++) {
			icol = pjamvm[j];
			ind2col = list2ind2[icol];

			_icycle++;

			_imask[_iblkstr] = _icycle;
			listblk[0] = _iblkstr;
			nlistblk = 1;

			ibeg1 = ia[indrow];
			iend1 = ia[indrow+1]-1;
			ibeg2 = pia2t[ind2col];
			iend2 = pia2t[ind2col+1]-1;
			ip1 = ibeg1;
			ip2 = ibeg2;
			while (ip1 <= iend1 && ip2 <= iend2) {
				kk1 = ja[ip1];
				kk2 = pja2t[ip2];
				if (kk1 == kk2) {
					ind = ja2tind[ip2];
					for (k=pia2str[ind];k<pia2str[ind+1];k++) {
						kblk = pja2str[k];
						if (_imask[kblk] != _icycle) {
							listblk[nlistblk] = kblk;
							nlistblk++;
							_imask[kblk] = _icycle;
						};
					};
					ip1++;
					ip2++;
				} else if (kk1 < kk2) {
					ip1++;
				} else if (kk1 > kk2) {
					ip2++;
				};
			};
			sort (listblk, listblk + nlistblk);
			ibeg = iablkstr[j];
			for (k=0;k<nlistblk;k++) jablkstr[ibeg+k] = listblk[k];
		};
	};

// Store the result

	CSMatrix tempstr (nzjamvm,nzjablkstr);

	for (i=0;i<nzjamvm;i++) tempstr.list[i] = i;
	for (i=0;i<=nzjamvm;i++) tempstr.ia[i] = iablkstr[i];
	for (i=0;i<nzjablkstr;i++) tempstr.ja[i] = jablkstr[i];

	tempstr.m      = nzjamvm;
	tempstr.n      = nzjamvm;
	tempstr.nsupc  = nzjamvm;
	tempstr.nsupr  = nzjamvm;
	tempstr.nlist  = nzjamvm;
	tempstr.nzja   = nzjablkstr;

	_mvmstr = tempstr;

// Free work memory

	delete [] list2ind;
	delete [] list2ind2;
	delete [] listblk;
	delete [] ja2tind;
	delete [] iablkstr;
	delete [] jablkstr;

};

// Author: Kharchenko S.A.
// CSMatrix: Add the sparsity in supersparse format with block lists
//========================================================================================
void CSMatrix::SuperSparseAddWithBlockList (int &_icycle, int *_imask, // Add the sparsity in supersparse format with block lists
																int &_icycle1, int *_imask1, 
																CSMatrix &_ablkstr,
																CSMatrix &_amtradd, CSMatrix &_ablkstradd) const {

	const char *funcname = "SuperSparseAddWithBlockList";

// Main cycle over the rows

	int nlist2loc = _amtradd.GetNlist ();
	int *plist2 = _amtradd.GetList ();
	int *pia2 = _amtradd.GetIa ();
	int *pja2 = _amtradd.GetJa ();

	int ntot = n;

	int *listloc;
	int *list1loc;
	int *list2loc;
	int *list3loc;

	listloc = new int [ntot];
	if (!listloc) MemoryFail (funcname);
	list1loc = new int [ntot];
	if (!list1loc) MemoryFail (funcname);
	list2loc = new int [ntot];
	if (!list2loc) MemoryFail (funcname);
	list3loc = new int [ntot];
	if (!list3loc) MemoryFail (funcname);

	int nlistloc = 0;

	_icycle++;

	int i, jj;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		listloc[nlistloc] = jj;
		_imask[jj] = _icycle;
		nlistloc++;
	};

	for (i=0;i<nlist2loc;i++) {
		jj = plist2[i];
		if (_imask[jj] != _icycle) {
			listloc[nlistloc] = jj;
			_imask[jj] = _icycle;
			nlistloc++;
		};
	};

	sort (listloc, listloc + nlistloc);

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		list1loc[jj] = i;
	};

	int *ialoc;

	ialoc = new int [nlistloc+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	int j, ind1, k, kk, ibeg, ii1, ii2, nlist3loc;

	int ip1 = 0;
	int ip2 = 0;
	int iend1 = nlist-1;
	int iend2 = nlist2loc-1;

	nlistloc = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			ii1 = list[ip1];
			ii2 = plist2[ip2];
			if (ii1 == ii2) {
				_icycle1++;
				nlist3loc = 0;
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					jj = ja[j];
					_imask1[jj] = _icycle1;
					nlist3loc++;
				};
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					jj = pja2[j];
					if (_imask1[jj] != _icycle1) {
						_imask1[jj] = _icycle1;
						nlist3loc++;
					};
				};
				ialoc[nlistloc+1] = ialoc[nlistloc] + nlist3loc;
				nlistloc++;
				ip1++;
				ip2++;
			} else if (ii1 < ii2) {
				ialoc[nlistloc+1] = ialoc[nlistloc] + (ia[ip1+1]-ia[ip1]);
				nlistloc++;
				ip1++;
			} else if (ii1 > ii2) {
				ialoc[nlistloc+1] = ialoc[nlistloc] + (pia2[ip2+1]-pia2[ip2]);
				nlistloc++;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			ialoc[nlistloc+1] = ialoc[nlistloc] + (ia[ip1+1]-ia[ip1]);
			nlistloc++;
			ip1++;
		} else if (ip2 <= iend2) {
			ialoc[nlistloc+1] = ialoc[nlistloc] + (pia2[ip2+1]-pia2[ip2]);
			nlistloc++;
			ip2++;
		};
	};

	int nzjaloc = ialoc[nlistloc];

	int *jaloc;

	jaloc = new int [nzjaloc];
	if (!jaloc) MemoryFail (funcname);

	ip1 = 0;
	ip2 = 0;
	iend1 = nlist-1;
	iend2 = nlist2loc-1;

	nlistloc = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			ii1 = list[ip1];
			ii2 = plist2[ip2];
			if (ii1 == ii2) {
				_icycle1++;
				nlist3loc = 0;
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					jj = ja[j];
					_imask1[jj] = _icycle1;
					list2loc[nlist3loc] = jj;
					nlist3loc++;
				};
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					jj = pja2[j];
					if (_imask1[jj] != _icycle1) {
						_imask1[jj] = _icycle1;
						list2loc[nlist3loc] = jj;
						nlist3loc++;
					};
				};
				sort (list2loc, list2loc + nlist3loc);
				ibeg = ialoc[nlistloc];
				for (j=0;j<nlist3loc;j++) jaloc[ibeg+j] = list2loc[j];
				nlistloc++;
				ip1++;
				ip2++;
			} else if (ii1 < ii2) {
				ibeg = ialoc[nlistloc];
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					jaloc[ibeg+j-ia[ip1]] = ja[j];
				};
				nlistloc++;
				ip1++;
			} else if (ii1 > ii2) {
				ibeg = ialoc[nlistloc];
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					jaloc[ibeg+j-pia2[ip2]] = pja2[j];
				};
				nlistloc++;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			ibeg = ialoc[nlistloc];
			for (j=ia[ip1];j<ia[ip1+1];j++) {
				jaloc[ibeg+j-ia[ip1]] = ja[j];
			};
			nlistloc++;
			ip1++;
		} else if (ip2 <= iend2) {
			ibeg = ialoc[nlistloc];
			for (j=pia2[ip2];j<pia2[ip2+1];j++) {
				jaloc[ibeg+j-pia2[ip2]] = pja2[j];
			};
			nlistloc++;
			ip2++;
		};
	};

// Create new matrix and store new data

	CSMatrix temp (nlistloc,nzjaloc);

	for (i=0;i<nlistloc;i++) temp.list[i] = listloc[i];
	for (i=0;i<=nlistloc;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) temp.ja[i] = jaloc[i];

	temp.m      = m;
	temp.n      = n;
	temp.nsupc  = nsupc;
	temp.nsupr  = nsupr;
	temp.nlist  = nlistloc;
	temp.nzja   = nzjaloc;

// Add the block lists of each index

	int *iablkstr;

	iablkstr = new int [nzjaloc+1];
	if (!iablkstr) MemoryFail (funcname);

	for (i=0;i<=nzjaloc;i++) iablkstr[i] = 0;

	int *piablkstr = _ablkstr.GetIa ();
	int *pjablkstr = _ablkstr.GetJa ();
	int *piablkstradd = _ablkstradd.GetIa ();
	int *pjablkstradd = _ablkstradd.GetJa ();

	int ind;

	ip1 = 0;
	ip2 = 0;
	iend1 = nlist-1;
	iend2 = nlist2loc-1;

	nlistloc = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			ii1 = list[ip1];
			ii2 = plist2[ip2];
			if (ii1 == ii2) {
				for (j=ialoc[nlistloc];j<ialoc[nlistloc+1];j++) {
					jj = jaloc[j];
					list2loc[jj] = j;
					listloc[jj] = -1;
				};
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					jj = ja[j];
					ind = list2loc[jj];
					iablkstr[ind+1] = (piablkstr[j+1]-piablkstr[j]);
					listloc[jj] = j;
				};
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					jj = pja2[j];
					ind = list2loc[jj];
					ind1 = listloc[jj];
					if (ind1 >= 0) {
						_icycle1++;
						nlist3loc = 0;
						for (k=piablkstr[ind1];k<piablkstr[ind1+1];k++) {
							kk = pjablkstr[k];
							_imask1[kk] = _icycle1;
							nlist3loc++;
						};
						for (k=piablkstradd[j];k<piablkstradd[j+1];k++) {
							kk = pjablkstradd[k];
							if (_imask1[kk] != _icycle1) {
								_imask1[kk] = _icycle1;
								nlist3loc++;
							};
						};
						iablkstr[ind+1] = nlist3loc;
					} else {
						iablkstr[ind+1] = (piablkstradd[j+1]-piablkstradd[j]);
					};
				};
				nlistloc++;
				ip1++;
				ip2++;
			} else if (ii1 < ii2) {
				ibeg = ialoc[nlistloc];
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					iablkstr[ibeg+j-ia[ip1]+1] = (piablkstr[j+1]-piablkstr[j]);
				};
				nlistloc++;
				ip1++;
			} else if (ii1 > ii2) {
				ibeg = ialoc[nlistloc];
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					iablkstr[ibeg+j-pia2[ip2]+1] = (piablkstradd[j+1]-piablkstradd[j]);
				};
				nlistloc++;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			ibeg = ialoc[nlistloc];
			for (j=ia[ip1];j<ia[ip1+1];j++) {
				iablkstr[ibeg+j-ia[ip1]+1] = (piablkstr[j+1]-piablkstr[j]);
			};
			nlistloc++;
			ip1++;
		} else if (ip2 <= iend2) {
			ibeg = ialoc[nlistloc];
			for (j=pia2[ip2];j<pia2[ip2+1];j++) {
				iablkstr[ibeg+j-pia2[ip2]+1] = (piablkstradd[j+1]-piablkstradd[j]);
			};
			nlistloc++;
			ip2++;
		};
	};

	for (i=0;i<nzjaloc;i++) iablkstr[i+1] = iablkstr[i]+iablkstr[i+1];

	int nzjaadd = iablkstr[nzjaloc];

	int *jablkstr;

	jablkstr = new int [nzjaadd];
	if (!jablkstr) MemoryFail (funcname);

	int ibeg1;

	ip1 = 0;
	ip2 = 0;
	iend1 = nlist-1;
	iend2 = nlist2loc-1;

	nlistloc = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			ii1 = list[ip1];
			ii2 = plist2[ip2];
			if (ii1 == ii2) {
				for (j=ialoc[nlistloc];j<ialoc[nlistloc+1];j++) {
					jj = jaloc[j];
					list2loc[jj] = j;
					listloc[jj] = -1;
				};
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					jj = ja[j];
					ind = list2loc[jj];
					ibeg1 = iablkstr[ind];
					for (k=piablkstr[j];k<piablkstr[j+1];k++) {
						jablkstr[ibeg1+k-piablkstr[j]] = pjablkstr[k];
					};
					listloc[jj] = j;
				};
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					jj = pja2[j];
					ind = list2loc[jj];
					ind1 = listloc[jj];
					if (ind1 >= 0) {
						_icycle1++;
						nlist3loc = 0;
						for (k=piablkstr[ind1];k<piablkstr[ind1+1];k++) {
							kk = pjablkstr[k];
							_imask1[kk] = _icycle1;
							list3loc[nlist3loc] = kk;
							nlist3loc++;
						};
						for (k=piablkstradd[j];k<piablkstradd[j+1];k++) {
							kk = pjablkstradd[k];
							if (_imask1[kk] != _icycle1) {
								_imask1[kk] = _icycle1;
								list3loc[nlist3loc] = kk;
								nlist3loc++;
							};
						};
						sort (list3loc, list3loc + nlist3loc);
						ibeg1 = iablkstr[ind];
						for (k=0;k<nlist3loc;k++) jablkstr[ibeg1+k] = list3loc[k];
					} else {
						ibeg1 = iablkstr[ind];
						for (k=piablkstradd[j];k<piablkstradd[j+1];k++) {
							jablkstr[ibeg1+k-piablkstradd[j]] = pjablkstradd[k];
						};
					};
				};
				nlistloc++;
				ip1++;
				ip2++;
			} else if (ii1 < ii2) {
				ibeg = ialoc[nlistloc];
				for (j=ia[ip1];j<ia[ip1+1];j++) {
					ibeg1 = iablkstr[ibeg+j-ia[ip1]];
					for (k=piablkstr[j];k<piablkstr[j+1];k++) {
						jablkstr[ibeg1+k-piablkstr[j]] = pjablkstr[k];
					};
				};
				nlistloc++;
				ip1++;
			} else if (ii1 > ii2) {
				ibeg = ialoc[nlistloc];
				for (j=pia2[ip2];j<pia2[ip2+1];j++) {
					ibeg1 = iablkstr[ibeg+j-pia2[ip2]];
					for (k=piablkstradd[j];k<piablkstradd[j+1];k++) {
						jablkstr[ibeg1+k-piablkstradd[j]] = pjablkstradd[k];
					};
				};
				nlistloc++;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			ibeg = ialoc[nlistloc];
			for (j=ia[ip1];j<ia[ip1+1];j++) {
				ibeg1 = iablkstr[ibeg+j-ia[ip1]];
				for (k=piablkstr[j];k<piablkstr[j+1];k++) {
					jablkstr[ibeg1+k-piablkstr[j]] = pjablkstr[k];
				};
			};
			nlistloc++;
			ip1++;
		} else if (ip2 <= iend2) {
			ibeg = ialoc[nlistloc];
			for (j=pia2[ip2];j<pia2[ip2+1];j++) {
				ibeg1 = iablkstr[ibeg+j-pia2[ip2]];
				for (k=piablkstradd[j];k<piablkstradd[j+1];k++) {
					jablkstr[ibeg1+k-piablkstradd[j]] = pjablkstradd[k];
				};
			};
			nlistloc++;
			ip2++;
		};
	};

// Create new matrix and store new data

	CSMatrix tempstr (nzjaloc,nzjaadd);

	for (i=0;i<nzjaloc;i++) tempstr.list[i] = i;
	for (i=0;i<=nzjaloc;i++) tempstr.ia[i] = iablkstr[i];
	for (i=0;i<nzjaadd;i++) tempstr.ja[i] = jablkstr[i];

	tempstr.m      = nzjaloc;
	tempstr.n      = nzjaloc;
	tempstr.nsupc  = nzjaloc;
	tempstr.nsupr  = nzjaloc;
	tempstr.nlist  = nzjaloc;
	tempstr.nzja   = nzjaadd;

// Free work memory

	delete [] listloc;
	delete [] list1loc;
	delete [] list2loc;
	delete [] list3loc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] iablkstr;
	delete [] jablkstr;

	_amtradd = temp;
	_ablkstradd = tempstr;

};

// Author: Kharchenko S.A.
// CSMatrix: Filter matrix according to the indices
//========================================================================================
CSMatrix CSMatrix::FilterMatrix (int _isupbeg, int _isupend) const { // Filter matrix according to the indices

	const char *funcname = "FilterMatrix";

// Create new list of supernode rows

	int *listcount;

	int nlistloc = nlist;

	listcount = new int [nlistloc];
	if (!listcount) MemoryFail (funcname);

	int ilist;
	for (ilist=0;ilist<nlistloc;ilist++) {
		int isup = list[ilist];
		int nz = 0;
		if (isup >= _isupbeg && isup <= _isupend) {
			for (int j=ia[ilist];j<ia[ilist+1];j++) {
				int jj = ja[j];
				if (jj >= _isupbeg && jj <= _isupend) nz++;
			};
		};
		listcount[ilist] = nz;
	};

	int nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) nlistnew++;
	};

	int *listnew;

	listnew = new int [nlistnew];
	if (!listnew) MemoryFail (funcname);

	nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) {
			listnew[nlistnew] = list[ilist]-_isupbeg;
			nlistnew++;
		};
	};

// Create ia and ja arrays

	int *ialoc;

	ialoc = new int [nlistnew+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	nlistnew = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		if (listcount[ilist]>0) {
			ialoc[nlistnew+1] = ialoc[nlistnew] + listcount[ilist];
			nlistnew++;
		};
	};

	int nzjanew = ialoc[nlistnew];

	int *jaloc;

	jaloc = new int [nzjanew];
	if (!jaloc) MemoryFail (funcname);

	int nzjaloc = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		int isup = list[ilist];
		if (isup >= _isupbeg && isup <= _isupend) {
			for (int j=ia[ilist];j<ia[ilist+1];j++) {
				int jsup = ja[j];
				if (jsup >= _isupbeg && jsup <= _isupend) {
					jaloc[nzjaloc] = jsup-_isupbeg;
					nzjaloc++;
				};
			};
		};
	};

// Create new matrix

	CSMatrix temp (nlistnew,nzjanew);

// Store new data

	int i;
	for (i=0;i<nlistnew;i++) temp.list[i] = listnew[i];
	for (i=0;i<=nlistnew;i++) temp.ia[i] = ialoc[i];
	for (i=0;i<nzjaloc;i++) temp.ja[i] = jaloc[i];

	int nsuploc = _isupend-_isupbeg+1;

	temp.m      = nsuploc;
	temp.n      = nsuploc;
	temp.nsupc  = nsuploc;
	temp.nsupr  = nsuploc;
	temp.nlist  = nlistnew;
	temp.nzja   = nzjaloc;

// Free work memory

	delete [] listcount;
	delete [] listnew;
	delete [] ialoc;
	delete [] jaloc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute estimate of the Schur complement matrix
//========================================================================================
void CSMatrix::ComputeSchur (int _nbord, int _ncycle, CSMatrix &_ashur) const { // Compute estimate of the Schur complement matrix

	const char *funcname = "ComputeSchur";

// Compute the sparsity of the bordering

	int *imask, *listloc;
	int *ibord, *jbord;

	imask = new int [n+1];
	if (!imask) MemoryFail (funcname);
	listloc = new int [n+1];
	if (!listloc) MemoryFail (funcname);
	ibord = new int [n+1];
	if (!ibord) MemoryFail (funcname);

	int icycle = -1;
	int kii;
	for (kii=0;kii<n;kii++) imask[kii] = icycle;

	int ibegbrd = n - _nbord;
	int nloc = n;
	int nnloc = _nbord;

	ibord[0] = 0;

	int nzloc = 0;

	int i, j, jj;

	for (i=0;i<ibegbrd;i++) {
		icycle++;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
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

	for (i=0;i<ibegbrd;i++) {
		icycle++;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj >= ibegbrd) {
				if (imask[jj] != icycle) {
					jbord[nzloc] = jj;
					nzloc++;
					imask[jj] = icycle;
				};
			};
		};
	};

	for (i=0;i<ibegbrd;i++) {
		nzloc = ibord[i+1]-ibord[i];
		int ibs = ibord[i];
		qsort (jbord+ibs, nzloc, sizeof(int), compint);
	};

	for (kii=ibegbrd;kii<n;kii++) ibord[kii+1] = ibord[kii];

// Compute transposed bordering matrix

	nzloc = ibord[n];

	int *ibordt, *jbordt;

	ibordt = new int [n+1];
	if (!ibordt) MemoryFail (funcname);
	jbordt = new int [nzloc];
	if (!jbordt) MemoryFail (funcname);

	TransposeMatrix (n, ibord, jbord, ibordt, jbordt);

// Recompute aloc without bordering

	CSMatrix aloc;

	aloc = *this;

	nzloc = 0;

	for (i=0;i<ibegbrd;i++) {
		icycle++;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj < ibegbrd) {
				aloc.ja[nzloc] = jj;
				nzloc++;
			};
		};
		ibord[i+1] = nzloc;
	};

	for (i=ibegbrd;i<n;i++) {
		icycle++;
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			aloc.ja[nzloc] = jj;
			nzloc++;
		};
		ibord[i+1] = nzloc;
	};

	for (i=0;i<=n;i++) aloc.ia[i] = ibord[i];

// Compute extended bordering via ncycle parameter

	for (i=0;i<=_nbord;i++) ibordt[i] = ibordt[i+ibegbrd];

	int **jbordarr;

	jbordarr = new int * [_nbord];
	if (!jbordarr) MemoryFail (funcname);

	for (int icyc=0;icyc<_ncycle;icyc++) {

// Extend

		nzloc = 0;

		ibord[i] = 0;

		for (i=0;i<_nbord;i++) {

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
	for (i=0;i<=ibegbrd;i++) ibordt[i] = 0;

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

// Free work arrays

	delete [] imask;
	delete [] listloc;
	delete [] ibord;
	delete [] jbord;
	delete [] ibordt;
	delete [] jbordt;
	delete [] jbordarr;

	_ashur = shloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute estimate of the Schur complement matrix
//========================================================================================
void CSMatrix::ComputeSchur (int _nblks, int *_blks, // Compute estimate of the Schur complement matrix
										int _ibegbrd, int _iendbrd, int _ncycle,
										CSMatrix &_ashur, CSMatrix &_ashur_blk) const {

	const char *funcname = "ComputeSchur";

// Create the block numbers for all matrix elements

	int *ja2blk;

	ja2blk = new int [nzja];
	if(!ja2blk) MemoryFail(funcname);

	int i, jj;
	int iblk, iblkbeg, iblkend, iblkprev;

	iblkprev = 0;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
			ja2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = _nblks-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= _blks[iblk] && jj < _blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < _blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= _blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			ja2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *list2blk;

	list2blk = new int [nlist];
	if(!list2blk) MemoryFail(funcname);

	iblkprev = 0;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
			list2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = _nblks-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= _blks[iblk] && jj < _blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < _blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= _blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			list2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

// Create the list of blocks

	int *imaskblk;
	int *listblk;

	imaskblk = new int [_nblks];
	if(!imaskblk) MemoryFail(funcname);
	listblk = new int [_nblks];
	if(!listblk) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) imaskblk[i] = -1;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		imaskblk[iblk] = 1;
	};

	for (i=0;i<nzja;i++) {
		iblk = ja2blk[i];
		imaskblk[iblk] = 1;
	};

	int nlistblk = 0;

	for (i=0;i<_nblks;i++) {
		if (imaskblk[i] >= 0) {
			listblk[nlistblk] = i;
			imaskblk[i] = nlistblk;
			nlistblk++;
		};
	};

// Create the local mask and list arrays

	int *ibsmask;

	ibsmask = new int [_nblks];
	if (!ibsmask) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) ibsmask[i] = -1;

	int nloc = 0;
	int niloc;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		niloc = _blks[iblk+1]-_blks[iblk];
		nloc += niloc;
	};

	int *imask;
	int *indloc;
	int *ind2loc;
	int *ind2ia;

	imask = new int [nloc];
	if (!imask) MemoryFail (funcname);
	indloc = new int [nloc];
	if (!indloc) MemoryFail (funcname);
	ind2loc = new int [nloc];
	if (!ind2loc) MemoryFail (funcname);
	ind2ia = new int [nloc];
	if (!ind2ia) MemoryFail (funcname);

	for (i=0;i<nloc;i++) imask[i] = -1;
	for (i=0;i<nloc;i++) ind2ia[i] = -1;

	int jjblk, jjloc, ibs, iblkloc, j, irow, k;

	for (i=0;i<nlist;i++) {
		jjblk = list2blk[i];
		jj = list[i];
		jjloc = jj-_blks[jjblk];
		ibs = ibsmask[jjblk];
		imask[ibs+jjloc] = 1;
	};

	for (i=0;i<nzja;i++) {
		jjblk = ja2blk[i];
		jj = ja[i];
		jjloc = jj-_blks[jjblk];
		ibs = ibsmask[jjblk];
		imask[ibs+jjloc] = 1;
	};

	int nlistloc = 0;
	int ind = 0;

	for (iblkloc=0;iblkloc<nlistblk;iblkloc++) {
		iblk = listblk[iblkloc];
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (imask[ind] >=0) {
				indloc[nlistloc] = j;
				ind2loc[nlistloc] = iblk;
				imask[ind] = nlistloc;
				nlistloc++;
			};
			ind++;
		};
	};

// Split list data into the main block and boundary

	int nlist_main = 0;

	for (i=0;i<nlistloc;i++) {
		ind = indloc[i];
		if (ind < _ibegbrd) nlist_main++;
	};

	int nlist_brd = nlistloc-nlist_main;

// Compute transposed main block

	int *iat_main;
	int *iptr_main;

	iat_main = new int [nlist_main+1];
	if (!iat_main) MemoryFail (funcname);
	iptr_main = new int [nlist_main];
	if (!iptr_main) MemoryFail (funcname);

	for (i=0;i<=nlist_main;i++) {
		iat_main[i] = 0;
	};

	for (i=0;i<nlist;i++) {
		irow = list[i];
		if (irow < _ibegbrd) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj < _ibegbrd) {
					jjloc = jj-_blks[jjblk];
					ibs = ibsmask[jjblk];
					ind = imask[ibs+jjloc];
					iat_main[ind+1]++;
				};
			};
		};
	};

	for (i=0;i<nlist_main;i++) iat_main[i+1] = iat_main[i]+iat_main[i+1];
	for (i=0;i<nlist_main;i++) iptr_main[i] = iat_main[i];

	int nzja_main = iat_main[nlist_main];

	int *jat_main;
	int *jat2_main;

	jat_main = new int [nzja_main];
	if (!jat_main) MemoryFail (funcname);
	jat2_main = new int [nzja_main];
	if (!jat2_main) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		if (irow < _ibegbrd) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj < _ibegbrd) {
					jjloc = jj-_blks[jjblk];
					ibs = ibsmask[jjblk];
					ind = imask[ibs+jjloc];
					k = iptr_main[ind];
					jat_main[k] = irow;
					jat2_main[k] = iblk;
					iptr_main[ind]++;
				};
			};
		};
	};

// Compute the transposed sparsity of the bordering

	int *iat_brd;
	int *iptr_brd;

	iat_brd = new int [nlist_brd+1];
	if (!iat_brd) MemoryFail (funcname);
	iptr_brd = new int [nlist_brd];
	if (!iptr_brd) MemoryFail (funcname);

	for (i=0;i<=nlist_brd;i++) {
		iat_brd[i] = 0;
	};

	for (i=0;i<nlist;i++) {
		irow = list[i];
		if (irow < _ibegbrd) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj >= _ibegbrd) {
					jjloc = jj-_blks[jjblk];
					ibs = ibsmask[jjblk];
					ind = imask[ibs+jjloc]-nlist_main;
					iat_brd[ind+1]++;
				};
			};
		};
	};

	for (i=0;i<nlist_brd;i++) iat_brd[i+1] = iat_brd[i]+iat_brd[i+1];
	for (i=0;i<nlist_brd;i++) iptr_brd[i] = iat_brd[i];

	int nzja_brd = iat_brd[nlist_brd];

	int *jat_brd;
	int *jat2_brd;

	jat_brd = new int [nzja_brd];
	if (!jat_brd) MemoryFail (funcname);
	jat2_brd = new int [nzja_brd];
	if (!jat2_brd) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		if (irow < _ibegbrd) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj >= _ibegbrd) {
					jjloc = jj-_blks[jjblk];
					ibs = ibsmask[jjblk];
					ind = imask[ibs+jjloc]-nlist_main;
					k = iptr_brd[ind];
					jat_brd[k] = irow;
					jat2_brd[k] = iblk;
					iptr_brd[ind]++;
				};
			};
		};
	};

// Perform multiplication

	int *iat_brd_mvm;

	iat_brd_mvm = new int [nlist_brd+1];
	if (!iat_brd_mvm) MemoryFail (funcname);

	for (i=0;i<=nlist_brd;i++) {
		iat_brd_mvm[i] = 0;
	};

	int *imask_mvm;
	int *list_mvm;
	int *list2_mvm;

	imask_mvm = new int [nlistloc];
	if (!imask_mvm) MemoryFail (funcname);
	list_mvm = new int [nlistloc];
	if (!list_mvm) MemoryFail (funcname);
	list2_mvm = new int [nlistloc];
	if (!list2_mvm) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nlistloc;i++) imask_mvm[i] = icycle;

	int nlistmvm, indnew, kk, kkblk, kkloc, kbs, indk;

	for (i=0;i<nlist_brd;i++) {
		icycle++;
		nlistmvm = 0;
		for (j=iat_brd[i];j<iat_brd[i+1];j++) {
			jj = jat_brd[j];
			jjblk = jat2_brd[j];
			jjloc = jj-_blks[jjblk];
			ibs = ibsmask[jjblk];
			ind = ibs+jjloc;
			if (imask_mvm[ind] != icycle) {
				list_mvm[nlistmvm] = jj;
				list2_mvm[nlistmvm] = jjblk;
				imask_mvm[ind] = icycle;
				nlistmvm++;
			};
		};
		for (j=iat_brd[i];j<iat_brd[i+1];j++) {
			jj = jat_brd[j];
			jjblk = jat2_brd[j];
			jjloc = jj-_blks[jjblk];
			ibs = ibsmask[jjblk];
			ind = ibs+jjloc;
			indnew = imask[ind];
			for (k=iat_main[indnew];k<iat_main[indnew+1];k++) {
				kk = jat_main[k];
				kkblk = jat2_main[k];
				kkloc = kk-_blks[kkblk];
				kbs = ibsmask[kkblk];
				indk = kbs+kkloc;
				if (imask_mvm[indk] != icycle) {
					list_mvm[nlistmvm] = kk;
					list2_mvm[nlistmvm] = kkblk;
					imask_mvm[indk] = icycle;
					nlistmvm++;
				};
			};
		};
		iat_brd_mvm[i+1] = nlistmvm;
	};

	for (i=0;i<nlist_brd;i++) iat_brd_mvm[i+1] = iat_brd_mvm[i]+iat_brd_mvm[i+1];

	int nimax = 0;

	for (i=0;i<nlist_brd;i++) {
		niloc = iat_brd_mvm[i+1] - iat_brd_mvm[i];
		if (niloc > nimax) nimax = niloc;
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nimax];
	if (!iiarr) MemoryFail (funcname);

	int nzjat_brd_mvm = iat_brd_mvm[nlist_brd];

	int *jat_brd_mvm;
	int *jat2_brd_mvm;

	jat_brd_mvm = new int [nzjat_brd_mvm];
	if (!jat_brd_mvm) MemoryFail (funcname);
	jat2_brd_mvm = new int [nzjat_brd_mvm];
	if (!jat2_brd_mvm) MemoryFail (funcname);

	int ishift, irowloc, jbs, indj;

	for (i=0;i<nlist_brd;i++) {
		icycle++;
		nlistmvm = 0;
		for (j=iat_brd[i];j<iat_brd[i+1];j++) {
			jj = jat_brd[j];
			jjblk = jat2_brd[j];
			jjloc = jj-_blks[jjblk];
			ibs = ibsmask[jjblk];
			ind = ibs+jjloc;
			if (imask_mvm[ind] != icycle) {
				list_mvm[nlistmvm] = jj;
				list2_mvm[nlistmvm] = jjblk;
				imask_mvm[ind] = icycle;
				nlistmvm++;
			};
		};
		for (j=iat_brd[i];j<iat_brd[i+1];j++) {
			jj = jat_brd[j];
			jjblk = jat2_brd[j];
			jjloc = jj-_blks[jjblk];
			ibs = ibsmask[jjblk];
			ind = ibs+jjloc;
			indnew = imask[ind];
			for (k=iat_main[indnew];k<iat_main[indnew+1];k++) {
				kk = jat_main[k];
				kkblk = jat2_main[k];
				kkloc = kk-_blks[kkblk];
				kbs = ibsmask[kkblk];
				indk = kbs+kkloc;
				if (imask_mvm[indk] != icycle) {
					list_mvm[nlistmvm] = kk;
					list2_mvm[nlistmvm] = kkblk;
					imask_mvm[indk] = icycle;
					nlistmvm++;
				};
			};
		};
		for (j=0;j<nlistmvm;j++) {
			iiarr[j].intvalue = list_mvm[j];
			iiarr[j].int2value = list2_mvm[j];
		};
		qsort (iiarr, nlistmvm, sizeof(CIntInt), CIntInt::CompareIntInt);
		ishift = iat_brd_mvm[i];
		for (j=0;j<nlistmvm;j++) {
			jat_brd_mvm[ishift+j] = iiarr[j].intvalue;
			jat2_brd_mvm[ishift+j] = iiarr[j].int2value;
		};
	};

// Compute transposed mvm data

	int *ia_brd_mvm;
	int *iptr_brd_mvm;
	int *ja_brd_mvm;
	int *ja2_brd_mvm;

	ia_brd_mvm = new int [nlist_main+1];
	if (!ia_brd_mvm) MemoryFail (funcname);
	iptr_brd_mvm = new int [nlist_main];
	if (!iptr_brd_mvm) MemoryFail (funcname);
	ja_brd_mvm = new int [nzjat_brd_mvm];
	if (!ja_brd_mvm) MemoryFail (funcname);
	ja2_brd_mvm = new int [nzjat_brd_mvm];
	if (!ja2_brd_mvm) MemoryFail (funcname);

	for (i=0;i<=nlist_main;i++) ia_brd_mvm[i] = 0;

	for (i=0;i<nlist_brd;i++) {
		irow = indloc[nlist_main+i];
		iblk = ind2loc[nlist_main+i];
		for (j=iat_brd_mvm[i];j<iat_brd_mvm[i+1];j++) {
			jj = jat_brd_mvm[j];
			jjblk = jat2_brd_mvm[j];
			jjloc = jj-_blks[jjblk];
			ibs = ibsmask[jjblk];
			ind = imask[ibs+jjloc];
			ia_brd_mvm[ind+1]++;
		};
	};

	for (i=0;i<nlist_main;i++) ia_brd_mvm[i+1] = ia_brd_mvm[i]+ia_brd_mvm[i+1];

	for (i=0;i<nlist_main;i++) iptr_brd_mvm[i] = ia_brd_mvm[i];

	for (i=0;i<nlist_brd;i++) {
		irow = indloc[nlist_main+i];
		iblk = ind2loc[nlist_main+i];
		for (j=iat_brd_mvm[i];j<iat_brd_mvm[i+1];j++) {
			jj = jat_brd_mvm[j];
			jjblk = jat2_brd_mvm[j];
			jjloc = jj-_blks[jjblk];
			ibs = ibsmask[jjblk];
			ind = imask[ibs+jjloc];
			k = iptr_brd_mvm[i];
			ja_brd_mvm[k] = irow;
			ja2_brd_mvm[k] = iblk;
			iptr_brd_mvm[i]++;
		};
	};

// Compute the Schur complement

	if (_ncycle < 1 || _ncycle > 2) throw " CSMatrix::ComputeSchur: incorrect value of Ncycle parameter ";

	int *iat_ref;
	int *jat_ref;
	int *jat2_ref;

	if (_ncycle == 1) {
		iat_ref = iat_brd;
		jat_ref = jat_brd;
		jat2_ref = jat2_brd;
	} else {
		iat_ref = iat_brd_mvm;
		jat_ref = jat_brd_mvm;
		jat2_ref = jat2_brd_mvm;
	};

	int *iamvm;
	int *iptrmvm;

	iamvm = new int [nlistloc+1];
	if (!iamvm) MemoryFail (funcname);
	iptrmvm = new int [nlistloc];
	if (!iptrmvm) MemoryFail (funcname);

	for (i=0;i<=nlistloc;i++) iamvm[i] = 0;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		irowloc = irow - _blks[iblk];
		ibs = ibsmask[iblk];
		ind = imask[ibs+irowloc];
		if (irow >= _ibegbrd) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj >= _ibegbrd) {
					iamvm[ind+1]++;
				};
			};
		};
	};

	int index;

	for (i=0;i<nlist_brd;i++) {
		icycle++;
		irow = indloc[nlist_main+i];
		iblk = ind2loc[nlist_main+i];
		ibs = ibsmask[iblk];
		irowloc = irow-_blks[iblk];
		ind = imask[ibs+irowloc];
		nlistmvm = 0;
		index = ind2ia[ibs+irowloc];
		if (index >= 0) {
			for (j=ia[index];j<ia[index+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj >= _ibegbrd) {
					jbs = ibsmask[jjblk];
					jjloc = jj-_blks[jjblk];
					indj = imask[jbs+jjloc];
					imask_mvm[indj] = icycle;
				};
			};
		};
		for (j=iat_ref[i];j<iat_ref[i+1];j++) {
			jj = jat_ref[j];
			jjblk = jat2_ref[j];
			jbs = ibsmask[jjblk];
			jjloc = jj-_blks[jjblk];
			indj = imask[jbs+jjloc];
			for (k=ia_brd_mvm[indj];k<ia_brd_mvm[indj];k++) {
				kk = ja_brd_mvm[k];
				kkblk = ja2_brd_mvm[k];
				kbs = ibsmask[kkblk];
				kkloc = kk-_blks[kkblk];
				indk = imask[kbs+kkloc];
				if (imask_mvm[indk] != icycle) {
					imask_mvm[indk] = icycle;
					list_mvm[nlistmvm] = kk;
					list2_mvm[nlistmvm] = kkblk;
					nlistmvm++;
				};
			};
		};
		iamvm[ind+1] += nlistmvm;
	};

	for (i=0;i<nlistloc;i++) iamvm[i+1] = iamvm[i]+iamvm[i+1];
	for (i=0;i<nlistloc;i++) iptrmvm[i] = iamvm[i];

	int nzjamvm = iamvm[nlistloc];

	int *jamvm;
	int *ja2mvm;

	jamvm = new int [nzjamvm];
	if (!jamvm) MemoryFail (funcname);
	ja2mvm = new int [nzjamvm];
	if (!ja2mvm) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		irowloc = irow - _blks[iblk];
		ibs = ibsmask[iblk];
		ind = imask[ibs+irowloc];
		if (irow >= _ibegbrd) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj >= _ibegbrd) {
					k = iptrmvm[ind];
					jamvm[k] = jj;
					ja2mvm[k] = jjblk;
					iptrmvm[ind]++;
				};
			};
		};
	};

	for (i=0;i<nlist_brd;i++) {
		icycle++;
		irow = indloc[nlist_main+i];
		iblk = ind2loc[nlist_main+i];
		ibs = ibsmask[iblk];
		irowloc = irow-_blks[iblk];
		ind = imask[ibs+irowloc];
		nlistmvm = 0;
		index = ind2ia[ibs+irowloc];
		if (index >= 0) {
			for (j=ia[index];j<ia[index+1];j++) {
				jj = ja[j];
				jjblk = ja2blk[j];
				if (jj >= _ibegbrd) {
					jbs = ibsmask[jjblk];
					jjloc = jj-_blks[jjblk];
					indj = imask[jbs+jjloc];
					imask_mvm[indj] = icycle;
				};
			};
		};
		for (j=iat_ref[i];j<iat_ref[i+1];j++) {
			jj = jat_ref[j];
			jjblk = jat2_ref[j];
			jbs = ibsmask[jjblk];
			jjloc = jj-_blks[jjblk];
			indj = imask[jbs+jjloc];
			for (k=ia_brd_mvm[indj];k<ia_brd_mvm[indj];k++) {
				kk = ja_brd_mvm[k];
				kkblk = ja2_brd_mvm[k];
				kbs = ibsmask[kkblk];
				kkloc = kk-_blks[kkblk];
				indk = imask[kbs+kkloc];
				if (imask_mvm[indk] != icycle) {
					imask_mvm[indk] = icycle;
					list_mvm[nlistmvm] = kk;
					list2_mvm[nlistmvm] = kkblk;
					nlistmvm++;
				};
			};
		};
		k = iptrmvm[ind];
		for (j=0;j<nlistmvm;j++) {
			jamvm[k] = list_mvm[j];
			ja2mvm[k] = list2_mvm[j];
			k++;
		};
		iptrmvm[ind] += nlistmvm;
	};

// Store the Schur complement submatrix

	int nnloc = nlist_brd;
	int nzloc = iamvm[nlistloc];

	CSMatrix shloc (nnloc,nzloc);

	for (i=0;i< nnloc;i++) shloc.list[i] = indloc[nlist_main+i];
	for (i=0;i<=nnloc;i++) shloc.ia[i] = iamvm[nlist_main+i]-iamvm[nlist_main];
	for (i=0;i< nzloc;i++) shloc.ja[i] = jamvm[i];

	shloc.m = m;
	shloc.n = n;
	shloc.nsupr = nsupr;
	shloc.nsupc = nsupc;
	shloc.nzja = nzloc;
	shloc.nlist = nnloc;

	for (i=0;i<nnloc;i++) {
		nzloc = shloc.ia[i+1]-shloc.ia[i];
		int ibs = shloc.ia[i];
		qsort (shloc.ja+ibs, nzloc, sizeof(int), compint);
	};

// Perform filtering of the Schur complement matrix

	int *ianew;

	ianew = new int [nnloc+1];
	if (!ianew) MemoryFail (funcname);

	int nnnew = 0;
	int nznew = 0;

	ianew[0] = 0;

	for (i=0;i<nnloc;i++) {
		irow = shloc.list[i];
		if (irow <= _iendbrd) {
			nnnew++;
			for (j=shloc.ia[i];j<shloc.ia[i+1];j++) {
				jj = shloc.ja[j];
				if (jj <= _iendbrd) {
					shloc.ja[nznew] = jj;
					nznew++;
				};
			};
		};
		ianew[i+1] = nznew;
	};

	for (i=0;i<=nnnew;i++) shloc.ia[i] = ianew[i];

	shloc.nzja = nznew;
	shloc.nlist = nnnew;

// For each row of the Schur complement create the list of blocks that generate it

	int *imask_mvm_blk;
	int *list_mvm_blk;
	int *ia_mvm_blk;

	imask_mvm_blk = new int [_nblks];
	if (!imask_mvm_blk) MemoryFail (funcname);
	list_mvm_blk = new int [_nblks];
	if (!list_mvm_blk) MemoryFail (funcname);
	ia_mvm_blk = new int [nlist_brd+1];
	if (!ia_mvm_blk) MemoryFail (funcname);

	for (i=0;i<_nblks;i++) imask_mvm_blk[i] = icycle;
	for (i=0;i<=nlist_brd;i++) ia_mvm_blk[i] = 0;

	for (i=0;i<nlist_brd;i++) {
		icycle++;
		iblk = ind2loc[nlist_main+i];
		nlistmvm = 0;
		for (j=iat_ref[i];j<iat_ref[i+1];j++) {
			jjblk = jat2_ref[j];
			if (imask_mvm_blk[jjblk] != icycle) {
				imask_mvm_blk[jjblk] = icycle;
				list_mvm_blk[nlistmvm] = jjblk;
				nlistmvm++;
			};
		};
		if (_ncycle == 1) {
			for (j=iat_brd[i];j<iat_brd[i+1];j++) {
				jjblk = jat2_brd[j];
				if (imask_mvm_blk[jjblk] != icycle) {
					imask_mvm_blk[jjblk] = icycle;
					list_mvm_blk[nlistmvm] = jjblk;
					nlistmvm++;
				};
			};
		};
		ia_mvm_blk[i+1] = nlistmvm;
	};

	for (i=0;i<nlist_brd;i++) ia_mvm_blk[i+1] = ia_mvm_blk[i]+ia_mvm_blk[i+1];

	int nzja_blk = ia_mvm_blk[nlist_brd];

	int *ja_mvm_blk;

	ja_mvm_blk = new int [nzja_blk];
	if (!ja_mvm_blk) MemoryFail (funcname);

	for (i=0;i<nlist_brd;i++) {
		icycle++;
		iblk = ind2loc[nlist_main+i];
		nlistmvm = 0;
		for (j=iat_ref[i];j<iat_ref[i+1];j++) {
			jjblk = jat2_ref[j];
			if (imask_mvm_blk[jjblk] != icycle) {
				imask_mvm_blk[jjblk] = icycle;
				list_mvm_blk[nlistmvm] = jjblk;
				nlistmvm++;
			};
		};
		if (_ncycle == 1) {
			for (j=iat_brd[i];j<iat_brd[i+1];j++) {
				jjblk = jat2_brd[j];
				if (imask_mvm_blk[jjblk] != icycle) {
					imask_mvm_blk[jjblk] = icycle;
					list_mvm_blk[nlistmvm] = jjblk;
					nlistmvm++;
				};
			};
		};
		qsort (list_mvm_blk, nlistmvm, sizeof(int), compint);
		ibs = ia_mvm_blk[i];
		for (j=0;j<nlistmvm;j++) {
			ja_mvm_blk[ibs+j] = list_mvm_blk[j];
		};
	};

	nzja_blk = ia_mvm_blk[nnnew];

	CSMatrix shblk (nnnew,nzja_blk);

	for (i=0;i< nnnew;i++) shblk.list[i] = indloc[nlist_main+i];
	for (i=0;i<=nnnew;i++) shblk.ia[i] = ia_mvm_blk[i];
	for (i=0;i<nzja_blk;i++) shblk.ja[i] = ja_mvm_blk[i];

	shblk.m = m;
	shblk.n = n;
	shblk.nsupr = nsupr;
	shblk.nsupc = nsupc;
	shblk.nzja = nzja_blk;
	shblk.nlist = nnnew;

// Free work arrays

	delete [] ja2blk;
	delete [] list2blk;
	delete [] imaskblk;
	delete [] listblk;
	delete [] ibsmask;
	delete [] imask;
	delete [] indloc;
	delete [] ind2loc;
	delete [] ind2ia;
	delete [] iat_main;
	delete [] iptr_main;
	delete [] jat_main;
	delete [] jat2_main;
	delete [] iat_brd;
	delete [] iptr_brd;
	delete [] jat_brd;
	delete [] jat2_brd;
	delete [] iat_brd_mvm;
	delete [] imask_mvm;
	delete [] list_mvm;
	delete [] list2_mvm;
	delete [] iiarr;
	delete [] jat_brd_mvm;
	delete [] jat2_brd_mvm;
	delete [] ia_brd_mvm;
	delete [] iptr_brd_mvm;
	delete [] ja_brd_mvm;
	delete [] ja2_brd_mvm;
	delete [] iamvm;
	delete [] iptrmvm;
	delete [] jamvm;
	delete [] ja2mvm;
	delete [] ianew;
	delete [] imask_mvm_blk;
	delete [] list_mvm_blk;
	delete [] ia_mvm_blk;
	delete [] ja_mvm_blk;

// Return the result

	_ashur = shloc;
	_ashur_blk = shblk;

};

// Author: Kharchenko S.A.
// CSMatrix: Modify partitioning and ordering according to the list
//========================================================================================
void CSMatrix::ModifyPartitioning (int *_order, int _nblks, int *_blks, // Modify partitioning and ordering according to the list
											 int _nlist, int *_list) {

	const char *funcname = "ModifyPartitioning";

// Allocate the memory

	int nloc = _blks[_nblks];

	int *blksnw;
	int *iorder, *iordernew;

	blksnw = new int [_nblks+1];
	if (!blksnw) MemoryFail (funcname);
	iorder = new int [nloc];
	if (!iorder) MemoryFail (funcname);
	iordernew = new int [nloc];
	if (!iordernew) MemoryFail (funcname);

// Compute new ordering

	int i, j, jj, jjnew;

	for (i=0;i<nloc;i++) {
		iorder[_order[i]] = i;
	};

	for (i=0;i<_nlist;i++) {
		jj = _list[i];
		jjnew = _order[jj];
		iorder[jjnew] = -1;
	};

	for (i=0;i<=_nblks;i++) blksnw[i] = 0;

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			if (iorder[j] >= 0) blksnw[i+1]++;
		};
	};

	for (i=0;i<_nblks;i++) blksnw[i+1] = blksnw[i]+blksnw[i+1];

	int ipord, iplist, k;

	ipord = 0;
	iplist = 0;

	for (i=0;i<nloc;i++) {
		if (iorder[i] >= 0) {
			iordernew[ipord] = iorder[i];
			ipord++;
		} else {
			k = blksnw[_nblks];
			iordernew[k] = _list[iplist];
			iplist++;
			blksnw[_nblks]++;
		};
	};

	for (i=0;i<=_nblks;i++) _blks[i] = blksnw[i];

	for (i=0;i<nloc;i++) _order[iordernew[i]] = i;

// Free work arrays

	delete [] blksnw;
	delete [] iorder;
	delete [] iordernew;

};

// Author: Kharchenko S.A.
// Description: Sort list indices and column indices
// CSMatrix::SortListAndColumns()
//========================================================================================
void CSMatrix::SortListAndColumns () { // Sort list indices and column indices

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

	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzja];
	if (!jaloc) MemoryFail (funcname);

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
		qsort (jaloc+ibeg,ni,sizeof(int),compint);
	};

	for (i=0;i<nlist;i++) list[i] = listloc[i];
	for (i=0;i<=nlist;i++) ia[i] = ialoc[i];
	for (i=0;i<nzja;i++) ja[i] = jaloc[i];

	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;

	delete [] order;
	delete [] iorder;

};

// Author: Kharchenko S.A.
// Description: Sort list indices and column indices (2index version)
// CSMatrix::SortListAndColumns2Index()
//========================================================================================
void CSMatrix::SortListAndColumns2Index () { // Sort list indices and column indices (2index version)

	const char *funcname = "SortListAndColumns2Index";

// Reorder the rows

	CInd2Int *iiarr;

	iiarr = new CInd2Int [nlist];
	if (!iiarr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<nlist;i++) {
		iiarr[i].indx = list2[i];
		iiarr[i].indy = list[i];
		iiarr[i].intvalue = i;
	};

	qsort (iiarr,nlist,sizeof(CInd2Int),CInd2Int::CompareInd2Int);

	int *order, *iorder;

	order = new int [nlist];
	if (!order) MemoryFail (funcname);
	iorder = new int [nlist];
	if (!iorder) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		iorder[i] = iiarr[i].intvalue;
	};

	for (i=0;i<nlist;i++) {
		order[iorder[i]] = i;
	};

	delete [] iiarr;

	int *listloc;
	int *list2loc;
	int *ialoc;
	int *jaloc;
	int *ja2loc;

	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlist];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzja];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzja];
	if (!ja2loc) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		j = order[i];
		listloc[j] = list[i];
		list2loc[j] = list2[i];
		ialoc[j+1] = ia[i+1]-ia[i];
	};
	ialoc[0] = 0;

	for (i=0;i<nlist;i++) ialoc[i+1] = ialoc[i+1] + ialoc[i];

	int ni;
	int nimax = 0;
	for (i=0;i<nlist;i++) {
		ni = (ialoc[i+1] - ialoc[i]);
		if (ni > nimax) nimax = ni;
	};

	iiarr = new CInd2Int [nimax];
	if (!iiarr) MemoryFail (funcname);

	int ibeg;

	for (i=0;i<nlist;i++) {
		j = order[i];
		ibeg = ialoc[j];
		ni = ialoc[j+1]-ialoc[j];
		for (k=ia[i];k<ia[i+1];k++) {
			iiarr[k-ia[i]].indx = ja2[k];
			iiarr[k-ia[i]].indy = ja[k];
			iiarr[k-ia[i]].intvalue = k-ia[i];
		};
		qsort (iiarr,ni,sizeof(CInd2Int),CInd2Int::CompareInd2Int);
		for (k=ia[i];k<ia[i+1];k++) {
			jaloc[k-ia[i]+ialoc[j]] = iiarr[k-ia[i]].indy;
			ja2loc[k-ia[i]+ialoc[j]] = iiarr[k-ia[i]].indx;
		};
	};

	for (i=0;i<nlist;i++) list[i] = listloc[i];
	for (i=0;i<nlist;i++) list2[i] = list2loc[i];
	for (i=0;i<=nlist;i++) ia[i] = ialoc[i];
	for (i=0;i<nzja;i++) ja[i] = jaloc[i];
	for (i=0;i<nzja;i++) ja2[i] = ja2loc[i];

	delete [] listloc;
	delete [] list2loc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] ja2loc;

	delete [] order;
	delete [] iorder;

	delete [] iiarr;

};

// Author: Kharchenko S.A.
// Description: Compute ja2blk array with the previous row block history
// CSMatrix::ComputeJa2Blk()
//========================================================================================
void CSMatrix::ComputeJa2Blk (int _nblks, int *_blks, // Compute ja2blk array with the previous row block history
										int _nlist, int *_ia, int *_ja,
										int *_ja2blk) {

	const char *funcname = "ComputeJa2Blk";

// Allocate work arrays

	int *blklist1;
	int *blklist2;

	blklist1 = new int [_nblks+2];
	if (!blklist1) MemoryFail (funcname);
	blklist2 = new int [_nblks+2];
	if (!blklist2) MemoryFail (funcname);

	blklist1[0] = 0;
	blklist1[1] = _nblks;

	blklist2[0] = 0;

	int ipblklist = 0;

	int *blklist = blklist1;
	int *blklistnew = blklist2;

	int nblklist = 1;
	int nblklistnew = 1;

// Main search cycle

	int iblkprev = 0;

	int i, j, jj, iblk, ibegblk, iendblk, ipblk, iblknext;

	for (i=0;i<_nlist;i++) {
		ipblk = 0;
		iblkprev = blklist[ipblk];
		for (j=_ia[i];j<_ia[i+1];j++) {
			jj = _ja[j];
			if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
				iblk = iblkprev;
			} else {
				iblknext = blklist[ipblk+1]-1;
				while (jj >= _blks[iblknext+1] && ipblk < nblklist-1) {
					ipblk++;
					iblknext = blklist[ipblk+1]-1;
				};
				iblkprev = blklist[ipblk];
				if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
					iblk = iblkprev;
				} else {
					ibegblk = blklist[ipblk];
					iendblk = blklist[ipblk+1]-1;
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
				blklistnew[nblklistnew] = iblk;
				nblklistnew++;
			};
			_ja2blk[j] = iblk;
			iblkprev = iblk;
		};
		ipblklist++;
		ipblklist = ipblklist%2;
		if (ipblklist == 0) {
			blklist = blklist1;
			blklistnew = blklist2;
		} else {
			blklist = blklist2;
			blklistnew = blklist1;
		};
		nblklist = nblklistnew;
		blklist[nblklist] = _nblks;
		nblklistnew = 1;
		blklistnew[0] = 0;
	};

// Free work arrays

	delete [] blklist1;
	delete [] blklist2;

};

// Author: Kharchenko S.A.
// CSMatrix: Pack matrix data
//========================================================================================
void CSMatrix::PackMatrix (int &_length, char *&_obj) { // Pack matrix data

	const char* funcname = "PackMatrix";

// Count the total size of the object

	_length = 11 * sizeof(int) + (2*nlist+1+nlist2+nzja+nzja2+nzja3) * sizeof(int);

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

	pHead = (int *) pLoc;
	pLoc += 11 * sizeof(int);

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

	pHead[0]  = m;
	pHead[1]  = n;
	pHead[2]  = nsupr;
	pHead[3]  = nsupc;
	pHead[4]  = nlist;
	pHead[5]  = nlist2;
	pHead[6]  = nzja;
	pHead[7]  = nzja2;
	pHead[8]  = nzja3;

// Fill arrays

	int j;

	for(j = 0; j < nlist;   j++) plist[j] = list[j];
	for(j = 0; j < nlist2;  j++) plist2[j] = list2[j];
	for(j = 0; j < nlist+1; j++) pia[j] = ia[j];
	for(j = 0; j < nzja;    j++) pja[j] = ja[j];
	for(j = 0; j < nzja2;   j++) pja2[j] = ja2[j];
	for(j = 0; j < nzja3;   j++) pja3[j] = ja3[j];

};

// Author: Kharchenko S.A.
// CSMatrix: UnPack matrix data
//========================================================================================
void CSMatrix::UnPackMatrix (int _length, char *_obj) { // UnPack matrix data

	const char* funcname = "UnPackMatrix";

// Free previous arrays

	delete [] list;
	delete [] list2;
	delete [] ia;
	delete [] ja;
	delete [] ja2;
	delete [] ja3;

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

	pHead = (int *) pLoc;
	pLoc += 11 * sizeof(int);

	m       = pHead[0];
	n       = pHead[1];
	nsupr   = pHead[2];
	nsupc   = pHead[3];
	nlist   = pHead[4];
	nlist2  = pHead[5];
	nzja    = pHead[6];
	nzja2   = pHead[7];
	nzja3   = pHead[8];

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

	int j;

	for(j = 0; j < nlist;   j++) list[j] = plist[j];
	for(j = 0; j < nlist2;  j++) list2[j] = plist2[j];
	for(j = 0; j < nlist+1; j++) ia[j] = pia[j];
	for(j = 0; j < nzja;    j++) ja[j] = pja[j];
	for(j = 0; j < nzja2;   j++) ja2[j] = pja2[j];
	for(j = 0; j < nzja3;   j++) ja3[j] = pja3[j];

};

// Author: Kharchenko S.A.
// CSMatrix: PostScript collapsed sparsity structure
//========================================================================================
void CSMatrix::A2Ps (int _collap, const char *_fname, int _Nbl, const int *_Blks) const { // PostScript collapsed sparsity structure

	const char *funcname = "A2Ps_00";

// Compute new blocks partitioning

	int *blksnw;

	blksnw = new int [_Nbl+2];
	if (!blksnw) MemoryFail (funcname);

	int nsupmx = nsupc;
	if (nsupr > nsupmx) nsupmx = nsupr;

	if (_Nbl == 0) {
		blksnw[0] = 0;
		blksnw[1] = (nsupmx+_collap-1) / _collap;
	} else {
		blksnw[0] = 0;
		for (int iblk=0;iblk<_Nbl;iblk++) {
			int ni = _Blks[iblk+1]-_Blks[iblk];
			int niloc = (ni+_collap-1) / _collap;
			blksnw[iblk+1] = blksnw[iblk]+niloc;
		};
	};

// Compute nd2sp and sp2nd arrays

	int nnew;

	if (_Nbl == 0) {
		nnew = blksnw[_Nbl+1];
	} else {
		nnew = blksnw[_Nbl];
	};

	int *nd2sp, *sp2nd;
	int i, j;

	nd2sp = new int [nsupr];
	if (!nd2sp) MemoryFail (funcname);
	sp2nd = new int [nnew+1];
	if (!sp2nd) MemoryFail (funcname);

	if (_Nbl == 0) {
		sp2nd[0] = 0;
		for (i=0;i<nnew-1;i++) {
			sp2nd[i+1] = sp2nd[i] + _collap;
		}
		sp2nd[nnew] = nsupmx;
	} else {
		sp2nd[0] = 0;
		for (int iblk=0;iblk<_Nbl;iblk++) {
			for (i=blksnw[iblk];i<blksnw[iblk+1]-1;i++) {
				sp2nd[i+1] = sp2nd[i] + _collap;
			};
			sp2nd[blksnw[iblk+1]] = _Blks[iblk+1];
		};
	};

	for (i=0;i<nnew;i++) {
		for (j=sp2nd[i];j<sp2nd[i+1];j++) {
			nd2sp[j] = i;
		};
	};

// Count the total number of elements in the collapsed matrix

	int isup, jsup, jj;

	int icycle=0;

	int *imask, *lstloc;

	imask = new int [nnew];
	if (!imask) MemoryFail (funcname);
	lstloc = new int [nnew];
	if (!lstloc) MemoryFail (funcname);

	for (i=0;i<nnew;i++) imask[i] = icycle;

	int nz = 0;
	int nlstloc;

	for (isup=0;isup<nnew;isup++) {
		icycle++;
		nlstloc = 0;
		for (i=sp2nd[isup];i<sp2nd[isup+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
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

// Compute collapsed matrix

	CSMatrix temp (nnew,nz);

	temp.ia[0] = 0;
	for (isup=0;isup<nnew;isup++) {
		icycle++;
		nlstloc = 0;
		for (i=sp2nd[isup];i<sp2nd[isup+1];i++) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
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
	for (i=0;i<nnew;i++) temp.list[i] = i;

// Print collapsed sparsity structure

	temp.A2Ps (_fname, _Nbl, blksnw);

// Delete work arrays

	delete [] blksnw;
	delete [] sp2nd;
	delete [] nd2sp;
	delete [] imask;
	delete [] lstloc;

};

// Author: Kharchenko S.A.
// CSMatrix: PostScript sparsity structure
//========================================================================================
void CSMatrix::A2Ps (const char *_fname, int _Nbl, const int *_Blks) const { // PostScript sparsity structure

//	const char *funcname = "A2Ps_01";

// Open output file

	ofstream fout (_fname); 

	if (!fout.is_open()) {
		cout << " Error: File named " << _fname << " is not opened !!!" << endl;
		return;
	};

// Write the header information

	fout << "%!PS-Adobe-2.0" << endl;
	fout << "%%BoundingBox: 0 0 600 600" << endl;
	fout << "/m {moveto} def % x y" << endl;
	fout << "/l {lineto} def % x y" << endl;
	fout << "/s {stroke} def % x y" << endl;
	fout << "/n {newpath} def % x y" << endl;
	fout << "/c {closepath} def % x y" << endl;

// Parameters of the local window

	int nsupmx = nsupc;
	if (nsupr > nsupmx) nsupmx = nsupr;

	double s, s1;

	s1 = 50.0e0;
	s = 500.0e0 / ((double) nsupmx);

// Print the bounding window

	double dx, dy;

	fout << SetPw << 0.03 << " setlinewidth" << endl;

	fout << " n" << endl;
	dx = s1+s*0.5; dy = s1+s*0.5;
	Round (dx,dy);
	fout << SetPw << dx << SetPw << dy << " m" << endl;
	dx = s1+s*(nsupmx+0.5); dy = s1+s*0.5;
	Round (dx,dy);
	fout << SetPw << dx << SetPw << dy << " l" << endl;
	dx = s1+s*(nsupmx+0.5); dy = s1+s*(nsupmx+0.5);
	Round (dx,dy);
	fout << SetPw << dx << SetPw << dy << " l" << endl;
	dx = s1+s*0.5; dy = s1+s*(nsupmx+0.5);
	Round (dx,dy);
	fout << SetPw << dx << SetPw << dy << " l" << endl;
	dx = s1+s*0.5; dy = s1+s*0.5;
	Round (dx,dy);
	fout << SetPw << dx << SetPw << dy << " l s c" << endl;

// Print blocks partitioning

	int i, j, k;
	double x, y;

	if (_Nbl > 0) {
		for (i=0; i<_Nbl-1; i++) {
			x = _Blks[i+1]+0.5;
			y = n-(_Blks[i+1]+0.5)+1;
			fout << " n" << endl;
			dx = s1+s*0.5; dy = s1+s*y;
			Round (dx,dy);
			fout << SetPw << dx << SetPw << dy << " m" << endl;
			dx = s1+s*(nsupmx+0.5); dy = s1+s*y;
			Round (dx,dy);
			fout << SetPw << dx << SetPw << dy << " l s c" << endl;
			dx = s1+s*x; dy = s1+s*0.5;
			Round (dx,dy);
			fout << SetPw << dx << SetPw << dy << " m" << endl;
			dx = s1+s*x; dy = s1+s*(nsupmx+0.5);
			Round (dx,dy);
			fout << SetPw << dx << SetPw << dy << " l s c" << endl;
		};
	};

// Compute radius of the circle

	double r;

	r = (log10(1.0e0*nsupmx)+1.0) / 4.0 ;
	if (r > 1.0e0) r=1.0e0;

	r = s * r / 2.0;

	if (r < 0.03e0) r=0.03e0;

	fout << SetPw << 2*r << " setlinewidth" << endl;
	fout << " /d {moveto currentpoint " << SetPw << r;
	fout << " 0 360 arc fill} def % x y" << endl;

// Print the sparsity

	int i1, j1;

	fout << " n" << endl;
	for (int ilist=0;ilist<nlist;ilist++) {
		i = list[ilist];
		for (k=ia[ilist];k<ia[ilist+1];k++) {
			j = ja[k];
			i1 = i+1;
			j1 = j+1;
			x = j1;
			y = n-i1+1;
			dx = s1+s*x-r; dy = s1+s*y;
			Round (dx,dy);
			fout << SetPw << dx << SetPw << dy << " m" << endl;
			dx = s1+s*x+r; dy = s1+s*y;
			Round (dx,dy);
			fout << SetPw << dx << SetPw << dy << " l" << endl;
		};
	};

// Write the footer

	fout << " s c" << endl;
	fout << " showpage" << endl;

// Close output file

	fout.close();
};

// Author: Kharchenko S.A.
// CSMatrix: PostScript sparsity structure
//========================================================================================
void CSMatrix::PrintSparsity (int _collap, const char *_fname, int _n, int *_ia, int *_ja) { // PostScript print of collapsed sparsity structure

//	const char *funcname = "PrintSparsity";

// Prepare new matrix

	int nzjaloc = _ia[_n];

	CSMatrix aloc (_n, nzjaloc);

	int *plist, *pia, *pja;

	plist = aloc.list;
	pia = aloc.ia;
	pja = aloc.ja;

	int i;

	for (i=0;i<_n;i++) plist[i] = i;
	for (i=0;i<=_n;i++) pia[i] = _ia[i];
	for (i=0;i<nzjaloc;i++) pja[i] = _ja[i];

	aloc.m = _n;
	aloc.n = _n;
	aloc.nlist = _n;
	aloc.nsupc = _n;
	aloc.nsupr = _n;
	aloc.nzja = nzjaloc;

	aloc.A2Ps (_collap,_fname, 0, &i);

};

// Author: Kharchenko S.A.
// CSMatrix: Read sparsity from disk
//========================================================================================
CSMatrix CSMatrix::ReadSparsity (istream &_fin) { // Read sparsity from disk

	const char *funcname = "ReadSparsity";

	int nloc, nzjalc;
	char buffer[20]="";

	_fin >> nloc;
	cout << " N = " << nloc << endl;

	int *ialoc;

	ialoc = new int [nloc+1];
	if (!ialoc) MemoryFail (funcname);

	int i;
	for (i=0; i<=nloc; i++) {
		_fin >> ialoc[i];
	};

	nzjalc = ialoc[nloc];

	cout << " Nzjalc = " << nzjalc << endl;

	CSMatrix temp (nloc,nzjalc);

	for (i=0; i<=nloc; i++) {
		temp.ia[i] = ialoc[i];
	};

	_fin.getline(buffer,10);

	for (i=0; i<nzjalc; i++) {
		_fin >> temp.ja[i];
	};

	for (i=0;i<nloc;i++) {
		int ibs = temp.ia[i];
		int nz = temp.ia[i+1]-temp.ia[i];
		qsort (temp.ja+ibs, nz, sizeof(int), compint);
	};

	for (i=0;i<nloc;i++) temp.list[i] = i;

	temp.m      = nloc;
	temp.n      = nloc;
	temp.nsupr  = nloc;
	temp.nsupc  = nloc;
	temp.nlist  = nloc;
	temp.nzja   = nzjalc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrix: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CSMatrix &_a) { // Output matrix

	_stream << " CSMatrix:" << endl;
	_stream << " M = " << _a.m << " N = " << _a.n;
	_stream << " NsupR = " << _a.nsupr << " NsupC = " << _a.nsupc;
	_stream << " Nlist = " << _a.nlist << " Nlist2 = " << _a.nlist2 << " Nzja = " << _a.nzja << " Nzja2 = " << _a.nzja2 << " Nzja3 = " << _a.nzja3 << endl;

	OutArr (_stream, " List =", _a.nlist, _a.list);
	OutArr (_stream, " List2 =", _a.nlist2, _a.list2);
	OutArr (_stream, " Ia =", _a.nlist+1, _a.ia);
	OutArr (_stream, " Ja =", _a.ia[_a.nlist], _a.ja);
	OutArr (_stream, " Ja2 =", _a.nzja2, _a.ja2);
	OutArr (_stream, " Ja3 =", _a.nzja3, _a.ja3);

	return _stream;

};
