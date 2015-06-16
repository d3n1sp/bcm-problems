//------------------------------------------------------------------------------------------------
// File: smatrixr.cpp
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

#include "smatrix.h"
#include "gsmatrix.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixR: Memory allocation zero data constructor
//========================================================================================
CSMatrixR::CSMatrixR (): CSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixR_00";

	nza = 0;
	nzatot = 0;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = 0.0e0;
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Memory allocation zero data constructor
//========================================================================================
CSMatrixR::CSMatrixR (int _nlist, int _nzja): CSMatrix (_nlist,_nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixR_01";

	nza = _nzja;
	nzatot = nzja;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = 0.0e0;
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Memory allocation zero data constructor
//========================================================================================
CSMatrixR::CSMatrixR (int _nlist, int _nzja, int _nza): CSMatrix (_nlist,_nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixR_02";

	nza = _nza;
	nzatot = _nza;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = 0.0e0;
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Memory allocation zero data constructor
//========================================================================================
CSMatrixR::CSMatrixR (int _nlist, int _nlist2, int _nzja, int _nzja2, int _nza): CSMatrix (_nlist,_nlist2,_nzja,_nzja2) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixR_03";

	nza = _nza;
	nzatot = _nza;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = 0.0e0;
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Constructor for Laplace
//========================================================================================
CSMatrixR::CSMatrixR (int _nn): CSMatrix (_nn*_nn,3*_nn*_nn) { // Constructor for Laplace

	const char *funcname = "CSMatrixR_04";

	nza = nzja;

	a = new double [nza];
	if (!a) MemoryFail (funcname);

// Init ia, ja and a
	nzja = -1;
	ia[0] = 0;
	int idiag=1;

	for (int i=0;i<_nn;i++) {
		for (int j=0;j<_nn;j++) {
			ja[++nzja] = idiag-1;
			a[nzja] = 4.0e0;
			if (j<_nn-1) {
				ja[++nzja] = idiag;
				a[nzja] = -1.0e0;
			};
			if (i<_nn-1) {
				ja[++nzja] = idiag+_nn-1;
				a[nzja] = -1.0e0;
//				a[nzja] = 0.0e0;
			};
			ia[idiag++] = nzja+1;
		};
	};
	nzja++;

	nza = nzja;

	nzatot = nza;

};


// Author: Kharchenko S.A.
// CSMatrixR: Constructor for Laplace^2 with diagonal shift
//========================================================================================
CSMatrixR::CSMatrixR (int _nn, double _dshft): CSMatrix (_nn*_nn,7*_nn*_nn) { // Constructor for Laplace^2 with diagonal shift

	const char *funcname = "CSSMatrix_05";

	nza = nzja;

	a = new double [nzja];
	if (!a) MemoryFail (funcname);

	// Init ia, ja and a
	nzja = -1;
	ia[0] = 0;

	int idiag=0;

	for (int i=0;i<_nn;i++) {
		for (int j=0;j<_nn;j++) {
			ja[++nzja] = idiag;
			a[nzja] = 20.0e0+_dshft;
			if (j<_nn-1) {
				ja[++nzja] = idiag+1;
				a[nzja] = -8.0e0;
			};
			if (j<_nn-2) {
				ja[++nzja] = idiag+2;
				a[nzja] = 1.0e0;
			};
			if (i<_nn-1) {
				if (j>0) {
					ja[++nzja] = idiag+_nn-1;
					a[nzja] = 2.0e0;
				};
				ja[++nzja] = idiag+_nn;
				a[nzja] = -8.0e0;
				if (j<_nn-1) {
					ja[++nzja] = idiag+_nn+1;
					a[nzja] = 2.0e0;
				};
			};
			if (i<_nn-2) {
				ja[++nzja] = idiag+2*_nn;
				a[nzja] = 1.0e0;
			};
			ia[++idiag] = nzja+1;
		};
	};
	nzja++;

	nza = nzja;

	nzatot = nzja;

};

// Author: Kharchenko S.A.
// CSMatrixR: Laplace^2 structure unsymmetric matrix constructor with diagonal shift
//========================================================================================
CSMatrixR::CSMatrixR (int _nn, double _dshft, double _dunsy) { // Laplace^2 structure unsymmetric matrix constructor with diagonal shift

	const char *funcname = "CSMatrix_06";

	delete [] list;
	delete [] ia;
	delete [] ja;

	n = _nn*_nn;
	m = n;
	nlist = n;
	nzja = 14*n;
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) {
		list[i] = i;
	};
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	a = new double [nzja];
	if (!a) MemoryFail (funcname);

	// Init ia, ja and a
	nzja = -1;
	ia[0] = 0;
	int j,idiag=0;
	for (i=0;i<_nn;i++) {
		for (j=0;j<_nn;j++) {
			if (i>1) {
				ja[++nzja] = idiag-2*_nn;
				a[nzja] = 1.0e0-_dunsy;
			};
			if (i>0) {
				if (j>0) {
					ja[++nzja] = idiag-_nn-1;
					a[nzja] = 2.0e0-_dunsy;
				};
				ja[++nzja] = idiag-_nn;
				a[nzja] = -8.0e0-_dunsy;
				if (j<_nn-1) {
					ja[++nzja] = idiag-_nn+1;
					a[nzja] = 2.0e0-_dunsy;
				};
			};
			if (j>1) {
				ja[++nzja] = idiag-2;
				a[nzja] = 1.0e0-_dunsy;
			};
			if (j>0) {
				ja[++nzja] = idiag-1;
				a[nzja] = -8.0e0-_dunsy;
			};
			ja[++nzja] = idiag;
			a[nzja] = 20.0e0+_dshft;
			if (j<_nn-1) {
				ja[++nzja] = idiag+1;
				a[nzja] = -8.0e0+_dunsy;
			};
			if (j<_nn-2) {
				ja[++nzja] = idiag+2;
				a[nzja] = 1.0e0+_dunsy;
			};
			if (i<_nn-1) {
				if (j>0) {
					ja[++nzja] = idiag+_nn-1;
					a[nzja] = 2.0e0+_dunsy;
				};
				ja[++nzja] = idiag+_nn;
				a[nzja] = -8.0e0+_dunsy;
				if (j<_nn-1) {
					ja[++nzja] = idiag+_nn+1;
					a[nzja] = 2.0e0+_dunsy;
				};
			};
			if (i<_nn-2) {
				ja[++nzja] = idiag+2*_nn;
				a[nzja] = 1.0e0+_dunsy;
			};
			ia[++idiag] = nzja+1;
		};
	};
	nzja++;
	nza = nzja;
	nzatot = nza;
	nsupr = n;
	nsupc = n;
};

// Author: Kharchenko S.A.
// CSMatrixR: Laplace^3 structure unsymmetric matrix constructor with diagonal shift
//========================================================================================
CSMatrixR::CSMatrixR (int _nx, int _ny, int _nz, double _dshft, double _dunsy) { // Laplace^3 structure unsymmetric matrix constructor with diagonal shift

	const char *funcname = "CSMatrix_07";

	CSMatrix astr = CSMatrix::Sparsity3D (_nx, _ny, _nz);

	((CSMatrix *)this)->operator=(astr);

	a = new double [nzja];
	if (!a) MemoryFail (funcname);

	// Init a

	int i, j, jj;

	for (i=0;i<n;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (jj < i) {
				a[j] = -1.0e0-_dunsy;
			} else if (jj == i) {
				a[j] = 8.0e0+_dshft;
			} else if (jj > i) {
				a[j] = -1.0e0+_dunsy;
			};
		};
	};

	nza = nzja;
	nzatot = nza;
	nsupr = n;
	nsupc = n;

};

// Author: Kharchenko S.A.
// CSMatrixR: Laplace structure 3D unsymmetric matrix constructor with diagonal shift
//========================================================================================
CSMatrixR::CSMatrixR (int _nx, int _ny, int _nz, int _ibegz, int _iendz, double _dshft, double _dunsy) { // Laplace structure 3D unsymmetric matrix constructor with diagonal shift

	const char *funcname = "CSMatrix_07";

	CSMatrix astr = CSMatrix::Sparsity3D (_nx, _ny, _nz, _ibegz, _iendz);

	((CSMatrix *)this)->operator=(astr);

	a = new double [nzja];
	if (!a) MemoryFail (funcname);

	// Init a

	int i, j, jj, ilist;

	for (ilist=0;ilist<nlist;ilist++) {
		i = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj < i) {
				a[j] = -1.0e0-_dunsy;
			} else if (jj == i) {
				a[j] = 8.0e0+_dshft;
			} else if (jj > i) {
				a[j] = -1.0e0+_dunsy;
			};
		};
	};

	nza = nzja;
	nzatot = nza;
	nsupr = n;
	nsupc = n;

};

// Author: Kharchenko S.A.
// CSMatrixR: Highly ill-conditioned 3D unsymmetric matrix constructor
//========================================================================================
CSMatrixR::CSMatrixR (bool _ivar, int _nx, int _ny, int _nz, int _ibegz, int _iendz, double _eps, double _dunsy) { // Highly ill-conditioned 3D unsymmetric matrix constructor

	const char *funcname = "CSMatrix_08";

	CSMatrix astr = CSMatrix::Sparsity3D (_nx, _ny, _nz, _ibegz, _iendz);

	((CSMatrix *)this)->operator=(astr);

	a = new double [nzja];
	if (!a) MemoryFail (funcname);

// Init a

	double eps2 = _eps*_eps;

	double s[] = {-1.0e0,-_eps,-eps2-_dunsy/2.0e0,2.0e0*(1.0e0+_eps+eps2),-eps2+_dunsy/2.0e0,-_eps,-1.0e0};

	int is[] = {0,0,-1,0,1,0,0};
	int js[] = {0,-1,0,0,0,1,0};
	int ks[] = {-1,0,0,0,0,0,1};

	double d;

	int ii, jj, kk;
	int i, j, k;
	int iii, jjj, t, indd;
	int nxy = _nx*_ny;
	int nz = 0;

	for (kk=_ibegz;kk<=_iendz;kk++) {
		for (jj=0;jj<_ny;jj++) {
			for (ii=0;ii<_nx;ii++) {
				iii = kk*nxy+jj*_nx+ii;
				d = 0.0e0;
				for (t=0;t<7;t++) {
					i = ii + is[t];
					j = jj + js[t];
					k = kk + ks[t];
					jjj = k*nxy+j*_nx+i;
					if (iii == jjj) indd = nz;
					if (k<0) d += s[0];
					if (j<0) d += s[1];
					if (j>=_ny) d += s[5];
					if (k>=_nz) d += s[6];
					if (i>=0 && i<_nx && j>=0 && j<_ny && k>=0 && k<_nz) {
						a[nz] = s[t];
						nz++;
					};
				};
				a[indd] += d;
			};
		};
	};

	nza = nzja;
	nzatot = nza;
	nsupr = n;
	nsupc = n;

};

// Author: Kharchenko S.A.
// CSMatrixR: Unsymmetric banded matrix
//========================================================================================
CSMatrixR::CSMatrixR (char _ch, int _n, int _lb, int _rb) { // Unsymmetric banded matrix

	const char *funcname = "CSMatrixR_08";

	delete [] list;
	delete [] ia;
	delete [] ja;

	n = _n;
	m = n;
	nlist = n;
	nzja = n*(_lb+1+_rb);
	list = new int [nlist];
	if (!list) MemoryFail (funcname);
	int i;
	for (i=0;i<nlist;i++) list[i] = i;
	ia = new int [nlist+1];
	if (!ia) MemoryFail (funcname);
	ja = new int [nzja];
	if (!ja) MemoryFail (funcname);
	a = new double [nzja];
	if (!a) MemoryFail (funcname);

	// Init ia, ja and a

	int nz = 0;

	ia[0] = 0;

	int j,jloc;

	for (i=0;i<n;i++) {
		for (j=-_lb;j<=_rb;j++) {
			jloc = i+j;
			if (jloc >= 0 && jloc < n) {
				ja[nz] = jloc;
				a[nz] = 1.0e0 / ((double) (i-jloc) + 0.4e0);
				nz++;
			};
		};
		ia[i+1] = nz;
	};
	nza = nz;
	nzatot = nz;
	nsupr = n;
	nsupc = n;
};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate the memory
//========================================================================================
CSMatrixR::CSMatrixR (int _nlist, int _nlist2, int _nzja, int _nzja2) :CSMatrix (_nlist, _nlist2, _nzja, _nzja2) { // Unsymmetric banded matrix

	const char *funcname = "CSMatrixR_09";

	a = new double [nzja];
	if (!a) MemoryFail (funcname);

	nza = nzja;
	nzatot = nzja;

};

// Author: Kharchenko S.A.
// CSMatrixR: Copy constructor
//========================================================================================
CSMatrixR::CSMatrixR (const CSMatrixR &_aa) :CSMatrix (_aa) { // Copy constructor

	const char *funcname = "CSMatrixR_copy";

	nza = _aa.nza;
	nzatot = _aa.nzatot;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = _aa.a[i];
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Equality operator
//========================================================================================
CSMatrixR &CSMatrixR::operator= (const CSMatrixR &_aa) { // Equality operator

	const char *funcname = "CSMatrixR_=";

	delete [] a;

	this->CSMatrix::operator= ((const CSMatrix &) _aa);

	nza = _aa.nza;
	nzatot = _aa.nzatot;
	a = new double [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = _aa.a[i];
	};

	return *this;

};

// Author: Kharchenko S.A.
// CSMatrixR: Add two matrices
//========================================================================================
CSMatrixR CSMatrixR::operator+ (const CSMatrixR &_mtr2) const { // Add two matrices

//	const char *funcname = "CSMatrixR_+";

// Init dummy matrix

	CSMatrixR mtrdummy;

// Count the number of elements

	int nz = ia[nlist] + _mtr2.ia[nlist];

// Allocate and init temporary U

	CSMatrixR aloc (nlist,nz);

	aloc.m = n;
	aloc.n = n;

	int ilist, irow, ibeg1, iend1, ibeg2, iend2, ind1, ind2, jj1, jj2;

	for (ilist=0;ilist<nlist;ilist++) aloc.list[ilist] = list[ilist];

	aloc.ia[0] = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		ibeg1 = ia[ilist];
		iend1 = ia[ilist+1]-1;
		ibeg2 = _mtr2.ia[ilist];
		iend2 = _mtr2.ia[ilist+1]-1;
		ind1 = ibeg1;
		ind2 = ibeg2;
		while (ind1 <= iend1 || ind2 <= iend2) {
			if (ind1 <= iend1 && ind2 <= iend2) {
				jj1 = ja[ind1];
				jj2 = _mtr2.ja[ind2];
				if (jj1 == jj2) {
					aloc.ja[nz] = jj1;
					aloc.a[nz] = a[ind1];
					ind1++;
					ind2++;
					nz++;
				} else if (jj1 < jj2) {
					aloc.ja[nz] = jj1;
					aloc.a[nz] = a[ind1];
					ind1++;
					nz++;
				} else if (jj1 > jj2) {
					aloc.ja[nz] = jj2;
					aloc.a[nz] = _mtr2.a[ind2];
					ind2++;
					nz++;
				};
			} else if (ind1 <= iend1) {
				jj1 = ja[ind1];
				aloc.ja[nz] = jj1;
				aloc.a[nz] = a[ind1];
				ind1++;
				nz++;
			} else if (ind2 <= iend2) {
				jj2 = _mtr2.ja[ind2];
				aloc.ja[nz] = jj2;
				aloc.a[nz] = _mtr2.a[ind2];
				ind2++;
				nz++;
			};
		};
		aloc.ia[ilist+1] = nz;
	};

	aloc.nzja = nz;
	aloc.nza = nz;
	aloc.nzatot = nz;

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Add two matrices in super sparse two index format
//========================================================================================
CSMatrixR CSMatrixR::AddMatrSS2Index (bool _is_index3, const CSMatrixR &_mtr2) const { // Add two matrices in super sparse two index format

	const char *funcname = "AddMatricesSS2Index";

// Create the list of rows

	int nlistmax = nlist+_mtr2.nlist;

	int *listloc, *list2loc, *ialoc;

	listloc = new int [nlistmax];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlistmax];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlistmax+1];
	if (!ialoc) MemoryFail (funcname);

	int isup1, iblk1, isup2, iblk2;

	int nlistloc = 0;

	int iend1 = nlist-1;
	int iend2 = _mtr2.nlist-1;

	int ip1 = 0;
	int ip2 = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			isup1 = list[ip1];
			iblk1 = list2[ip1];
			isup2 = _mtr2.list[ip2];
			iblk2 = _mtr2.list2[ip2];
			if (iblk1 == iblk2) {
				if (isup1 == isup2) {
					listloc[nlistloc] = isup1;
					list2loc[nlistloc] = iblk1;
					ip1++;
					ip2++;
				} else if (isup1 < isup2) {
					listloc[nlistloc] = isup1;
					list2loc[nlistloc] = iblk1;
					ip1++;
				} else {
					listloc[nlistloc] = isup2;
					list2loc[nlistloc] = iblk2;
					ip2++;
				};
			} else if (iblk1 < iblk2) {
				listloc[nlistloc] = isup1;
				list2loc[nlistloc] = iblk1;
				ip1++;
			} else {
				listloc[nlistloc] = isup2;
				list2loc[nlistloc] = iblk2;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			isup1 = list[ip1];
			iblk1 = list2[ip1];
			listloc[nlistloc] = isup1;
			list2loc[nlistloc] = iblk1;
			ip1++;
		} else {
			isup2 = _mtr2.list[ip2];
			iblk2 = _mtr2.list2[ip2];
			listloc[nlistloc] = isup2;
			list2loc[nlistloc] = iblk2;
			ip2++;
		};
		nlistloc++;
	};

// Fill ja, ja2 (ja3) and a arrays

	int nzjamax = nzja + _mtr2.nzja;

	int *jaloc, *ja2loc, *ja3loc;
	double *aloc;

	jaloc = new int [nzjamax];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjamax];
	if (!ja2loc) MemoryFail (funcname);
	if (_is_index3) {
		ja3loc = new int [nzjamax];
		if (!ja3loc) MemoryFail (funcname);
	} else {
		ja3loc = new int [0];
		if (!ja3loc) MemoryFail (funcname);
	};
	aloc = new double [nzjamax];
	if (!aloc) MemoryFail (funcname);

	int ilist, irow, jj1, jj2, jp1, jp2;
	int iblk, jjblk1, jjblk2, irow1, irow2, jend1, jend2, nz, jlev, jlev1, jlev2;

	ialoc[0] = 0;
	nz = 0;

	ip1 = 0;
	ip2 = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		irow = listloc[ilist];
		iblk = list2loc[ilist];
		if (ip1 <= iend1) {
			irow1 = list[ip1];
			iblk1 = list2[ip1];
			if (irow1 != irow || iblk1 != iblk) iblk1 = -1;
		} else {
			iblk1 = -1;
		};
		if (ip2 <= iend2) {
			irow2 = _mtr2.list[ip2];
			iblk2 = _mtr2.list2[ip2];
			if (irow2 != irow || iblk2 != iblk) iblk2 = -1;
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
					jjblk1 = ja2[jp1];
					jj2 = _mtr2.ja[jp2];
					jjblk2 = _mtr2.ja2[jp2];
					if (jj1 == jj2 && jjblk1 == jjblk2) {
						jaloc[nz] = jj1;
						ja2loc[nz] = jjblk1;
						if (_is_index3) {
							jlev1 = ja3[jp1];
							jlev2 = _mtr2.ja3[jp2];
							jlev = jlev1;
							if (jlev > jlev2) jlev = jlev2;
							ja3loc[nz] = jlev;
						};
						aloc[nz] = a[jp1]+_mtr2.a[jp2];
						nz++;
						jp1++;
						jp2++;
					} else if (jjblk1 < jjblk2 || (jjblk1 == jjblk2 && jj1 < jj2)) {
						jaloc[nz] = jj1;
						ja2loc[nz] = jjblk1;
						if (_is_index3) {
							ja3loc[nz] = ja3[jp1];
						};
						aloc[nz] = a[jp1];
						nz++;
						jp1++;
					} else if (jjblk1 > jjblk2 || (jjblk1 == jjblk2 && jj1 > jj2)) {
						jaloc[nz] = jj2;
						ja2loc[nz] = jjblk2;
						if (_is_index3) {
							ja3loc[nz] = _mtr2.ja3[jp2];
						};
						aloc[nz] = _mtr2.a[jp2];
						nz++;
						jp2++;
					};
				} else if (jp1 <= jend1) {
					jaloc[nz] = ja[jp1];
					ja2loc[nz] = ja2[jp1];
					if (_is_index3) {
						ja3loc[nz] = ja3[jp1];
					};
					aloc[nz] = a[jp1];
					nz++;
					jp1++;
				} else if (jp2 <= jend2) {
					jaloc[nz] = _mtr2.ja[jp2];
					ja2loc[nz] = _mtr2.ja2[jp2];
					if (_is_index3) {
						ja3loc[nz] = _mtr2.ja3[jp2];
					};
					aloc[nz] = _mtr2.a[jp2];
					nz++;
					jp2++;
				};
			};
			ip1++;
			ip2++;
		} else if (iblk1 >= 0) {
			for (jp1=ia[ip1];jp1<ia[ip1+1];jp1++) {
				jaloc[nz] = ja[jp1];
				ja2loc[nz] = ja2[jp1];
				if (_is_index3) {
					ja3loc[nz] = ja3[jp1];
				};
				aloc[nz] = a[jp1];
				nz++;
			};
			ip1++;
		} else {
			for (jp2=_mtr2.ia[ip2];jp2<_mtr2.ia[ip2+1];jp2++) {
				jaloc[nz] = _mtr2.ja[jp2];
				ja2loc[nz] = _mtr2.ja2[jp2];
				if (_is_index3) {
					ja3loc[nz] = _mtr2.ja3[jp2];
				};
				aloc[nz] = _mtr2.a[jp2];
				nz++;
			};
			ip2++;
		};
		ialoc[ilist+1] = nz;
	};

// Store computed data

	CSMatrixR mtr;

	delete [] mtr.list;
	delete [] mtr.list2;
	delete [] mtr.ia;
	delete [] mtr.ja;
	delete [] mtr.ja2;
	delete [] mtr.ja3;
	delete [] mtr.a;

	mtr.list  = listloc;
	mtr.list2 = list2loc;
	mtr.ia    = ialoc;
	mtr.ja    = jaloc;
	mtr.ja2   = ja2loc;
	mtr.ja3   = ja3loc;
	mtr.a     = aloc;

	int ntot = _mtr2.m;

	mtr.m      = ntot;
	mtr.n      = ntot;
	mtr.nsupr  = ntot;
	mtr.nsupc  = ntot;
	mtr.nlist  = nlistloc;
	mtr.nlist2 = nlistloc;
	mtr.nzja   = nz;
	mtr.nzja2  = nz;
	if (_is_index3) {
		mtr.nzja3  = nz;
	} else {
		mtr.nzja3  = 0;
	};
	mtr.nza    = nz;
	mtr.nzatot = nz;

	return mtr;

};

// Author: Kharchenko S.A.
// CSMatrixR: Double add two matrices in super sparse two index format
//========================================================================================
void CSMatrixR::DoubleAddMatrSS2Index (bool _is_index3, // Double add two matrices in super sparse two index format
														const CSMatrixR &_mtrl, const CSMatrixR &_mtrl2, CSMatrixR &_mtrlsum,
														const CSMatrixR &_mtru, const CSMatrixR &_mtru2, CSMatrixR &_mtrusum) {

	const char *funcname = "DoubleAddMatricesSS2Index";

// Create the list of rows

	int nlistmax = _mtru.nlist+_mtru2.nlist;

	int *listloc, *list2loc, *ialoc;

	listloc = new int [nlistmax];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlistmax];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlistmax+1];
	if (!ialoc) MemoryFail (funcname);

	int isup1, iblk1, isup2, iblk2;

	int nlistloc = 0;

	int iend1 = _mtru.nlist-1;
	int iend2 = _mtru2.nlist-1;

	int ip1 = 0;
	int ip2 = 0;

	while (ip1 <= iend1 || ip2 <= iend2) {
		if (ip1 <= iend1 && ip2 <= iend2) {
			isup1 = _mtru.list[ip1];
			iblk1 = _mtru.list2[ip1];
			isup2 = _mtru2.list[ip2];
			iblk2 = _mtru2.list2[ip2];
			if (iblk1 == iblk2) {
				if (isup1 == isup2) {
					listloc[nlistloc] = isup1;
					list2loc[nlistloc] = iblk1;
					ip1++;
					ip2++;
				} else if (isup1 < isup2) {
					listloc[nlistloc] = isup1;
					list2loc[nlistloc] = iblk1;
					ip1++;
				} else {
					listloc[nlistloc] = isup2;
					list2loc[nlistloc] = iblk2;
					ip2++;
				};
			} else if (iblk1 < iblk2) {
				listloc[nlistloc] = isup1;
				list2loc[nlistloc] = iblk1;
				ip1++;
			} else {
				listloc[nlistloc] = isup2;
				list2loc[nlistloc] = iblk2;
				ip2++;
			};
		} else if (ip1 <= iend1) {
			isup1 = _mtru.list[ip1];
			iblk1 = _mtru.list2[ip1];
			listloc[nlistloc] = isup1;
			list2loc[nlistloc] = iblk1;
			ip1++;
		} else {
			isup2 = _mtru2.list[ip2];
			iblk2 = _mtru2.list2[ip2];
			listloc[nlistloc] = isup2;
			list2loc[nlistloc] = iblk2;
			ip2++;
		};
		nlistloc++;
	};

// Fill ja, ja2 (ja3) and a arrays

	int nzjamax = _mtru.nzja + _mtru2.nzja;

	int *jaloc, *ja2loc, *ja3loc;
	double *alloc;
	double *auloc;

	jaloc = new int [nzjamax];
	if (!jaloc) MemoryFail (funcname);
	ja2loc = new int [nzjamax];
	if (!ja2loc) MemoryFail (funcname);
	if (_is_index3) {
		ja3loc = new int [nzjamax];
		if (!ja3loc) MemoryFail (funcname);
	} else {
		ja3loc = new int [0];
		if (!ja3loc) MemoryFail (funcname);
	};
	alloc = new double [nzjamax];
	if (!alloc) MemoryFail (funcname);
	auloc = new double [nzjamax];
	if (!auloc) MemoryFail (funcname);

	int ilist, irow, jj1, jj2, jp1, jp2;
	int iblk, jjblk1, jjblk2, irow1, irow2, jend1, jend2, nz, jlev, jlev1, jlev2;

	ialoc[0] = 0;
	nz = 0;

	ip1 = 0;
	ip2 = 0;

	for (ilist=0;ilist<nlistloc;ilist++) {
		irow = listloc[ilist];
		iblk = list2loc[ilist];
		if (ip1 <= iend1) {
			irow1 = _mtru.list[ip1];
			iblk1 = _mtru.list2[ip1];
			if (irow1 != irow || iblk1 != iblk) iblk1 = -1;
		} else {
			iblk1 = -1;
		};
		if (ip2 <= iend2) {
			irow2 = _mtru2.list[ip2];
			iblk2 = _mtru2.list2[ip2];
			if (irow2 != irow || iblk2 != iblk) iblk2 = -1;
		} else {
			iblk2 = -1;
		};
		if (iblk1 >= 0 && iblk2 >= 0) {
			jend1 = _mtru.ia[ip1+1]-1;
			jend2 = _mtru2.ia[ip2+1]-1;
			jp1 = _mtru.ia[ip1];
			jp2 = _mtru2.ia[ip2];
			while (jp1 <= jend1 || jp2 <= jend2) {
				if (jp1 <= jend1 && jp2 <= jend2) {
					jj1 = _mtru.ja[jp1];
					jjblk1 = _mtru.ja2[jp1];
					jj2 = _mtru2.ja[jp2];
					jjblk2 = _mtru2.ja2[jp2];
					if (jj1 == jj2 && jjblk1 == jjblk2) {
						jaloc[nz] = jj1;
						ja2loc[nz] = jjblk1;
						if (_is_index3) {
							jlev1 = _mtru.ja3[jp1];
							jlev2 = _mtru2.ja3[jp2];
							jlev = jlev1;
							if (jlev > jlev2) jlev = jlev2;
							ja3loc[nz] = jlev;
						};
						alloc[nz] = _mtrl.a[jp1]+_mtrl2.a[jp2];
						auloc[nz] = _mtru.a[jp1]+_mtru2.a[jp2];
						nz++;
						jp1++;
						jp2++;
					} else if (jjblk1 < jjblk2 || (jjblk1 == jjblk2 && jj1 < jj2)) {
						jaloc[nz] = jj1;
						ja2loc[nz] = jjblk1;
						if (_is_index3) {
							ja3loc[nz] = _mtru.ja3[jp1];
						};
						alloc[nz] = _mtrl.a[jp1];
						auloc[nz] = _mtru.a[jp1];
						nz++;
						jp1++;
					} else if (jjblk1 > jjblk2 || (jjblk1 == jjblk2 && jj1 > jj2)) {
						jaloc[nz] = jj2;
						ja2loc[nz] = jjblk2;
						if (_is_index3) {
							ja3loc[nz] = _mtru2.ja3[jp2];
						};
						alloc[nz] = _mtrl2.a[jp2];
						auloc[nz] = _mtru2.a[jp2];
						nz++;
						jp2++;
					};
				} else if (jp1 <= jend1) {
					jaloc[nz] = _mtru.ja[jp1];
					ja2loc[nz] = _mtru.ja2[jp1];
					if (_is_index3) {
						ja3loc[nz] = _mtru.ja3[jp1];
					};
					alloc[nz] = _mtrl.a[jp1];
					auloc[nz] = _mtru.a[jp1];
					nz++;
					jp1++;
				} else if (jp2 <= jend2) {
					jaloc[nz] = _mtru2.ja[jp2];
					ja2loc[nz] = _mtru2.ja2[jp2];
					if (_is_index3) {
						ja3loc[nz] = _mtru2.ja3[jp2];
					};
					alloc[nz] = _mtrl2.a[jp2];
					auloc[nz] = _mtru2.a[jp2];
					nz++;
					jp2++;
				};
			};
			ip1++;
			ip2++;
		} else if (iblk1 >= 0) {
			for (jp1=_mtru.ia[ip1];jp1<_mtru.ia[ip1+1];jp1++) {
				jaloc[nz] = _mtru.ja[jp1];
				ja2loc[nz] = _mtru.ja2[jp1];
				if (_is_index3) {
					ja3loc[nz] = _mtru.ja3[jp1];
				};
				alloc[nz] = _mtrl.a[jp1];
				auloc[nz] = _mtru.a[jp1];
				nz++;
			};
			ip1++;
		} else {
			for (jp2=_mtru2.ia[ip2];jp2<_mtru2.ia[ip2+1];jp2++) {
				jaloc[nz] = _mtru2.ja[jp2];
				ja2loc[nz] = _mtru2.ja2[jp2];
				if (_is_index3) {
					ja3loc[nz] = _mtru2.ja3[jp2];
				};
				alloc[nz] = _mtrl2.a[jp2];
				auloc[nz] = _mtru2.a[jp2];
				nz++;
			};
			ip2++;
		};
		ialoc[ilist+1] = nz;
	};

// Store computed data

	delete [] _mtrusum.list;
	delete [] _mtrusum.list2;
	delete [] _mtrusum.ia;
	delete [] _mtrusum.ja;
	delete [] _mtrusum.ja2;
	delete [] _mtrusum.ja3;
	delete [] _mtrusum.a;

	_mtrusum.list  = listloc;
	_mtrusum.list2 = list2loc;
	_mtrusum.ia    = ialoc;
	_mtrusum.ja    = jaloc;
	_mtrusum.ja2   = ja2loc;
	_mtrusum.ja3   = ja3loc;
	_mtrusum.a     = auloc;

	int ntot = _mtru2.m;

	_mtrusum.m      = ntot;
	_mtrusum.n      = ntot;
	_mtrusum.nsupr  = ntot;
	_mtrusum.nsupc  = ntot;
	_mtrusum.nlist  = nlistloc;
	_mtrusum.nlist2 = nlistloc;
	_mtrusum.nzja   = nz;
	_mtrusum.nzja2  = nz;
	if (_is_index3) {
		_mtrusum.nzja3  = nz;
	} else {
		_mtrusum.nzja3  = 0;
	};
	_mtrusum.nza    = nz;
	_mtrusum.nzatot = nz;

	_mtrlsum = _mtrusum;

	int i;

	for (i=0;i<nz;i++) _mtrlsum.a[i] = alloc[i];

	delete [] alloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute L and U data for given matrix
//========================================================================================
void CSMatrixR::CombineLU (CSMatrixR &_mtrl, CSMatrixR &_mtru) const { // Compute L and U data for given matrix

//	const char *funcname = "CombineLU";

// Init dummy matrix

	CSMatrixR mtrdummy;

// Prepare U data

// Count the number of elements

	int nz, ilist, irow, j, jj;

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj >= irow) nz++;
		};
	};

// Allocate and init temporary U

	CSMatrixR mtru (nlist,nz);

	mtru.m = n;
	mtru.n = n;

	for (ilist=0;ilist<nlist;ilist++) mtru.list[ilist] = list[ilist];

	mtru.ia[0] = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj >= irow) {
				mtru.ja[nz] = jj;
				mtru.a [nz] = a[j];
				nz++;
			};
		};
		mtru.ia[ilist+1] = nz;
	};

// Store U data

	_mtru = mtru;

	mtru = mtrdummy;

// Prepare L data

// Count the number of elements

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj <= irow) nz++;
		};
	};

// Allocate and init temporary L

	CSMatrixR mtrl (nlist,nz);

	mtrl.m = n;
	mtrl.n = n;

	for (ilist=0;ilist<nlist;ilist++) mtrl.list[ilist] = list[ilist];

	mtrl.ia[0] = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj <= irow) {
				mtrl.ja[nz] = jj;
				mtrl.a [nz] = a[j];
				nz++;
			};
		};
		mtrl.ia[ilist+1] = nz;
	};

// Transpose L

	CSMatrixR mtrlt;

	mtrlt = mtrl.TranspMtr ();

// Copy the result

	_mtrl = mtrlt;

	mtrl = mtrdummy;
	mtrlt = mtrdummy;

// Extend the sparsity structures of L and U 

	_mtrl.ExtendSparsity (_mtru);
	_mtru.ExtendSparsity (_mtrl);

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute L and U data for given matrix in parallel
//========================================================================================
void CSMatrixR::CombineLU (CMPIComm &_comm, int _nblks, int *_blks, int *_blks2cpu, // Compute L and U data for given matrix in parallel
									CSMatrixR &_mtrl, CSMatrixR &_mtru) const {

	const char *funcname = "CombineLU";

	int nproc = _comm.GetNproc ();
//	int myid = _comm.GetMyid ();

// Init dummy matrix

	CSMatrixR mtrdummy;

// Prepare U data

// Count the number of elements

	int nz, ilist, irow, j, jj;

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj >= irow) nz++;
		};
	};

// Allocate and init temporary U

	CSMatrixR mtru (nlist,nz);

	mtru.m = n;
	mtru.n = n;

	for (ilist=0;ilist<nlist;ilist++) mtru.list[ilist] = list[ilist];

	mtru.ia[0] = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj >= irow) {
				mtru.ja[nz] = jj;
				mtru.a [nz] = a[j];
				nz++;
			};
		};
		mtru.ia[ilist+1] = nz;
	};

// Store U data

	_mtru = mtru;

	mtru = mtrdummy;

// Prepare local L data

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj <= irow) nz++;
		};
	};

// Allocate and init temporary L

	CSMatrixR mtrl (nlist,nz);

	mtrl.m = n;
	mtrl.n = n;

	for (ilist=0;ilist<nlist;ilist++) mtrl.list[ilist] = list[ilist];

	mtrl.ia[0] = 0;
	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		irow = list[ilist];
		for (j=ia[ilist];j<ia[ilist+1];j++) {
			jj = ja[j];
			if (jj <= irow) {
				mtrl.ja[nz] = jj;
				mtrl.a [nz] = a[j];
				nz++;
			};
		};
		mtrl.ia[ilist+1] = nz;
	};

// Mark block number for each sparsity element

	int *jablk;

	jablk = new int [nz];
	if (!jablk) MemoryFail (funcname);

	int iblkini = 0;
	int iblkprev = 0;

	int i, ip, iblkbeg, iblkend;

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
	int *iptrcpu;
	int *jacpuind;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);
	jacpuind = new int [nz];
	if (!jacpuind) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int iblk, iproc;

	for (j=0;j<nz;j++) {
		iblk = jablk[j];
		iproc = _blks2cpu[iblk];
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	int k;

	for (j=0;j<nz;j++) {
		iblk = jablk[j];
		iproc = _blks2cpu[iblk];
		k = iptrcpu[iproc];
		jacpuind[k] = j;
		iptrcpu[iproc]++;
	};

	int *ja2irow;

	ja2irow = new int [nz];
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
	jablkcpu = new int [_nblks];
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

	CSMatrixR *mtrarr;

	mtrarr = new CSMatrixR [nproc];
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

	int jjloc, ibs, nzloc, ipos;

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

		nzloc = iacpu[iproc+1]-iacpu[iproc];

		CSMatrixR temp (nlistloc,nzloc);

		for (i=0;i<nlistloc;i++) temp.list[i] = listloc[i];

		for (i=0;i<=nlistloc;i++) temp.ia[i] = 0;

		for (j=iacpu[iproc];j<iacpu[iproc+1];j++) {
			jind = jacpuind[j];
			jj = mtrl.ja[jind];
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
			temp.a[k] = mtrl.a[jind];
			iptr[ilist]++;
		};

		temp.m = n;
		temp.n = n;
		temp.nsupc = n;
		temp.nsupr = n;
		temp.nlist = nlistloc;
		temp.nzja = nzloc;
		temp.nza = nzloc;
		temp.nzatot = nzloc;

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
	if(info) throw " CSMatrixR::CombineLU: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

	if (NObjRecv != nproc) throw " CSMatrixR::CombineLU: wrong number of the received objects ";

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

	CGSMatrixR::CombineMatricesSort (n, nproc, listblkloc,
												mtrarr,
												_mtrl);

	for (i=0;i<nproc;i++) mtrarr[i] = mtrdummy;

	_mtrl.RemoveEquivalentRowsNumbers ();

// Free work arrays

	delete [] jablk;
	delete [] iacpu;
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
// CSMatrixR: Remove from the list of rows equivalent rows numbers
//========================================================================================
void CSMatrixR::RemoveEquivalentRowsNumbers () { // Remove from the list of rows equivalent rows numbers

	const char *funcname = "RemoveEquivalentRowsNumbers";

	if (nlist == 0) return;

	CSMatrixR mtrnew = *this;

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

	for (i=0;i<nlistnew;i++) {
		nz = mtrnew.ia[i+1]-mtrnew.ia[i];
		if (nz > nzmax) nzmax = nz;
	};

	CIntDouble *idarr;

	idarr = new CIntDouble [nzmax];
	if(!idarr) MemoryFail(funcname);

	int j;

	for (i=0;i<nlistnew;i++) {
		nz = 0;
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			idarr[nz].intvalue = mtrnew.ja[j];
			idarr[nz].dvalue = mtrnew.a[j];
			nz++;
		};
		qsort (idarr, nz, sizeof(CIntDouble), CIntDouble::CompareIntDouble);
		nz = 0;
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			mtrnew.ja[j] = idarr[nz].intvalue;
			mtrnew.a[j] = idarr[nz].dvalue;
			nz++;
		};
	};

	delete [] idarr;

	mtrnew.SetNlist (nlistnew);

// Return the result

	*this = mtrnew;

};

// Author: Kharchenko S.A.
// Description: Sort list indices and column indices
// CSMatrixR::SortListAndColumns()
//========================================================================================
void CSMatrixR::SortListAndColumns () { // Sort list indices and column indices

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
	double *aloc;

	listloc = new int [nlist];
	if (!listloc) MemoryFail (funcname);
	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);
	jaloc = new int [nzja];
	if (!jaloc) MemoryFail (funcname);
	aloc = new double [nzja];
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
			aloc[k-ia[i]+ialoc[j]] = a[k];
		};
	};

	int nzmax = 0;
	int nz;

	for (i=0;i<nlist;i++) {
		nz = ialoc[i+1]-ialoc[i];
		if (nz > nzmax) nzmax = nz;
	};

	CIntDouble *idarr;

	idarr = new CIntDouble [nzmax];
	if(!idarr) MemoryFail(funcname);

	for (i=0;i<nlist;i++) {
		nz = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			idarr[nz].intvalue = jaloc[j];
			idarr[nz].dvalue = aloc[j];
			nz++;
		};
		qsort (idarr, nz, sizeof(CIntDouble), CIntDouble::CompareIntDouble);
		nz = 0;
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jaloc[j] = idarr[nz].intvalue;
			aloc[j] = idarr[nz].dvalue;
			nz++;
		};
	};

	delete [] idarr;

	for (i=0;i<nlist;i++) list[i] = listloc[i];
	for (i=0;i<=nlist;i++) ia[i] = ialoc[i];
	for (i=0;i<nzja;i++) ja[i] = jaloc[i];
	for (i=0;i<nzja;i++) a[i] = aloc[i];

	delete [] listloc;
	delete [] ialoc;
	delete [] jaloc;
	delete [] aloc;

	delete [] order;
	delete [] iorder;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute new matrix as a sum of L and U data for given matrix
//========================================================================================
void CSMatrixR::AddLU (CSMatrixR &_mtrl, CSMatrixR &_mtru) { // Compute new matrix as a sum of L and U data for given matrix

//	const char *funcname = "AddLU";

// Transpose L matrix

	CSMatrixR mtrlt;

	*this = mtrlt;

	mtrlt = _mtrl.TranspMtr ();

	*this = mtrlt + _mtru;

};

// Author: Kharchenko S.A.
// CSMatrixR: Symmetrize the sparsity of a general rectangular matrix, replace added alements by zeroes
//========================================================================================
CSMatrixR CSMatrixR::SymmMtrRect () const { // Symmetrize the sparsity of a general rectangular matrix, replace added alements by zeroes

//	const char *funcname = "SymmMtrRect";

// Symmetrize the sparsity only

	CSMatrix asstr;

	asstr = ((CSMatrix *)this)->SymmMtrRect ();

// Create and fill elements data

	int nliststr = asstr.GetNlist ();
	int nzjastr = asstr.GetNzja ();

	CSMatrixR aloc (nliststr,nzjastr);

	((CSMatrix &)aloc) = asstr;

	int *plist = aloc.GetList ();
	int *pia = aloc.GetIa ();
	int *pja = aloc.GetJa ();
	double *pa = aloc.GetA ();

	int i, jj;

	for (i=0;i<nzjastr;i++) pa[i] = 0.0e0;

	int ilist1 = 0;
	int ilist2 = 0;

	int nz = 0;
	int nzloc, irow, irow2, ip, ip2, jj2;

	while (ilist1 < nliststr) {
		irow = plist[ilist1];
		nzloc = pia[ilist1+1]-pia[ilist1];
		if (ilist2 < nlist) {
			irow2 = list[ilist2];
			if (irow == irow2) {
				ip = pia[ilist1];
				ip2 = ia[ilist2];
				while (ip < pia[ilist1+1]) {
					jj = pja[ip];
					if (ip2 < ia[ilist2+1]) {
						jj2 = ja[ip2];
						if (jj == jj2) {
							pa[ip] = a[ip2];
							ip++;
							nz++;
							ip2++;
						} else if (jj < jj2) {
							ip++;
							nz++;
						} else {
							throw " CSMatrixR::SymmRectMtr: extra column index found ";
						};
					} else {
						ip++;
						nz++;
					};
				};
				ilist1++;
				ilist2++;
			} else if (irow < irow2) {
				ilist1++;
				nz += nzloc;
			} else {
				throw " CSMatrixR::SymmRectMtr: extra row found ";
			};
		} else {
			ilist1++;
			nz += nzloc;
		};
	};

	if (ilist2 != nlist) throw " CSMatrixR::SymmRectMtr: not all rows processed ";

	aloc.SetNza (nzjastr);
	aloc.SetNzatot (nzjastr);

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Symmetrize the sparsity of a general rectangular matrix, replace added alements by zeroes
//========================================================================================
CSMatrixR CSMatrixR::SymmMtrRect2Index () const { // Symmetrize the sparsity of a general rectangular matrix, replace added alements by zeroes

//	const char *funcname = "SymmMtrRect2Index";

// Symmetrize the sparsity only

	CSMatrix asstr;

	asstr = ((CSMatrix *)this)->SymmMtrRect2Index ();

// Create and fill elements data

	int nliststr = asstr.GetNlist ();
	int nzjastr = asstr.GetNzja ();

	CSMatrixR aloc (nliststr,nliststr,nzjastr,nzjastr);

	((CSMatrix &)aloc) = asstr;

	int *plist = aloc.GetList ();
	int *plist2 = aloc.GetList2 ();
	int *pia = aloc.GetIa ();
	int *pja = aloc.GetJa ();
	int *pja2 = aloc.GetJa2 ();
	double *pa = aloc.GetA ();

	int i, jj;

	for (i=0;i<nzjastr;i++) pa[i] = 0.0e0;

	int ilist1 = 0;
	int ilist2 = 0;

	int nz = 0;
	int nzloc, irow, irow2, iblk, iblk2, ip, ip2, jj2, jjblk, jjblk2;

	while (ilist1 < nliststr) {
		irow = plist[ilist1];
		iblk = plist2[ilist1];
		nzloc = pia[ilist1+1]-pia[ilist1];
		if (ilist2 < nlist) {
			irow2 = list[ilist2];
			iblk2 = list2[ilist2];
			if (iblk == iblk2 && irow == irow2) {
				ip = pia[ilist1];
				ip2 = ia[ilist2];
				while (ip < pia[ilist1+1]) {
					jj = pja[ip];
					jjblk = pja2[ip];
					if (ip2 < ia[ilist2+1]) {
						jj2 = ja[ip2];
						jjblk2 = ja2[ip2];
						if (jjblk == jjblk2 && jj == jj2) {
							pa[ip] = a[ip2];
							ip++;
							nz++;
							ip2++;
						} else if (jjblk < jjblk2 || (jjblk == jjblk2 && jj < jj2)) {
							ip++;
							nz++;
						} else {
							throw " CSMatrixR::SymmRectMtr2Index: extra column index found ";
						};
					} else {
						ip++;
						nz++;
					};
				};
				ilist1++;
				ilist2++;
			} else if (iblk < iblk2 || (iblk == iblk2 && irow < irow2)) {
				ilist1++;
				nz += nzloc;
			} else {
				throw " CSMatrixR::SymmRectMtr2Index: extra row found ";
			};
		} else {
			ilist1++;
			nz += nzloc;
		};
	};

	if (ilist2 != nlist) throw " CSMatrixR::SymmRectMtr2Index: not all rows processed ";

	aloc.SetNza (nzjastr);
	aloc.SetNzatot (nzjastr);

	return aloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Extend sparsity of the matrix according to the second matrix
//========================================================================================
void CSMatrixR::ExtendSparsity (const CSMatrixR &_mtr2) { // Extend sparsity of the matrix according to the second matrix

	const char *funcname = "ExtendSparsity";

// Check the data on entry

	if (nlist != _mtr2.nlist) throw " Error in sparsity extend: the numbers of nonzero rows do not coincide";

	int ilist;
	for (ilist=0;ilist<nlist;ilist++) {
		if (list[ilist] != _mtr2.list[ilist]) throw " Error in sparsity extend: the lists of nonzero rows do not coincide";
	};

	bool sparsity_coincide = true;

	int i;

	if (nzja != _mtr2.nzja) {
		sparsity_coincide = false;
	} else {
		for (i=0;i<=nlist;i++) {
			if (ia[i] != _mtr2.ia[i]) sparsity_coincide = false;
		};
		if (sparsity_coincide) {
			for (i=0;i<nzja;i++) {
				if (ja[i] != _mtr2.ja[i]) sparsity_coincide = false;
			};
		};
	};

	if (sparsity_coincide) return;

// Count the new extended number of elements in each row

	int *nelem;

	nelem = new int [nlist];
	if (!nelem) MemoryFail (funcname);

	int nz, ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;

	for (ilist=0;ilist<nlist;ilist++) {

// Scan both rows at the same time

		nz = 0;

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
					nz++;
					jind1++;
				} else if (jj1 > jj2) {
					nz++;
					jind2++;
				} else if (jj1 == jj2) {
					nz++;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				nz++;
				jind1++;
			} else if (jind2 <= iend2) {
				nz++;
				jind2++;
			};
		};

		nelem[ilist] = nz;

	};

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) nz += nelem[ilist];

// Allocate work memory for extended matrix

	CSMatrixR temp (nlist,nz);

	temp.m = m;
	temp.n = n;

// Init list and ia arrays

	temp.ia[0] = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		temp.list[ilist] = list[ilist];
		temp.ia[ilist+1] = temp.ia[ilist]+nelem[ilist];
	};

// Fill ja and a arrays

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {

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
					temp.ja[nz] = jj1;
					temp.a [nz] = a[jind1];
					nz++;
					jind1++;
				} else if (jj1 > jj2) {
					temp.ja[nz] = jj2;
					temp.a [nz] = 0.0e0;
					nz++;
					jind2++;
				} else if (jj1 == jj2) {
					temp.ja[nz] = jj1;
					temp.a [nz] = a[jind1];
					nz++;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				temp.ja[nz] = ja[jind1];
				temp.a [nz] = a[jind1];
				nz++;
				jind1++;
			} else if (jind2 <= iend2) {
				temp.ja[nz] = _mtr2.ja[jind2];
				temp.a [nz] = 0.0e0;
				nz++;
				jind2++;
			};
		};
	};

// Delete work arrays

	delete [] nelem;

// Return obtained matrix

	*this = temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Copy matrix into extended sparsity
//========================================================================================
void CSMatrixR::ExtendedCopy2Index (CSMatrix &_mtrsp, const CSMatrixR &_mtr2) { // Copy matrix into extended sparsity

//	const char *funcname = "ExtendedCopy2Index";

// Check the data on entry

	if (_mtrsp.nlist < _mtr2.nlist) {
		cout << " Wait " << endl;
		throw " CSMatrixR::ExtendedCopy2Index: too small number of nonzero rows in the sparsity";
	};

// Allocate work memory for extended matrix

	int nlistloc = _mtrsp.GetNlist ();
	int nzjaloc = _mtrsp.GetNzja ();

	CSMatrixR temp (nlistloc, nlistloc, nzjaloc, nzjaloc, nzjaloc);

	((CSMatrix &)temp) = _mtrsp;

// Assign the elements

	int i;

	for (i=0;i<nzjaloc;i++) temp.a[i] = 0.0e0;

	int ilist1, irow1, iblk1, jind1, iend1, jj1, jjblk1;
	int ilist2, irow2, iblk2, jind2, iend2, jj2, jjblk2;

	ilist1 = 0;
	ilist2 = 0;

	while (ilist1 < _mtrsp.nlist || ilist2 < _mtr2.nlist) {
		if (ilist1 < _mtrsp.nlist && ilist2 < _mtr2.nlist) {
			irow1 = _mtrsp.list[ilist1];
			iblk1 = _mtrsp.list2[ilist1];
			irow2 = _mtr2.list[ilist2];
			iblk2 = _mtr2.list2[ilist2];
			if (irow1 == irow2 && iblk1 == iblk2) {
				jind1 = _mtrsp.ia[ilist1]; 
				jind2 = _mtr2.ia[ilist2]; 
				iend1 = _mtrsp.ia[ilist1+1]-1; 
				iend2 = _mtr2.ia[ilist2+1]-1; 
				while (jind1 <= iend1 || jind2 <= iend2) {
					if (jind1 <= iend1 && jind2 <= iend2) {
						jj1 = _mtrsp.ja[jind1];
						jjblk1 = _mtrsp.ja2[jind1];
						jj2 = _mtr2.ja[jind2];
						jjblk2 = _mtr2.ja2[jind2];
						if (jj1 == jj2 && jjblk1 == jjblk2) {
							temp.a[jind1] = _mtr2.a[jind2];
							jind1++;
							jind2++;
						} else if (jjblk1 < jjblk2 || (jj1 < jj2 && jjblk1 == jjblk2)) {
							jind1++;
						} else {
							cout << " Wait " << endl;
							throw " CSMatrixR::ExtendedCopy2Index: smaller matrix have extra columns";;
						};
					} else if (jind1 <= iend1) {
						jind1++;
					} else if (jind2 <= iend2) {
						cout << " Wait " << endl;
						throw " CSMatrixR::ExtendedCopy2Index: smaller matrix have extra columns";
					};
				};
				ilist1++;
				ilist2++;
			} else if (iblk1 < iblk2 || (iblk1 == iblk2 && irow1 < irow2)) {
				ilist1++;
			} else {
				cout << " Wait " << endl;
				throw " CSMatrixR::ExtendedCopy2Index: smaller matrix have extra rows";
			};
		} else if (ilist1 < _mtrsp.nlist) {
			ilist1++;
		} else {
			cout << " Wait " << endl;
			throw " CSMatrixR::ExtendedCopy2Index: smaller matrix have extra rows";;
		};
	};

	temp.nza = nzjaloc;
	temp.nzatot = nzjaloc;

// Return obtained matrix

	*this = temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Cut the sparsity of the matrix according to the prescribed one
//========================================================================================
CSMatrixR CSMatrixR::CutSparsity (const CSMatrix &_mtr2) { // Cut the sparsity of the matrix according to the prescribed one

	const char *funcname = "CutSparsity";

// Check the data on entry

	int nlistloc = _mtr2.GetNlist ();

	if (nlist != nlistloc) throw " Error in sparsity cut: the numbers of nonzero rows do not coincide";

	int ilist;
	int *plistloc = _mtr2.GetList ();

	for (ilist=0;ilist<nlist;ilist++) {
		if (list[ilist] != plistloc[ilist]) throw " Error in sparsity cut: the lists of nonzero rows do not coincide";
	};

// Count the new extended number of elements in each row

	int *nelem;

	nelem = new int [nlist];
	if (!nelem) MemoryFail (funcname);

	int nz, ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;

	int *pialoc = _mtr2.GetIa ();
	int *pjaloc = _mtr2.GetJa ();

	for (ilist=0;ilist<nlist;ilist++) {

// Scan both rows at the same time

		nz = 0;

		ibeg1 = ia[ilist]; 
		iend1 = ia[ilist+1]-1; 
		ibeg2 = pialoc[ilist]; 
		iend2 = pialoc[ilist+1]-1; 
		jind1 = ibeg1;
		jind2 = ibeg2;

		while (jind1 <= iend1 || jind2 <= iend2) {
			if (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = ja[jind1];
				jj2 = pjaloc[jind2];
				if (jj1 < jj2) {
					jind1++;
				} else if (jj1 > jj2) {
					jind2++;
				} else if (jj1 == jj2) {
					nz++;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				jind1++;
			} else if (jind2 <= iend2) {
				jind2++;
			};
		};

		nelem[ilist] = nz;

	};

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) nz += nelem[ilist];

// Allocate work memory for extended matrix

	CSMatrixR temp (nlist,nz);

	temp.m = m;
	temp.n = n;

// Init list and ia arrays

	temp.ia[0] = 0;

	for (ilist=0;ilist<nlist;ilist++) {
		temp.list[ilist] = list[ilist];
		temp.ia[ilist+1] = temp.ia[ilist]+nelem[ilist];
	};

// Fill ja and a arrays

	nz = 0;

	for (ilist=0;ilist<nlist;ilist++) {

// Scan both rows at the same time

		ibeg1 = ia[ilist]; 
		iend1 = ia[ilist+1]-1; 
		ibeg2 = pialoc[ilist]; 
		iend2 = pialoc[ilist+1]-1; 
		jind1 = ibeg1;
		jind2 = ibeg2;

		while (jind1 <= iend1 || jind2 <= iend2) {
			if (jind1 <= iend1 && jind2 <= iend2) {
				jj1 = ja[jind1];
				jj2 = pjaloc[jind2];
				if (jj1 < jj2) {
					jind1++;
				} else if (jj1 > jj2) {
					jind2++;
				} else if (jj1 == jj2) {
					temp.ja[nz] = jj1;
					temp.a [nz] = a[jind1];
					nz++;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				jind1++;
			} else if (jind2 <= iend2) {
				jind2++;
			};
		};
	};

// Delete work arrays

	delete [] nelem;

// Return obtained matrix

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Cut the sparsity of the matrix according to the list of index numbers
//========================================================================================
CSMatrixR CSMatrixR::ExcludeList (int _nlist, int *_jalist) { // Cut the sparsity of the matrix according to the list of index numbers

//	const char *funcname = "ExcludeList";

// Allocate work memory for extended matrix

	CSMatrixR temp (nlist,nzja);

	temp.m = m;
	temp.n = n;
	temp.nsupc = n;
	temp.nsupr = n;

// Init new matrix

	temp.ia[0] = 0;

	int i, j, jj, jp;
	bool keep;

	int nz = 0;
	int ip = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			keep = true;
			if (ip < _nlist) {
				jp = _jalist[ip];
				while (ip < _nlist-1 && jp < j) {
					ip++;
					jp = _jalist[ip];
				};
				if (jp < j) {
					ip++;
				} else if (jp == j) {
					if (jj != i) keep = false;
				};
			};
			if (keep) {
				temp.ja[nz] = jj;
				temp.a[nz] = a[j];
				nz++;
			};
		};
		temp.ia[i+1] = nz;
		temp.list[i] = list[i];
	};

	temp.nzja   = nz;
	temp.nza    = nz;
	temp.nzatot = nz;

// Return obtained matrix

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Include into the sparsity of the new matrix only elements from the list of index numbers
//========================================================================================
CSMatrixR CSMatrixR::IncludeList (int _nlist, int *_jalist) { // Include into the sparsity of the new matrix only elements from the list of index numbers

	const char *funcname = "IncludeList";

// Count number of elements in each row of the list

	int *ianew, *iptr;

	ianew = new int [nlist+1];
	if (!ianew) MemoryFail (funcname);
	iptr = new int [nlist+1];
	if (!iptr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<=nlist;i++) ianew[i] = 0;

	int ip = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			if (ip < _nlist) {
				if (_jalist[ip] == j) {
					ianew[i+1]++;
					ip++;
				};
			};
		};
	};

	for (i=0;i<nlist;i++) ianew[i+1] = ianew[i]+ianew[i+1];

	for (i=0;i<nlist;i++) iptr[i] = ianew[i];

	int nzjanew = ianew[nlist];

	int *janew;
	double *anew;

	janew = new int [nzjanew];
	if (!janew) MemoryFail (funcname);
	anew = new double [nzjanew];
	if (!anew) MemoryFail (funcname);

	ip = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			if (ip < _nlist) {
				if (_jalist[ip] == j) {
					k = iptr[i];
					janew[k] = ja[j];
					anew[k] = a[j];
					iptr[i]++;
					ip++;
				};
			};
		};
	};

// Filter the list and reassign ianew

	int *listnew;
	int *ianew1;

	listnew = new int [nlist];
	if (!listnew) MemoryFail (funcname);
	ianew1 = new int [nlist+1];
	if (!ianew1) MemoryFail (funcname);

	int nlistnew = 0;

	ianew1[0] = 0;

	for (i=0;i<nlist;i++) {
		if (ianew[i+1] > ianew[i]) {
			listnew[nlistnew] = list[i];
			ianew1[nlistnew+1] = ianew[i+1];
			nlistnew++;
		};
	};

// Allocate work memory for extended matrix

	CSMatrixR temp (nlistnew,nzjanew);

// Copy elements

	for (i=0;i<nlistnew;i++) temp.list[i] = listnew[i];
	for (i=0;i<=nlistnew;i++) temp.ia[i] = ianew1[i];
	for (i=0;i<nzjanew;i++) temp.ja[i] = janew[i];
	for (i=0;i<nzjanew;i++) temp.a[i] = anew[i];

	temp.m = m;
	temp.n = n;
	temp.nsupc = n;
	temp.nsupr = n;
	temp.nlist  = nlistnew;
	temp.nzja   = nzjanew;
	temp.nza    = nzjanew;
	temp.nzatot = nzjanew;

// Free work arrays

	delete [] ianew;
	delete [] iptr;
	delete [] janew;
	delete [] anew;
	delete [] listnew;
	delete [] ianew1;

// Return obtained matrix

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Include into the sparsity of the new matrix only elements from the list of index numbers excluding diagonal ones
//========================================================================================
CSMatrixR CSMatrixR::IncludeListNoDiag (int _nlist, int *_jalist) { // Include into the sparsity of the new matrix only elements from the list of index numbers excluding diagonal ones

	const char *funcname = "IncludeListNoDiag";

// Count number of elements in each row of the list

	int *ianew, *iptr;

	ianew = new int [nlist+1];
	if (!ianew) MemoryFail (funcname);
	iptr = new int [nlist+1];
	if (!iptr) MemoryFail (funcname);

	int i, j, k;

	for (i=0;i<=nlist;i++) ianew[i] = 0;

	int ip = 0;

	int jj;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (ip < _nlist) {
				if (_jalist[ip] == j) {
					if (jj != i) ianew[i+1]++;
					ip++;
				};
			};
		};
	};

	for (i=0;i<nlist;i++) ianew[i+1] = ianew[i]+ianew[i+1];

	for (i=0;i<nlist;i++) iptr[i] = ianew[i];

	int nzjanew = ianew[nlist];

	int *janew;
	double *anew;

	janew = new int [nzjanew];
	if (!janew) MemoryFail (funcname);
	anew = new double [nzjanew];
	if (!anew) MemoryFail (funcname);

	ip = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			if (ip < _nlist) {
				if (_jalist[ip] == j) {
					if (jj != i) {
						k = iptr[i];
						janew[k] = ja[j];
						anew[k] = a[j];
						iptr[i]++;
					};
					ip++;
				};
			};
		};
	};

// Filter the list and reassign ianew

	int *listnew;
	int *ianew1;

	listnew = new int [nlist];
	if (!listnew) MemoryFail (funcname);
	ianew1 = new int [nlist+1];
	if (!ianew1) MemoryFail (funcname);

	int nlistnew = 0;

	ianew1[0] = 0;

	for (i=0;i<nlist;i++) {
//		if (ianew[i+1] > ianew[i]) {
			listnew[nlistnew] = list[i];
			ianew1[nlistnew+1] = ianew[i+1];
			nlistnew++;
//		};
	};

// Allocate work memory for extended matrix

	CSMatrixR temp (nlistnew,nzjanew);

// Copy elements

	for (i=0;i<nlistnew;i++) temp.list[i] = listnew[i];
	for (i=0;i<=nlistnew;i++) temp.ia[i] = ianew1[i];
	for (i=0;i<nzjanew;i++) temp.ja[i] = janew[i];
	for (i=0;i<nzjanew;i++) temp.a[i] = anew[i];

	temp.m = m;
	temp.n = n;
	temp.nsupc = n;
	temp.nsupr = n;
	temp.nlist  = nlistnew;
	temp.nzja   = nzjanew;
	temp.nza    = nzjanew;
	temp.nzatot = nzjanew;

// Free work arrays

	delete [] ianew;
	delete [] iptr;
	delete [] janew;
	delete [] anew;
	delete [] listnew;
	delete [] ianew1;

// Return obtained matrix

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute transposed matrix
//========================================================================================
void CSMatrixR::SplitAccordingBlockSparsity (int _nblks, int *_blks, // Split the matrix according to the specified block sparsity
															int *_iabmax, int *_jabmax, 
															CSMatrixR &_mtr2) {

	const char *funcname = "SplitAccordingBlockSparsity";

// Init the mask

	int *imaskja;

	imaskja = new int [nzja];
	if (!imaskja) MemoryFail (funcname);

	int i;

	for (i=0;i<nzja;i++) imaskja[i] = -1;

// Scan the matrix and mark first and second order elements

	int iblk = 0;
	int jblkprev = 0;

	int irow, jb, j, jcol, jblk, jbegblk, jendblk, jbbeg, k, jblksp;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		if (irow < _blks[iblk]) {
			throw " CSMatrixR::SplitAccordingBlockSparsity: nonmonotone list on entry ";
		} else if (irow >= _blks[iblk+1]) {
			while (irow >= _blks[iblk+1] && iblk < _nblks-1) iblk++;
		};
		jb = _iabmax[iblk];
		for (j=ia[i];j<ia[i+1];j++) {
			jcol = ja[j];
			if (jcol >= _blks[jblkprev] && jcol < _blks[jblkprev+1]) {
				jblk = jblkprev;
			} else {
				jbegblk = 0;
				jendblk = _nblks-1;
				while (jbegblk != jendblk) {
					jblk = (jbegblk+jendblk)/2;
					if (jcol >= _blks[jblk] && jcol < _blks[jblk+1]) {
						jbegblk = jblk;
						jendblk = jblk;
					} else if (jcol < _blks[jblk]) {
						jendblk = jblk-1;
					} else if (jcol >= _blks[jblk+1]) {
						jbegblk = jblk+1;
					};
				};
				jblk = jbegblk;
			};
			jblkprev = jblk;
			jbbeg = jb;
			for (k=jbbeg;k<_iabmax[iblk+1];k++) {
				jblksp = _jabmax[k];
				if (jblksp == jblk) {
					jb = k;
					imaskja[j] = 1;
					break;
				} else if (jblksp < jblk) {
					jb = k;
				} else if (jblksp > jblk) {
					break;
				};
			};
		};
	};

// Split the matrix according to the mask array

	int nzja2new = 0;

	for (i=0;i<nzja;i++) {
		if (imaskja[i] < 0) nzja2new++;
	};

	CSMatrixR mtr2loc (nlist,nzja2new);

	nzja2new = 0;
	mtr2loc.ia[0] = 0;

	for (i=0;i<nlist;i++) {
		mtr2loc.list[i] = list[i];
		for (j=ia[i];j<ia[i+1];j++) {
			if (imaskja[j] < 0) {
				mtr2loc.ja[nzja2new] = ja[j];
				mtr2loc.a[nzja2new] = a[j];
				nzja2new++;
			};
		};
		mtr2loc.ia[i+1] = nzja2new;
	};

	mtr2loc.m = m;
	mtr2loc.n = n;
	mtr2loc.nsupr = nsupr;
	mtr2loc.nsupc = nsupc;
	mtr2loc.nlist = nlist;
	mtr2loc.nzja = nzja2new;
	mtr2loc.nza = nzja2new;
	mtr2loc.nzatot = nzja2new;

	_mtr2 = mtr2loc;

	CSMatrixR mtrdummy;

	mtr2loc = mtrdummy;

// Filter current submatrix

	int *ialoc;

	ialoc = new int [nlist+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;
	int nzjaloc = 0;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			if (imaskja[j] >= 0) {
				ja[nzjaloc] = ja[j];
				a[nzjaloc] = a[j];
				nzjaloc++;
			};
		};
		ialoc[i+1] = nzjaloc;
	};

	for (i=0;i<=nlist;i++) ia[i] = ialoc[i];

	nzja = nzjaloc;
	nza = nzjaloc;
	nzatot = nzjaloc;

// Free work memory

	delete [] imaskja;
	delete [] ialoc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute transposed matrix
//========================================================================================
CSMatrixR CSMatrixR::TranspMtr () const { // Compute transposed matrix 

	const char *funcname = "TranspMtr";

	CSMatrixR at(n,nzja);

	int i, j, jj;

	for (i=0;i<=n;i++) at.ia[i] = 0;

	for (i=0;i<n;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			at.ia[jj+1]++;
		};
	};

	for (i=0;i<n;i++) at.ia[i+1] += at.ia[i];

	int *pat;

	pat = new int [n];

	if (!pat) MemoryFail (funcname);

	for (i=0;i<n;i++) pat[i] = at.ia[i];
	
	int k;
	for (i=0;i<n;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			pat[jj]++;
			k = pat[jj];
			at.ja[k-1] = i;
			at.a[k-1] = a[j];
		};
	};

	for (i=0;i<n;i++) at.list[i] = i;

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
// CSMatrixR: Compute symmetrized matrix
//========================================================================================
CSMatrixR CSMatrixR::SymmMtr () const { // Compute symmetrized matrix 

//	const char *funcname = "SymmMtr";

// Compute transposed matrix

	CSMatrixR at;

	at = TranspMtr();

// Compute the symmetric matrix

	int nzjas;
	nzjas = 2*nzja+n;

	CSMatrixR as(n,nzjas);

	int i, nz;
	int ibeg1, iend1, ibeg2, iend2, jind1, jind2, jj1, jj2;

	as.ia[0] = 0;
	nz = 0;

	for (i=0;i<n;i++) {
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
					as.a [nz] = a [jind1];
					nz++;
					jind1++;
				} else if (jj1 > jj2) {
					as.ja[nz] = at.ja[jind2];
					as.a [nz] = at.a [jind2];
					nz++;
					jind2++;
				} else if (jj1 == jj2) {
					as.ja[nz] = ja[jind1];
					as.a [nz] = a [jind1];
					nz++;
					jind1++;
					jind2++;
				};
			} else if (jind1 <= iend1) {
				as.ja[nz] = ja[jind1];
				as.a [nz] = a [jind1];
				nz++;
				jind1++;
			} else if (jind2 <= iend2) {
				as.ja[nz] = at.ja[jind2];
				as.a [nz] = at.a [jind2];
				nz++;
				jind2++;
			};
		};
		as.ia[i+1] = nz;
	};

	for (i=0;i<n;i++) as.list[i] = i;

	as.nzja = as.ia[n];
	as.nza = as.ia[n];
	as.nzatot = as.ia[n];

	return as;

};

// Author: Kharchenko S.A.
// CSMatrixR: Matrix ordering routine
//========================================================================================
CSMatrixR CSMatrixR::OrdMtr (const int *_order) const { // Matrix ordering routine

	const char *funcname = "OrdMtr";

// Compute inverse order

	int *iord;
	iord = new int [n];
	if (!iord) MemoryFail (funcname);

	int i;
	for (i=0;i<n;i++) iord[_order[i]] = i;

// Reorder the matrix

	int nzjas = nzja;

	CSMatrixR an(n,nzjas);

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
	double *pa;

	for (inew=0;inew<n;inew++) {
		i = iord[inew];
		nzloc = ia[i+1]-ia[i];
		ibs = nz;
		for (j=ia[i];j<ia[i+1];j++) {
			jold = ja[j];
			jnew = _order[jold];
			an.ja[nz] = jnew;
			an.a[nz] = a[j];
			nz++;
		};

// Sort elements

		pja = an.ja;
		pa = an.a;

		SortElm (nzloc, pja+ibs, pa+ibs);

	};

// Delete ordering array

	delete [] iord;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixR: Symmetric matrix ordering routine
//========================================================================================
CSMatrixR CSMatrixR::OrdMtrSymm (const int *_order) const { // Symmetric matrix ordering routine

	const char *funcname = "OrdMtrSymm";

// Compute inverse order

	int *iord;

	iord = new int [n];
	if (!iord) MemoryFail (funcname);

	int i;
	for (i=0;i<n;i++) iord[_order[i]] = i;

// Compute the symmetric matrix

	CSMatrixR at;

	at = TranspMtr ();

// Reorder the matrix

	CSMatrixR an(n,nzja);

	for (i=0;i<=n;i++) an.ia[i] = 0;

	int j, jj, jjnew, k;

	for (i=0;i<n;i++) {
		j = _order[i];
		for (k=ia[i];k<ia[i+1];k++) {
			jj = ja[k];
			jjnew = _order[jj];
			if (jjnew >= j) an.ia[j+1]++;
		};
		for (k=at.ia[i];k<at.ia[i+1];k++) {
			jj = at.ja[k];
			jjnew = _order[jj];
			if (jjnew > j) an.ia[j+1]++;
		};
	};

	for (i=0;i<n;i++) an.ia[i+1] += an.ia[i];

	int nz = 0;

	int inew;

	for (inew=0;inew<n;inew++) {
		i = iord[inew];
		for (k=ia[i];k<ia[i+1];k++) {
			jj = ja[k];
			jjnew = _order[jj];
			if (jjnew >= inew) {
				an.ja[nz] = jjnew;
				an.a[nz] = a[k];
				nz++;
			};
		};
		for (k=at.ia[i];k<at.ia[i+1];k++) {
			jj = at.ja[k];
			jjnew = _order[jj];
			if (jjnew > inew) {
				an.ja[nz] = jjnew;
				an.a[nz] = at.a[k];
				nz++;
			};
		};
	};

// Sort elements

	CIntDouble *idarr;

	idarr = new CIntDouble [n];
	if (!idarr) MemoryFail (funcname);

	int ibeg, jloc, nzloc;

	for (i=0;i<n;i++) {

		ibeg = an.ia[i];
		nzloc = an.ia[i+1]-an.ia[i];

		for (j=an.ia[i];j<an.ia[i+1];j++) {
			jloc = j-ibeg;
			idarr[jloc].intvalue = an.ja[j];
			idarr[jloc].dvalue = an.a[j];
		};

		qsort (idarr,nzloc,sizeof(CIntDouble),CIntDouble::CompareIntDouble);

		for (j=an.ia[i];j<an.ia[i+1];j++) {
			jloc = j-ibeg;
			an.ja[j] = idarr[jloc].intvalue;
			an.a[j] = idarr[jloc].dvalue;
		};

	};

// Delete work arrays

	delete [] iord;
	delete [] idarr;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixR: Sort matrix data by ordering the list (not indices)
//========================================================================================
CSMatrixR CSMatrixR::SortDataViaList (const int *_orderr) const { // Sort matrix data by ordering the list (not indices)

//	const char *funcname = "SortDataViaList";

// Reorder the matrix

	int i, j;

	CSMatrixR an(nlist,nzja);

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
			an.a[ibs+j-ia[i]] = a[j];
		};
	};

	an.n = n;
	an.m = m;
	an.nsupr = nsupr;
	an.nsupc = nsupc;
	an.nlist = nlist;
	an.nzja = nzja;
	an.nza = nza;
	an.nzatot = nzatot;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixR: Sort matrix data by ordering the list (not indices)
//========================================================================================
CSMatrixR CSMatrixR::Sort2IndexDataViaList (const int *_orderr) const { // Sort matrix data by ordering the list (not indices)

//	const char *funcname = "Sort2IndexDataViaList";

// Reorder the matrix

	int i, j;

	CSMatrixR an(nlist,nlist,nzja,nzja);

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
			an.a[ibs+j-ia[i]] = a[j];
		};
	};

	an.n = n;
	an.m = m;
	an.nsupr = nsupr;
	an.nsupc = nsupc;
	an.nlist = nlist;
	an.nlist2 = nlist;
	an.nzja = nzja;
	an.nzja2 = nzja;
	an.nza = nza;
	an.nzatot = nzatot;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixR: Super sparsify symmetrically the matrices L and U by numerically zero values
//========================================================================================
void CSMatrixR::SuperSparsifySymm2Index (bool _cpy3index, CSMatrixR &_mtrl, CSMatrixR &_mtru) { // Super sparsify symmetrically the matrices L and U by numerically zero values

	const char *funcname = "SuperSparsifySymm2Index";

	double dzero = 0.0e0;

// Create the mask of zero values

	int nzjal = _mtrl.nzja;

	int *jamask;

	jamask = new int [nzjal];
	if (!jamask) MemoryFail (funcname);

	int i, j;

	for (i=0;i<nzjal;i++) jamask[i] = 1;

	for (i=0;i<nzjal;i++) {
		if (_mtrl.a[i] == dzero && _mtru.a[i] == dzero) {
			jamask[i] = 0;
		};
	};

	int nzjanew = 0;

	for (i=0;i<nzjal;i++) {
		if (jamask[i] != 0) nzjanew++;
	};

// Replace ja, ja2 and a arrays

	int *janew;
	int *ja2new;
	int *ja3new;
	double *anew;

	janew = new int [nzjanew];
	if (!janew) MemoryFail (funcname);
	ja2new = new int [nzjanew];
	if (!ja2new) MemoryFail (funcname);
	if (_cpy3index) {
		ja3new = new int [nzjanew];
		if (!ja3new) MemoryFail (funcname);
	};
	anew = new double [nzjanew];
	if (!anew) MemoryFail (funcname);

	nzjanew = 0;

	for (i=0;i<nzjal;i++) {
		if (jamask[i] != 0) {
			janew[nzjanew] = _mtrl.ja[i];
			ja2new[nzjanew] = _mtrl.ja2[i];
			if (_cpy3index) {
				ja3new[nzjanew] = _mtrl.ja3[i];
			};
			anew[nzjanew] = _mtrl.a[i];
			nzjanew++;
		};
	};

	delete [] _mtrl.ja;
	delete [] _mtrl.ja2;
	if (_cpy3index) {
		delete [] _mtrl.ja3;
	};
	delete [] _mtrl.a;

	_mtrl.ja = janew;
	_mtrl.ja2 = ja2new;
	if (_cpy3index) {
		_mtrl.ja3 = ja3new;
	};
	_mtrl.a = anew;

	janew = new int [nzjanew];
	if (!janew) MemoryFail (funcname);
	ja2new = new int [nzjanew];
	if (!ja2new) MemoryFail (funcname);
	if (_cpy3index) {
		ja3new = new int [nzjanew];
		if (!ja3new) MemoryFail (funcname);
	};
	anew = new double [nzjanew];
	if (!anew) MemoryFail (funcname);

	nzjanew = 0;

	for (i=0;i<nzjal;i++) {
		if (jamask[i] != 0) {
			janew[nzjanew] = _mtru.ja[i];
			ja2new[nzjanew] = _mtru.ja2[i];
			if (_cpy3index) {
				ja3new[nzjanew] = _mtru.ja3[i];
			};
			anew[nzjanew] = _mtru.a[i];
			nzjanew++;
		};
	};

	delete [] _mtru.ja;
	delete [] _mtru.ja2;
	if (_cpy3index) {
		delete [] _mtru.ja3;
	};
	delete [] _mtru.a;

	_mtru.ja = janew;
	_mtru.ja2 = ja2new;
	if (_cpy3index) {
		_mtru.ja3 = ja3new;
	};
	_mtru.a = anew;

// Replace list, list2 and ia

	int nlistl = _mtrl.nlist;

	int *listmask;

	listmask = new int [nlistl];
	if (!listmask) MemoryFail (funcname);

	for (i=0;i<nlistl;i++) listmask[i] = 0;

	for (i=0;i<nlistl;i++) {
		for (j=_mtrl.ia[i];j<_mtrl.ia[i+1];j++) {
			if (jamask[j] == 1) listmask[i] = 1;
		};
	};

	int nlistnew = 0;

	for (i=0;i<nlistl;i++) {
		if (listmask[i] != 0) nlistnew++;
	};

	int *listloc, *list2loc, *ialoc;

	listloc = new int [nlistnew];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlistnew];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlistnew+1];
	if (!ialoc) MemoryFail (funcname);

	nlistnew = 0;

	for (i=0;i<nlistl;i++) {
		if (listmask[i] != 0) {
			listloc[nlistnew] = _mtrl.list[i];
			list2loc[nlistnew] = _mtrl.list2[i];
			nlistnew++;
		};
	};

	for (i=0;i<=nlistnew;i++) ialoc[i] = 0;

	nlistnew = 0;

	for (i=0;i<nlistl;i++) {
		if (listmask[i] != 0) {
			for (j=_mtrl.ia[i];j<_mtrl.ia[i+1];j++) {
				if (jamask[j] == 1) ialoc[nlistnew+1]++;
			};
			nlistnew++;
		};
	};

	for (i=0;i<nlistnew;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	delete [] _mtrl.list;
	delete [] _mtrl.list2;
	delete [] _mtrl.ia;

	_mtrl.list = new int [nlistnew];
	if (!_mtrl.list) MemoryFail (funcname);
	_mtrl.list2 = new int [nlistnew];
	if (!_mtrl.list2) MemoryFail (funcname);
	_mtrl.ia = new int [nlistnew+1];
	if (!_mtrl.ia) MemoryFail (funcname);

	for (i=0;i<nlistnew;i++) _mtrl.list[i] = listloc[i];
	for (i=0;i<nlistnew;i++) _mtrl.list2[i] = list2loc[i];
	for (i=0;i<=nlistnew;i++) _mtrl.ia[i] = ialoc[i];

	_mtrl.nlist = nlistnew;
	_mtrl.nlist2 = nlistnew;
	_mtrl.nzja = nzjanew;
	_mtrl.nzja2 = nzjanew;
	if (_cpy3index) {
		_mtrl.nzja3 = nzjanew;
	};
	_mtrl.nza = nzjanew;
	_mtrl.nzatot = nzjanew;

	delete [] _mtru.list;
	delete [] _mtru.list2;
	delete [] _mtru.ia;

	_mtru.list = listloc;
	_mtru.list2 = list2loc;
	_mtru.ia = ialoc;

	_mtru.nlist = nlistnew;
	_mtru.nlist2 = nlistnew;
	_mtru.nzja = nzjanew;
	_mtru.nzja2 = nzjanew;
	if (_cpy3index) {
		_mtru.nzja3 = nzjanew;
	};
	_mtru.nza = nzjanew;
	_mtru.nzatot = nzjanew;

// Free work arrays

	delete [] jamask;
	delete [] listmask;

};

// Author: Kharchenko S.A.
// CSMatrixR: Super sparsify the matrix L by numerically zero values
//========================================================================================
void CSMatrixR::SuperSparsifyL2Index (bool _cpy3index, CSMatrixR &_mtrl) { // Super sparsify the matrix L by numerically zero values

	const char *funcname = "SuperSparsifyL2Index";

	double dzero = 0.0e0;

// Create the mask of zero values

	int nzjal = _mtrl.nzja;

	int *jamask;

	jamask = new int [nzjal];
	if (!jamask) MemoryFail (funcname);

	int i, j;

	for (i=0;i<nzjal;i++) jamask[i] = 1;

	for (i=0;i<nzjal;i++) {
		if (_mtrl.a[i] == dzero) {
			jamask[i] = 0;
		};
	};

	int nzjanew = 0;

	for (i=0;i<nzjal;i++) {
		if (jamask[i] != 0) nzjanew++;
	};

// Replace ja, ja2 and a arrays

	int *janew;
	int *ja2new;
	int *ja3new;
	double *anew;

	janew = new int [nzjanew];
	if (!janew) MemoryFail (funcname);
	ja2new = new int [nzjanew];
	if (!ja2new) MemoryFail (funcname);
	if (_cpy3index) {
		ja3new = new int [nzjanew];
		if (!ja3new) MemoryFail (funcname);
	};
	anew = new double [nzjanew];
	if (!anew) MemoryFail (funcname);

	nzjanew = 0;

	for (i=0;i<nzjal;i++) {
		if (jamask[i] != 0) {
			janew[nzjanew] = _mtrl.ja[i];
			ja2new[nzjanew] = _mtrl.ja2[i];
			if (_cpy3index) {
				ja3new[nzjanew] = _mtrl.ja3[i];
			};
			anew[nzjanew] = _mtrl.a[i];
			nzjanew++;
		};
	};

	delete [] _mtrl.ja;
	delete [] _mtrl.ja2;
	if (_cpy3index) {
		delete [] _mtrl.ja3;
	};
	delete [] _mtrl.a;

	_mtrl.ja = janew;
	_mtrl.ja2 = ja2new;
	if (_cpy3index) {
		_mtrl.ja3 = ja3new;
	};
	_mtrl.a = anew;

// Replace list, list2 and ia

	int nlistl = _mtrl.nlist;

	int *listmask;

	listmask = new int [nlistl];
	if (!listmask) MemoryFail (funcname);

	for (i=0;i<nlistl;i++) listmask[i] = 0;

	for (i=0;i<nlistl;i++) {
		for (j=_mtrl.ia[i];j<_mtrl.ia[i+1];j++) {
			if (jamask[j] == 1) listmask[i] = 1;
		};
	};

	int nlistnew = 0;

	for (i=0;i<nlistl;i++) {
		if (listmask[i] != 0) nlistnew++;
	};

	int *listloc, *list2loc, *ialoc;

	listloc = new int [nlistnew];
	if (!listloc) MemoryFail (funcname);
	list2loc = new int [nlistnew];
	if (!list2loc) MemoryFail (funcname);
	ialoc = new int [nlistnew+1];
	if (!ialoc) MemoryFail (funcname);

	nlistnew = 0;

	for (i=0;i<nlistl;i++) {
		if (listmask[i] != 0) {
			listloc[nlistnew] = _mtrl.list[i];
			list2loc[nlistnew] = _mtrl.list2[i];
			nlistnew++;
		};
	};

	for (i=0;i<=nlistnew;i++) ialoc[i] = 0;

	nlistnew = 0;

	for (i=0;i<nlistl;i++) {
		if (listmask[i] != 0) {
			for (j=_mtrl.ia[i];j<_mtrl.ia[i+1];j++) {
				if (jamask[j] == 1) ialoc[nlistnew+1]++;
			};
			nlistnew++;
		};
	};

	for (i=0;i<nlistnew;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	delete [] _mtrl.list;
	delete [] _mtrl.list2;
	delete [] _mtrl.ia;

	_mtrl.list = new int [nlistnew];
	if (!_mtrl.list) MemoryFail (funcname);
	_mtrl.list2 = new int [nlistnew];
	if (!_mtrl.list2) MemoryFail (funcname);
	_mtrl.ia = new int [nlistnew+1];
	if (!_mtrl.ia) MemoryFail (funcname);

	for (i=0;i<nlistnew;i++) _mtrl.list[i] = listloc[i];
	for (i=0;i<nlistnew;i++) _mtrl.list2[i] = list2loc[i];
	for (i=0;i<=nlistnew;i++) _mtrl.ia[i] = ialoc[i];

	_mtrl.nlist = nlistnew;
	_mtrl.nlist2 = nlistnew;
	_mtrl.nzja = nzjanew;
	_mtrl.nzja2 = nzjanew;
	if (_cpy3index) {
		_mtrl.nzja3 = nzjanew;
	};
	_mtrl.nza = nzjanew;
	_mtrl.nzatot = nzjanew;

// Free work arrays

	delete [] jamask;
	delete [] listmask;
	delete [] listloc;
	delete [] list2loc;
	delete [] ialoc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute the list of small diagonal elements
//========================================================================================
void CSMatrixR::ListOfSmallDiagonalValues (double _thresh, int &_nlist, int *_list) const { // Compute the list of small diagonal elements

//	const char *funcname = "ListOfSmallDiagonalValues";

	int i, irow, j, jcol;
	double aux;

	_nlist = 0;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jcol = ja[j];
			if (irow == jcol) {
				aux = a[j];
				if (aux < 0.0e0) aux = -aux;
				if (aux <= _thresh) {
					_list[_nlist] = irow;
					_nlist++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Read matrix from disk in binary form
//========================================================================================
CSMatrixR CSMatrixR::ReadMtrBin (istream &_stream) { // Read matrix from disk in binary form

	int mloc, nloc, nlistl, nzjal;

	ReadArrAccV (_stream, 1, &mloc);
	ReadArrAccV (_stream, 1, &nloc);
	ReadArrAccV (_stream, 1, &nlistl);
	ReadArrAccV (_stream, 1, &nzjal);

	cout << " CSMatrixR:" << endl;
	cout << " M = " << mloc << " N = " << nloc << endl;
	cout << " Nlist = " << nlistl << " Nzja = " << nzjal << endl;

	CSMatrixR temp (nlistl,nzjal);

	temp.m      = mloc;
	temp.n      = nloc;
	temp.nsupr  = mloc;
	temp.nsupc  = nloc;
	temp.nlist  = nlistl;
	temp.nzja   = nzjal;
	temp.nza    = nzjal;
	temp.nzatot = nzjal;

	ReadArrAccV (_stream, nlistl, temp.list);
	ReadArrAccV (_stream, nlistl+1, temp.ia);
	ReadArrAccV (_stream, nzjal,    temp.ja);
	ReadArrAccV (_stream, nzjal,     temp.a);

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Read matrix from disk
//========================================================================================
CSMatrixR CSMatrixR::ReadMtr (istream &_fin) { // Read matrix from disk

	int nloc, nzjalc;
	char buffer[20]="";

	_fin >> nloc;
	cout << " N = " << nloc << endl;
	_fin >> nzjalc;
	cout << " Nzja = " << nzjalc << endl;

	CSMatrixR temp (nloc,nzjalc);

	_fin.getline(buffer,10);
	cout << buffer << endl;
	_fin.getline(buffer,10);
	cout << buffer << endl;

	int i;
	for (i=0; i<=nloc; i++)
	{
		_fin >> temp.ia[i];
	}

	_fin.getline(buffer,10);
	cout << buffer << endl;
	_fin.getline(buffer,10);
	cout << buffer << endl;

	for (i=0; i<nzjalc; i++)
	{
		_fin >> temp.ja[i];
	}

	_fin.getline(buffer,10);
	cout << buffer << endl;
	_fin.getline(buffer,10);
	cout << buffer << endl;

	for (i=0; i<nzjalc; i++)
	{
		_fin >> temp.a[i];
	}

	for (i=0;i<nloc;i++) temp.list[i] = i;

	temp.m      = nloc;
	temp.n      = nloc;
	temp.nsupr  = nloc;
	temp.nsupc  = nloc;
	temp.nlist  = nloc;
	temp.nzja   = nzjalc;
	temp.nza    = nzjalc;
	temp.nzatot = nzjalc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CSMatrixR &_a) { // Output matrix

	_stream << " CSMatrixR:" << endl;

	_stream << (const CSMatrix &) _a;

	_stream << " Nza = " << _a.nza << endl;
	_stream << " NzaTot = " << _a.nzatot << endl;

	OutArr (_stream, " A  =", _a.nza, _a.a);

	return _stream;

};

// Author: Kharchenko S.A.
// CSMatrixR: Output matrix to disk
//========================================================================================
void CSMatrixR::WriteMtr (ostream &_stream) { // Output matrix to disk

	_stream << m << endl;
	_stream << ia[nlist] << endl;
	_stream << "ia" << endl;
	int i;
	for (i=0;i<=nlist;i++) _stream << ia[i] << endl;
	_stream << "ja" << endl;
	for (i=0;i<=ia[nlist];i++) _stream << ja[i] << endl;
	_stream << "matrix" << endl;
	for (i=0;i<=ia[nlist];i++) _stream << a[i] << endl;

};

// Author: Kharchenko S.A.
// CSMatrixR: Output matrix to disk in binary form 
//========================================================================================
void CSMatrixR::WriteMtrBin (ostream &_stream) { // Output matrix to disk in binary form

	OutArrAccV (_stream, " M =", 1,          &m);
	OutArrAccV (_stream, " N =", 1,          &n);
	OutArrAccV (_stream, " Nlist =", 1,  &nlist);
	OutArrAccV (_stream, " Nzja =", 1, ia+nlist);
	OutArrAccV (_stream, " List =", nlist, list);
	OutArrAccV (_stream, " Ia =", nlist+1,   ia);
	OutArrAccV (_stream, " Ja =", ia[nlist], ja);
	OutArrAccV (_stream, " A =", ia[nlist],   a);

};

// Author: Kharchenko S.A.
// CSMatrixR: Read matrix, right hand side and norm vector from disk and write in binary form
//========================================================================================
void CSMatrixR::ConvertToBinary (ifstream &_fin, ofstream &_fout) { // Read matrix, right hand side and norm vector from disk and write in binary form

	if (!_fin.is_open()) throw " Error when opening file with data";

	CSMatrixR a;

	a = CSMatrixR::ReadMtr (_fin);

	int n = a.GetN ();

	CSVector b, norm, vectdummy;

	b = vectdummy.ReadVect (_fin,n);
	norm = vectdummy.ReadVect (_fin,n);

// Write matrix to the disk if necessary

	a.WriteMtrBin (_fout);

	b.WriteVectBin (_fout);
	norm.WriteVectBin (_fout);

	_fout.close ();

};

// Author: Kharchenko S.A.
// CSMatrixR: Read matrix, right hand side and norm vector from disk
//========================================================================================
void CSMatrixR::ReadFromBinary (ifstream &_fin, CSMatrix &_amatr, CSVector &_rhs, CSVector &_norm) { // Read matrix, right hand side and norm vector from disk

	if (!_fin.is_open()) throw " Error when opening file with data";

	_amatr = CSMatrixR::ReadMtrBin (_fin);

	_rhs = CSVector::ReadVectBin (_fin);
	_norm = CSVector::ReadVectBin (_fin);

};

// Author: Kharchenko S.A.
// CSMatrixR: Pack matrix data
//========================================================================================
void CSMatrixR::PackMatrix (int &_length, char *&_obj) { // Pack matrix data

	const char* funcname = "PackMatrix";

// Count the total size of the object

	_length = 11 * sizeof(int) + (2*nlist+1+nlist2+nzja+nzja2+nzja3) * sizeof(int) + (nza) * sizeof(double);

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
	double *pa;

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

	pa = (double *) pLoc;
	pLoc += nza * sizeof(double);

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

// Fill arrays

	int j;

	for(j = 0; j < nlist;   j++) plist[j] = list[j];
	for(j = 0; j < nlist2;  j++) plist2[j] = list2[j];
	for(j = 0; j < nlist+1; j++) pia[j] = ia[j];
	for(j = 0; j < nzja;    j++) pja[j] = ja[j];
	for(j = 0; j < nzja2;   j++) pja2[j] = ja2[j];
	for(j = 0; j < nzja3;   j++) pja3[j] = ja3[j];
	for(j = 0; j < nza;     j++) pa[j] = a[j];

};

// Author: Kharchenko S.A.
// CSMatrixR: UnPack matrix data
//========================================================================================
void CSMatrixR::UnPackMatrix (int _length, char *_obj) { // UnPack matrix data

	const char* funcname = "UnPackMatrix";

// Free previous arrays

	delete [] list;
	delete [] list2;
	delete [] ia;
	delete [] ja;
	delete [] ja2;
	delete [] ja3;
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
	double *pa;

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
	nzatot  = pHead[9];
	nza     = pHead[10];

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

	pa = (double *) pLoc;
	pLoc += nza * sizeof(double);

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
	a = new double [nza];
	if (!a) MemoryFail (funcname);

	int j;

	for(j = 0; j < nlist;   j++) list[j] = plist[j];
	for(j = 0; j < nlist2;  j++) list2[j] = plist2[j];
	for(j = 0; j < nlist+1; j++) ia[j] = pia[j];
	for(j = 0; j < nzja;    j++) ja[j] = pja[j];
	for(j = 0; j < nzja2;   j++) ja2[j] = pja2[j];
	for(j = 0; j < nzja3;   j++) ja3[j] = pja3[j];
	for(j = 0; j < nza;     j++) a[j] = pa[j];

};
