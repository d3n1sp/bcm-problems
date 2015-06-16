//------------------------------------------------------------------------------------------------
// File: smatrixc.cpp
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
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixC: Memory allocation zero data constructor
//========================================================================================
CSMatrixC::CSMatrixC (): CSMatrix () { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixC_00";

	dcmplx czero (0.0e0,0.0e0);

	nza = 0;
	nzatot = 0;
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = czero;
	};

};

// Author: Kharchenko S.A.
// CSMatrixC: Memory allocation zero data constructor
//========================================================================================
CSMatrixC::CSMatrixC (int _nlist, int _nzja): CSMatrix (_nlist,_nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixC_01";

	dcmplx czero (0.0e0,0.0e0);

	nza = _nzja;
	nzatot = nzja;
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = czero;
	};
};

// Author: Kharchenko S.A.
// CSMatrixC: Memory allocation zero data constructor
//========================================================================================
CSMatrixC::CSMatrixC (int _nlist, int _nzja, int _nza): CSMatrix (_nlist,_nzja) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixC_02";

	dcmplx czero (0.0e0,0.0e0);

	nza = _nza;
	nzatot = _nza;
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = czero;
	};
};

// Author: Kharchenko S.A.
// CSMatrixC: Memory allocation zero data constructor
//========================================================================================
CSMatrixC::CSMatrixC (int _nlist, int _nlist2, int _nzja, int _nzja2, int _nza): CSMatrix (_nlist,_nlist2,_nzja,_nzja2) { // Memory allocation zero data constructor

	const char *funcname = "CSMatrixR_03";

	dcmplx czero (0.0e0,0.0e0);

	nza = _nza;
	nzatot = _nza;
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = czero;
	};
};

// Author: Kharchenko S.A.
// CSMatrixC: Copy constructor
//========================================================================================
CSMatrixC::CSMatrixC (const CSMatrixC &_aa): CSMatrix (_aa) { // Copy constructor

	const char *funcname = "CSMatrixC_copy";

	nza = _aa.nza;
	nzatot = _aa.nzatot;
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = _aa.a[i];
	};

};

// Author: Kharchenko S.A.
// CSMatrixC: Equality operator
//========================================================================================
CSMatrixC &CSMatrixC::operator= (const CSMatrixC &_aa) { // Equality operator

	const char *funcname = "CSMatrixC_=";

	delete [] a;

	this->CSMatrix::operator= ((const CSMatrix &) _aa);

	nza = _aa.nza;
	nzatot = _aa.nzatot;
	a = new dcmplx [nza];
	if (!a) MemoryFail (funcname);
	for (int i=0;i<nza;i++) {
		a[i] = _aa.a[i];
	};

	return *this;

};

// Author: Kharchenko S.A.
// CSMatrixC: Ordering of the matrix columns
//========================================================================================
CSMatrixC CSMatrixC::OrdMtrRectCol (const int *_order) const { // Ordering of the matrix columns

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

	CSMatrixC an(nsupmx,nzja,nza);

// Copy list

	for (i=0;i<nsupc;i++) an.list[i] = list[i];

	int j;

// Create new ia

	for (i=0;i<=nsupc;i++) an.ia[i] = 0;

	for (i=0;i<nsupc;i++) {
		j = _order[i];
		an.ia[j+1] = ia[i+1]-ia[i];
	};

	for (i=0;i<nsupc;i++) an.ia[i+1] += an.ia[i];

// Create new ja and a

	int nz = 0;
	int isup;

	int isupnew;
	for (isupnew=0;isupnew<nsupc;isupnew++) {
		isup = iord[isupnew];
		for (j=ia[isup];j<ia[isup+1];j++) {
			an.ja[nz] = ja[j];
			an.a[nz] = a[j];
			nz++;
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

// Delete ordering array

	delete [] iord;

	return an;

};

// Author: Kharchenko S.A.
// CSMatrixC: Compute transposed matrix
//========================================================================================
CSMatrixC CSMatrixC::TranspMtr () const { // Compute transposed matrix 

	const char *funcname = "TranspMtr";

	CSMatrixC at(nsupr,nzja);

	int ilist;
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

	int k, ibs, kbs;

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			at.ja[k] = i;
			pat[jj]++;
		};
	};

//	int nz = 0;

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			ibs = j;
			kbs = k;
			at.a[kbs] = conj(a[ibs]);
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
	at.nza = nzja;
	at.nzatot = nzja;

	for (ilist=0;ilist<nsupr;ilist++) at.list[ilist] = ilist;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrixC: Compute transposed matrix without conjugation
//========================================================================================
CSMatrixC CSMatrixC::TranspMtrNoConj () const { // Compute transposed matrix without conjugation

	const char *funcname = "TranspMtrNoConj";

	CSMatrixC at(nsupr,nzja);

	int ilist;
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

	int k, ibs, kbs;

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			at.ja[k] = i;
			pat[jj]++;
		};
	};

//	int nz = 0;

	for (i=0;i<nsupr;i++) pat[i] = at.ia[i];

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = pat[jj];
			ibs = j;
			kbs = k;
			at.a[kbs] = a[ibs];
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
	at.nza = nzja;
	at.nzatot = nzja;

	for (ilist=0;ilist<nsupr;ilist++) at.list[ilist] = ilist;

	return at;

};

// Author: Kharchenko S.A.
// CSMatrixC: Perform matrix sparsification according to the column theshold (Tikhonov regularization block is included)
//========================================================================================
CSMatrixC CSMatrixC::SparsifyRectangularMatrix (double _alpha, double _thresh_col, // Perform matrix sparsification according to the column theshold (Tikhonov regularization block is included)
																int _m, int _n, dcmplx *_amatr) {

	const char *funcname = "SparsifyRectangularMatrix";

// Find the squared value of the largest element in each column

	double *valuearr;

	valuearr = new double [_n];
	if (!valuearr) MemoryFail (funcname);

	int i, j, ind;
	double aux, auxloc;

	for (i=0;i<_n;i++) {
		aux = 0.0e0;
		for (j=0;j<_m;j++) {
			ind = i*_m+j;
			auxloc = _amatr[ind].x*_amatr[ind].x+_amatr[ind].y*_amatr[ind].y;
			if (auxloc >= aux) aux = auxloc;
		};
		valuearr[i] = aux;
	};

// Count the number of elements in each column

	double thresh2 = _thresh_col * _thresh_col;

	int *nzcol;

	nzcol = new int [_n];
	if (!nzcol) MemoryFail (funcname);

	for (i=0;i<_n;i++) nzcol[i] = 1;

	int nzjatot = 0;

	for (i=0;i<_n;i++) {
		aux = thresh2 * valuearr[i];
		for (j=0;j<_m;j++) {
			ind = i*_m+j;
			auxloc = _amatr[ind].x*_amatr[ind].x+_amatr[ind].y*_amatr[ind].y;
			if (auxloc >= aux) nzcol[i]++;
		};
		nzjatot += nzcol[i];
	};

	int *ialoc;

	ialoc = new int [_n+1];
	if (!ialoc) MemoryFail (funcname);

	ialoc[0] = 0;

	for (i=0;i<_n;i++) ialoc[i+1] = ialoc[i]+nzcol[i];

// Create new matrix

	CSMatrixC temp (_n, nzjatot);

	int *pia = temp.GetIa ();
	int *plist = temp.GetList ();
	int *pja = temp.GetJa ();
	dcmplx *pa = temp.GetA ();

	for (i=0;i<_n;i++) plist[i] = i;

	int nz = 0;

	pia[0] = 0;

	for (i=0;i<_n;i++) {
		aux = thresh2 * valuearr[i];
		for (j=0;j<_m;j++) {
			ind = i*_m+j;
			auxloc = _amatr[ind].x*_amatr[ind].x+_amatr[ind].y*_amatr[ind].y;
			if (auxloc >= aux && nz < ialoc[i+1]-1) {
				pja[nz] = j;
				pa[nz] = _amatr[ind];
				nz++;
			};
		};
		pja[nz] = _m+i;
		pa[nz].x = _alpha;
		pa[nz].y = 0.0e0;
		nz++;
		pia[i+1] = nz;
	};

	temp.SetM     (_m+_n);
	temp.SetN     (_n);
	temp.SetNsupr (_m+_n);
	temp.SetNsupcBase (_n);
	temp.SetNlist (_n);
	temp.SetNzja  (nz);
	temp.SetNza   (nz);
	temp.SetNzatot(nz);

// Free work arrays

	delete [] valuearr;
	delete [] nzcol;
	delete [] ialoc;

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixC: Pack matrix data
//========================================================================================
void CSMatrixC::PackMatrix (int &_length, char *&_obj) { // Pack matrix data

	const char* funcname = "PackMatrix";

// Count the total size of the object

	_length = 12 * sizeof(int) + (2*nlist+1+nlist2+nzja+nzja2+nzja3) * sizeof(int) + (nza) * sizeof(dcmplx);

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
	pHead[11] = 0;

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
// CSMatrixC: UnPack matrix data
//========================================================================================
void CSMatrixC::UnPackMatrix (int _length, char *_obj) { // UnPack matrix data

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
	a = new dcmplx [nza];
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

// Author: Kharchenko S.A.
// CSMatrixC: Output matrix
//========================================================================================
ostream &operator<< (ostream &_stream, const CSMatrixC &_a) { // Output matrix

	_stream << " CSMatrixC:" << endl;

	_stream << (const CSMatrix &) _a;

	_stream << " Nza = " << _a.nza << endl;
	_stream << " NzaTot = " << _a.nzatot << endl;

//	double thresh = 1.0e-13;

	OutArr (_stream, " A =", _a.nza, _a.a);

	return _stream;

};
