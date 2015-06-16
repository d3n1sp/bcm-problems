//------------------------------------------------------------------------------------------------
// File: mvmc.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
//#include <complex>
#include <cmath>

#include "tree.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "mvm.h"
#include "corr.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CMvmC: Constructor
//========================================================================================
CMvmC::CMvmC (bool _is2index, const CTree &_tree, const CGSMatrixCS &_gmtra, int _nrhs): CMvm (_is2index, _tree, _gmtra) { // Constructor

	const char *funcname = "CMvmC_00";

// Allocate the data

	CSVectorC temp(nlocext,1,_nrhs);

	px = temp;
	pax = temp;

	CSVectorC vectdummy;

	temp = vectdummy;

// Allocate and init work arrays

	int *imask;
	int *listblk;

	imask = new int [_gmtra.nblksr];
	if (!imask) MemoryFail (funcname);
	listblk = new int [_gmtra.nblksr];
	if (!listblk) MemoryFail (funcname);

	int icycle = -1;

	for (int i=0;i<_gmtra.nblksr;i++) imask[i] = icycle;

// Determine the tree type

	int treetype = 0;

	int ilev;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = lev2node[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
			if (inodetype != compnode) treetype = 1;
		};
	};

// Create the list of blocks to be added

	icycle++;

	int nlistloc = 0;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = lev2node[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
//			int nodecpu = _tree.nodes[inode].nodecpu;
			int iblkbeg = _tree.nodes[inode].indbeg;
			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
			int nchildsloc = _tree.nodes[inode].nchilds;

			if (inodetype != exchangenode) {

// Next node if this is the degenerate node

				if (nchildsloc <= 1) goto NextNode;

// Combine all previous data into the receive list

				if (treetype == 0) {
					for (int jblk=iblkbeg;jblk<=iblkend;jblk++) {
						if (imask[jblk] != icycle) {
							listblk[nlistloc] = jblk;
							nlistloc++;
							imask[jblk] = icycle;
						};
					};
				} else {
					for (int jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
						if (imask[jblk] != icycle) {
							listblk[nlistloc] = jblk;
							nlistloc++;
							imask[jblk] = icycle;
						};
					};
				};
			};
		};

NextNode:;

	};

	qsort (listblk, nlistloc, sizeof(int), compint);

// Compute nlocext and bspx

	int nadd = 0;

	for (int ilist=0;ilist<nlistloc;ilist++) {
		int iblk = listblk[ilist];
		for (int isup=_gmtra.blksr[iblk];isup<_gmtra.blksr[iblk+1];isup++) {
			int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
			nadd += blai;
		};
	};

	vectsndrcv = new dcmplx [_nrhs*nadd];
	if (!vectsndrcv) MemoryFail (funcname);

	nrhs = _nrhs;

	delete [] imask;
	delete [] listblk;

};

// Author: Kharchenko S.A.
// CMvmC: Destructor
//========================================================================================
CMvmC::~CMvmC () { // Destructor
//	std::cout << " On entry to CMvmC destructor " << std::endl;

	delete [] vectsndrcv;

//	std::cout << " On return from CMvmC destructor " << std::endl;
};

// Author: Kharchenko S.A.
// CSMatrixCS: General matrix by vector multiplication
//========================================================================================
void CSMatrixCS::MvmA (const CSVectorC &_x, CSVectorC &_ax) const { // General matrix by vector multiplication

//	const char *funcname = "MvmA_01";

	dcmplx czero (0.0e0,0.0e0);

	int i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv, irhs;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	_ax.SetSVect (czero);

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0;i<nsupr;i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = sprndc[jj+1]-sprndc[jj];
				jbsv = sprndc[jj];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
//						*pax += *pa * *px++;
						(*pax).PlEq (*pa * *px++);
						pa += blai;
					};
					pax++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: General matrix by vector multiplication in out-of-core mode
//========================================================================================
void CSMatrixCS::MvmA (int _nblks, int *_blks, FILE **_mtrafiles, // General matrix by vector multiplication in out-of-core mode
						const CSVectorC &_x, CSVectorC &_ax) const {

	const char *funcname = "MvmA_02";

	dcmplx czero (0.0e0,0.0e0);

	int iblk,i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv, irhs;

	dcmplx *px, *pax, *pa, *pablk;

	dcmplx *celem;

	celem = new dcmplx [blamx*blamx];
	if (!celem) MemoryFail (funcname);

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	_ax.SetSVect (czero);

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				blaj = sprndc[jj+1]-sprndc[jj];

				int blaij = blai*blaj;
				jbsv = sprndc[jj];
				ibs = bsa[j];

				FGet (_mtrafiles[iblk],blaij,celem,ibs);

				for (irhs=0;irhs<nrhsx;irhs++) {
					pablk = celem;
					pax = _ax.vect+irhs*nvax+ibsv;
					for (kii=0;kii<blai;kii++) {
						px = _x.vect+irhs*nvx+jbsv;
						pa = pablk+kii;
						for (kjj=0;kjj<blaj;kjj++) {
//							_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
//							*pax += *pa * *px++;
							(*pax).PlEq (*pa * *px++);
							pa += blai;
						};
						pax++;
					};
				};
			};
		};
	};

	delete [] celem;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Mupliply by the current block row of Au
//========================================================================================
void CSMatrixCS::MvmBlockRowAu (int *_sprnds, // Mupliply by the current block row of Au
								int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const {

//	const char *funcname = "MvmBlockRowAu";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsax[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprnds[jj+1]-_sprnds[jj];
				jbsv = _bsx[jj];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
//						*pax += *pa * *px++;
						(*pax).PlEq (*pa * *px++);
						pa += blai;
					};
					pax++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Mupliply by the complex conjugate current block row of Au
//========================================================================================
void CSMatrixCS::MvmBlockRowAuConj (int *_sprnds, // Mupliply by the complex conjugate current block row of Au
								int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const {

//	const char *funcname = "MvmBlockRowAuConj";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsax[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprnds[jj+1]-_sprnds[jj];
				jbsv = _bsx[jj];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
//						*pax += *pa * *px++;
						(*pax).PlEq (conj(*pa) * *px++);
						pa += blai;
					};
					pax++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Mupliply by the current block column of Al without diagonal supernode
//========================================================================================
void CSMatrixCS::MvmBlockColumnAl (int *_sprnds, // Mupliply by the current block column of Al without diagonal supernode
									int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const {

//	const char *funcname = "MvmBlockColumnAl";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs, kiib;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			for (j=ia[ilist]+1;j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprnds[jj+1]-_sprnds[jj];
				jbsv = _bsax[jj];
				pax = _ax.vect+irhs*nvax+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_ax.vect[jbsv+kii] += a[ibs+kii*blai+kjj] * _x.vect[ibsv+kjj];
//						*pax += *pa++ * *px++;
						(*pax).PlEq (*pa++ * *px++);
					};
					pax++;
					kiib += blai;
				};
			};
		};
	};

};

//========================================================================
// CSMatrixCS: Mupliply by the complex conjugate current block column of Al without diagonal supernode
//========================================================================================
void CSMatrixCS::MvmBlockColumnAlConj (int *_sprnds, // Mupliply by the complex conjugate current block column of Al without diagonal supernode
										int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const {

//	const char *funcname = "MvmBlockColumnAlConj";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs, kiib;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			for (j=ia[ilist]+1;j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprnds[jj+1]-_sprnds[jj];
				jbsv = _bsax[jj];
				pax = _ax.vect+irhs*nvax+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_ax.vect[jbsv+kii] += a[ibs+kii*blai+kjj] * _x.vect[ibsv+kjj];
//						*pax += *pa++ * *px++;
						(*pax).PlEq (conj(*pa++) * *px++);
					};
					pax++;
					kiib += blai;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Mupliply by the current block column of A
//========================================================================================
void CSMatrixCS::MvmBlockColumnA (int *_sprndr, int *_sprndc, // Mupliply by the current block column of A
									int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const {

//	const char *funcname = "MvmBlockColumnA";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs, kiib;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprndc[i+1]-_sprndc[i];
			ibsv = _bsx[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprndr[jj+1]-_sprndr[jj];
				jbsv = _bsax[jj];
				pax = _ax.vect+irhs*nvax+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_ax.vect[jbsv+kii] += a[ibs+kii*blai+kjj] * _x.vect[ibsv+kjj];
//						*pax += *pa++ * *px++;
						(*pax).PlEq (*pa++ * *px++);
					};
					pax++;
					kiib += blai;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Mupliply by the complex conjugate transposed current block column of rectangular A
//========================================================================================
void CSMatrixCS::MvmBlockColumnAConj (int *_sprndr, int *_sprndc, // Mupliply by the complex conjugate current block column of rectangular A
													int *_bsx, const CSVectorC &_x, int *_bsax, CSVectorC &_ax) const {

//	const char *funcname = "MvmBlockColumnAConj";

// Attention: ax on entry should be initialized by zeroes if necessary

	int isup, j, jsup, kii, kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	dcmplx *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			isup = list[ilist];
			blai = _sprndc[isup+1]-_sprndc[isup];
			ibsv = _bsax[isup];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jsup = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprndr[jsup+1]-_sprndr[jsup];
				jbsv = _bsx[jsup];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
//						*pax += *pa * *px++;
						(*pax).PlEq (conj(*pa) * *px++);
						pa += blai;
					};
					pax++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General matrix by vector multiplication
//========================================================================================
void CGSMatrixCS::MvmA (const CGSMatrixCS &_gmtral, const CGSMatrixCS &_gmtrau, // General matrix by vector multiplication
						const CSVectorC &_x, CSVectorC &_ax) const {

//	const char *funcname = "MvmA_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize ax by zeroes

	int nrhsx = _x.nrhs;
	int nvax = _ax.nv;

	for (int kii=0;kii<nrhsx*nvax;kii++) _ax.vect[kii] = czero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<_gmtral.nblksr;iblk++) {

		_gmtrau.mtrarr[iblk].MvmBlockRowAu    (_gmtral.sprndr, _gmtral.sprndr, _x, _gmtral.sprndr, _ax);

		_gmtral.mtrarr[iblk].MvmBlockColumnAl (_gmtral.sprndr, _gmtral.sprndr, _x, _gmtral.sprndr, _ax);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Complex hermitian matrix by vector multiplication
//========================================================================================
void CGSMatrixCS::MvmASymm ( // Complex hermitian matrix by vector multiplication
						const CSVectorC &_x, CSVectorC &_ax) {

//	const char *funcname = "MvmA_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize ax by zeroes

	int nrhsx = _x.nrhs;
	int nvax = _ax.nv;

	for (int kii=0;kii<nrhsx*nvax;kii++) _ax.vect[kii] = czero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		ReadBlock (iblk);

		mtrarr[iblk].MvmBlockRowAu    (sprndr, sprndr, _x, sprndr, _ax);

		mtrarr[iblk].MvmBlockColumnAlConj (sprndr, sprndr, _x, sprndr, _ax);

		FreeBlock (iblk);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Complex hermitian matrix by vector multiplication
//========================================================================================
void CGSMatrixCS::MvmAhAOrdShiftSymm (double _dshift, int *_order,  // Complex hermitian matrix by vector multiplication
										const CSVectorC &_x, CSVectorC &_ax) {

//	const char *funcname = "MvmAhAOrdShiftSymm";

	dcmplx czero (0.0e0,0.0e0);

// Get sizes

	int mloc = this->m;
	int nloc = this->n;

// Create local work data

	CSVectorC xtemp (nloc);
	CSVectorC axtemp (mloc);

// Reorder initial data

	int i;

	for (i=0;i<nloc;i++) {
//		xtemp.vect[_order[i]] = _x.vect[i];
		xtemp.vect[i] = _x.vect[_order[i]];
	};

// Multiply twice

	this->MvmARect  (xtemp, axtemp);
	this->MvmAhRect (axtemp, xtemp);

// Reorder back and shift

	for (i=0;i<nloc;i++) {
//		_ax.vect[i] = xtemp.vect[_order[i]];
		_ax.vect[_order[i]] = xtemp.vect[i];
	};

	for (i=0;i<nloc;i++) {
		_ax.vect[i] += _x.vect[i]*_dshift;
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Complex hermitian matrix by vector multiplication
//========================================================================================
void CSVectorC::MvmAhAOrdShiftSymm (void *_obj, CGETM _getm, CGETN _getn, CMVMC _mvmc, CMVMHC _mvmhc, // Complex hermitian matrix by vector multiplication
												double _dshift, int *_order,
												const CSVectorC &_x, CSVectorC &_ax) const {

//	const char *funcname = "MvmAhAOrdShiftSymm";

	dcmplx czero (0.0e0,0.0e0);

// Get sizes

	int mloc = (_getm) (_obj);
	int nloc = (_getm) (_obj);

// Create local work data

	CSVectorC xtemp (nloc);
	CSVectorC axtemp (mloc);

// Reorder initial data

	int i;

	for (i=0;i<nloc;i++) {
//		xtemp.vect[_order[i]] = _x.vect[i];
		xtemp.vect[i] = _x.vect[_order[i]];
	};

// Multiply twice

	(_mvmc)  (_obj, xtemp, axtemp);
	(_mvmhc)  (_obj, axtemp, xtemp);

// Reorder back and shift

	for (i=0;i<nloc;i++) {
//		_ax.vect[i] = xtemp.vect[_order[i]];
		_ax.vect[_order[i]] = xtemp.vect[i];
	};

	for (i=0;i<nloc;i++) {
		_ax.vect[i] += _x.vect[i]*_dshift;
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixCS::MvmASymm (const CTree &_tree, CMvmC &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixCS &_gmtrau, 
							const CSVectorC &_x, CSVectorC &_ax) const {

//	const char *funcname = "MvmA_01";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.pax.vect[kii] = czero;

// Initialize px by x data

	int irhs;
	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_mvm.px.vect[irhs*_mvm.nlocext+kii] = _x.vect[irhs*_mvm.nloc+kii];
		};
	};

// Search the tree

	int inode, inodetype;
	int fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;
	int ilev;

	for (ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode1;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid, _tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'R', '=', _tree.myid, childcpuid, _tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

				};

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'r', _tree.myid, childcpuid, _tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid2, _tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

			};

		};

NextNode1:;

	};

// Perform all local computations

	for (int ilist=0;ilist<_gmtrau.nlist;ilist++) {
		int iblk = _gmtrau.listb[ilist];

		_gmtrau.ReadBlock (iblk);
		_gmtrau.mtrarr[iblk].MvmBlockRowAu    (_gmtrau.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtrau.mtrarr[iblk].MvmBlockColumnAlConj (_gmtrau.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtrau.FreeBlock (iblk);

	};

// Search the tree one more time

	nodebeg = _tree.cpuidbeg[_tree.myid];
	nodeend = _tree.cpuidend[_tree.myid];

	ilevbeg = _tree.nodes[nodebeg].nodelv;
	ilevend = _tree.nodes[nodeend].nodelv;

	for (ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode2;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'A', '+', _tree.myid, childcpuid, 2*_tree.nnodes+inode,
										_mvm.pax,
										nlistloc, listloc,
										_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid2, 2*_tree.nnodes+inode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid, 2*_tree.nnodes+fathernode,
								_mvm.pax,
								nlistloc, listloc,
								_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'R', _tree.myid, childcpuid, 2*_tree.nnodes+childnode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtrau.blksr, _gmtrau.sprndr, _gmtrau.bl2ndr);

				};

			};

		};

NextNode2:;

	};

// Store the result of computations in ax

	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ax.vect[irhs*_mvm.nloc+kii] = _mvm.pax.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixCS::MvmA (const CTree &_tree, CMvmC &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
							const CSVectorC &_x, CSVectorC &_ax) const {

//	const char *funcname = "MvmA_01";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.pax.vect[kii] = czero;

// Initialize px by x data

	int irhs;
	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_mvm.px.vect[irhs*_mvm.nlocext+kii] = _x.vect[irhs*_mvm.nloc+kii];
		};
	};

// Search the tree

	int inode, inodetype;
	int fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;
	int ilev;

	for (ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode1;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid, _tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'R', '=', _tree.myid, childcpuid, _tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'r', _tree.myid, childcpuid, _tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid2, _tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

NextNode1:;

	};

// Perform all local computations

	for (int ilist=0;ilist<_gmtral.nlist;ilist++) {
		int iblk = _gmtral.listb[ilist];

		_gmtrau.ReadBlock (iblk);
		_gmtrau.mtrarr[iblk].MvmBlockRowAu    (_gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtrau.FreeBlock (iblk);

		_gmtral.ReadBlock (iblk);
		_gmtral.mtrarr[iblk].MvmBlockColumnAl (_gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtral.FreeBlock (iblk);

	};

// Search the tree one more time

	nodebeg = _tree.cpuidbeg[_tree.myid];
	nodeend = _tree.cpuidend[_tree.myid];

	ilevbeg = _tree.nodes[nodebeg].nodelv;
	ilevend = _tree.nodes[nodeend].nodelv;

	for (ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode2;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'A', '+', _tree.myid, childcpuid, 2*_tree.nnodes+inode,
										_mvm.pax,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid2, 2*_tree.nnodes+inode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid, 2*_tree.nnodes+fathernode,
								_mvm.pax,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'R', _tree.myid, childcpuid, 2*_tree.nnodes+childnode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

NextNode2:;

	};

// Store the result of computations in ax

	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ax.vect[irhs*_mvm.nloc+kii] = _mvm.pax.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: General hermitian transposed matrix by vector multiplication 
//========================================================================================
void CSMatrixCS::MvmAh (const CSVectorC &_x, CSVectorC &_ax) const { // General hermitian transposed matrix by vector multiplication

//	const char *funcname = "MvmAh_00";

	dcmplx czero (0.0e0,0.0e0);

	int irhs,i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	dcmplx *px, *pax, *pa, *pablk;

	_ax.SetSVect (czero);

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0;i<nsupr;i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = sprndc[jj+1]-sprndc[jj];
				jbsv = sprndc[jj];
				px = _x.vect+irhs*nvx+ibsv;
				for (kii=0;kii<blai;kii++) {
					pax = _ax.vect+irhs*nvax+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[jbsv+kjj] += a[ibs+kjj*blai+kii] * _x.vect[ibsv+kii];
//						*pax++ += conj(*pa) * *px;
						(*pax++).PlEq (*px * conj(*pa));
						pa += blai;
					};
					px++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General rectangular matrix by vector multiplication
//========================================================================================
void CGSMatrixCS::MvmARect ( // General rectangular matrix by vector multiplication
									const CSVectorC &_x, CSVectorC &_ax) {

//	const char *funcname = "MvmARect_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize ax by zeroes

	int nrhsx = _x.nrhs;
	int nvax = _ax.nv;

	for (int kii=0;kii<nrhsx*nvax;kii++) _ax.vect[kii] = czero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksc;iblk++) {

		ReadBlock (iblk);

		mtrarr[iblk].MvmBlockColumnA (sprndr, sprndc, sprndc, _x, sprndr, _ax);

		FreeBlock (iblk);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General complex conjugate transposed rectangular matrix by vector multiplication
//========================================================================================
void CGSMatrixCS::MvmAhRect ( // General complex conjugate transposed rectangular matrix by vector multiplication
									const CSVectorC &_x, CSVectorC &_ax) {

//	const char *funcname = "MvmAhRect_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize ax by zeroes

	int nrhsx = _x.nrhs;
	int nvax = _ax.nv;

	for (int kii=0;kii<nrhsx*nvax;kii++) _ax.vect[kii] = czero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksc;iblk++) {

		ReadBlock (iblk);

		mtrarr[iblk].MvmBlockColumnAConj (sprndr, sprndc, sprndr, _x, sprndc, _ax);

		FreeBlock (iblk);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General complex conjugate transposed rectangular matrix by vector multiplication
//========================================================================================
void CGSMatrixCS::MvmAhRect (const CTree &_tree, // General complex conjugate transposed rectangular matrix by vector multiplication
									const CSVectorC &_x, CSVectorC &_ax) {

//	const char *funcname = "MvmAhRect_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize ax by zeroes

	int nrhsx = _x.nrhs;
	int nvax = _ax.nv;

	for (int kii=0;kii<nrhsx*nvax;kii++) _ax.vect[kii] = czero;

// Cycle over the block columns 

	int nnodestree = _tree.GetNnodes ();
	CNode *pnode = _tree.GetNode ();

//	int nproc = _tree.GetNproc ();
	int myid = _tree.GetMyid ();

	int iblk, inode, iproc, indbeg, indend, j;

	for (inode=0;inode<nnodestree;inode++) {
		iproc = pnode[inode].GetNodecpu ();
		indbeg = pnode[inode].GetIndbeg ();
		indend = pnode[inode].GetIndend ();
		if (iproc == myid) {
			for (j=indbeg;j<=indend;j++) {
				iblk = j;

				ReadBlock (iblk);

				mtrarr[iblk].MvmBlockColumnAConj (sprndr, sprndc, sprndr, _x, sprndc, _ax);

				FreeBlock (iblk);
			};
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: General matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixCS::MvmAh (const CTree &_tree, CMvmC &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixCS &_gmtral, CGSMatrixCS &_gmtrau, 
							const CSVectorC &_x, CSVectorC &_ax) const {

//	const char *funcname = "MvmAh_01";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.pax.vect[kii] = czero;

// Initialize px by x data

	int irhs;
	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_mvm.px.vect[irhs*_mvm.nlocext+kii] = _x.vect[irhs*_mvm.nloc+kii];
		};
	};

// Search the tree

	int inode, inodetype;
	int fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;
	int ilev;

	for (ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode1;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid, _tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'R', '=', _tree.myid, childcpuid, _tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'r', _tree.myid, childcpuid, _tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid2, _tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

NextNode1:;

	};

// Perform all local computations

	for (int ilist=0;ilist<_gmtral.nlist;ilist++) {
		int iblk = _gmtral.listb[ilist];

		_gmtral.ReadBlock (iblk);
		_gmtral.mtrarr[iblk].MvmBlockRowAuConj    (_gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtral.FreeBlock (iblk);

		_gmtrau.ReadBlock (iblk);
		_gmtrau.mtrarr[iblk].MvmBlockColumnAlConj (_gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtrau.FreeBlock (iblk);

	};

// Search the tree one more time

	nodebeg = _tree.cpuidbeg[_tree.myid];
	nodeend = _tree.cpuidend[_tree.myid];

	ilevbeg = _tree.nodes[nodebeg].nodelv;
	ilevend = _tree.nodes[nodeend].nodelv;

	for (ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode2;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'A', '+', _tree.myid, childcpuid, 2*_tree.nnodes+inode,
										_mvm.pax,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid2, 2*_tree.nnodes+inode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid, 2*_tree.nnodes+fathernode,
								_mvm.pax,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'R', _tree.myid, childcpuid, 2*_tree.nnodes+childnode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

NextNode2:;

	};

// Store the result of computations in ax

	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ax.vect[irhs*_mvm.nloc+kii] = _mvm.pax.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with L stored by columns
//========================================================================================
void CSMatrixCS::SolveL (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by columns

	const char *funcname = "SolveL_02";

	int i,kii,kjj,j;
	int nloc=n;

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	int blai, blaj, ibs, ibsv, jbsv, kiib, irhs;

	int nrhsx = _x.nrhs;
//	int nvx = _x.nv;
	int nvlx = _lx.nv;

	dcmplx *px, *plx, *pl, *pablk;

	px = _x.vect;
	plx = _lx.vect;

	for (i=0;i<nloc*nrhsx;i++) *plx++ = *px++;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0; i<nsupr; i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			int i1 = ia[i]+1;
			int i2 = ia[i+1];
			j = ia[i];
			ibs = bsa[j];
			pablk = a+ibs;
			plx = vloc;
			for (kii=0;kii<blai;kii++) *plx++ = czero;
			plx = vloc;
			kiib = 0;
			for (kii=0;kii<blai;kii++) {
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk+kiib;
				for (kjj=0;kjj<blai;kjj++) {
//					vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
//					*plx += *pl++ * *px++;
					(*plx).PlEq (*pl++ * *px++);
				};
				plx++;
				kiib += blai;
			};
			px = vloc;
			plx = _lx.vect+irhs*nvlx+ibsv;
			for (kii=0;kii<blai;kii++) *plx++ = *px++;
			for (int j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = sprndr[k+1]-sprndr[k];
				jbsv = sprndr[k];
				ibs = bsa[j];
				pablk = a+ibs;
				plx = _lx.vect+irhs*nvlx+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
//						*plx -= *pl++ * *px++;
						(*plx).MnEq (*pl++ * *px++);
					};
					plx++;
					kiib += blai;
				};
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with L stored by columns in out-of-core mode
//========================================================================================
void CSMatrixCS::SolveL (int _nblks, int *_blks, FILE **_mtrlfiles, // Solve triangular system with L stored by columns in out-of-core mode
							const CSVectorC &_x, CSVectorC &_lx) const {

	const char *funcname = "SolveL_03";

	dcmplx czero (0.0e0,0.0e0);

	int i,kii,kjj,j;
	int nloc=n;

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Allocate work arrays for reading 

	int nsupmx = 0;

	for (i=0;i<nsupc;i++) {
		int nsupi = ia[i+1]-ia[i];
		if (nsupi>nsupmx) nsupmx = nsupi;
	};

	dcmplx *aloc;

	aloc = new dcmplx [nsupmx*blamx*blamx];
	if (!aloc) MemoryFail (funcname);

// Direct solve

	int blai, blaj, ibs, ibs0, ibsv, jbsv, kiib;
	int iblk, irhs;

	dcmplx *px, *plx, *pl, *pablk;

	int nrhsx = _x.nrhs;
//	int nvx = _x.nv;
	int nvlx = _lx.nv;

	px = _x.vect;
	plx = _lx.vect;

	for (i=0;i<nloc*nrhsx;i++) *plx++ = *px++;

// Main cycle over the block columns

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {

			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];

// Read data of the supernode

			ibs0 = bsa[ia[i]];

			int nzloc = 0;

			for (j=ia[i];j<ia[i+1];j++) {
				int jj = ja[j];
				blaj = sprndr[jj+1]-sprndr[jj];
				nzloc += blai*blaj;
			};

			FGet (_mtrlfiles[iblk],nzloc,aloc,ibs0);

// Solve

			for (irhs=0;irhs<nrhsx;irhs++) {
				int i1 = ia[i]+1;
				int i2 = ia[i+1];
				j = ia[i];
				ibs = bsa[j];
				pablk = aloc+ibs-ibs0;
				plx = vloc;
				for (kii=0;kii<blai;kii++) *plx++ = czero;
				plx = vloc;
				kiib = 0;
				for (kii=0;kii<blai;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
//						*plx += *pl++ * *px++;
						(*plx).PlEq (*pl++ * *px++);
					};
					plx++;
					kiib += blai;
				};
				px = vloc;
				plx = _lx.vect+irhs*nvlx+ibsv;
				for (kii=0;kii<blai;kii++) *plx++ = *px++;
				for (int j=i1; j<i2; j++) {
					int k = ja[j];
					blaj = sprndr[k+1]-sprndr[k];
					jbsv = sprndr[k];
					ibs = bsa[j];
					pablk = aloc+ibs-ibs0;
					plx = _lx.vect+irhs*nvlx+jbsv;
					kiib = 0;
					for (kii=0;kii<blaj;kii++) {
						px = _lx.vect+irhs*nvlx+ibsv;
						pl = pablk+kiib;
						for (kjj=0;kjj<blai;kjj++) {
//							_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
//							*plx -= *pl++ * *px++;
							(*plx).MnEq (*pl++ * *px++);
						};
						plx++;
						kiib += blai;
					};
				};
			};
		};
	};

	delete [] vloc;
	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with block column of L
//========================================================================================
void CSMatrixCS::SolveBlockColumnL (int *_sprnds, // Solve triangular system with block column of L
									int *_bsx, const CSVectorC &_x, 
									int *_bslx, CSVectorC &_lx) const { 

	const char *funcname = "SolveBlockColumnL";

// Attention: lx on entry should be initialized (by zeroes or update values, not by x data) !!!

	int i,kii,kjj,j;

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	int ilist, blai, blaj, ibs, ibsv, jbsv, kiib, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvlx = _lx.nv;

	dcmplx *px, *plx, *pl, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			jbsv = _bslx[i];
			px = _x.vect+irhs*nvx+ibsv;
			plx = _lx.vect+irhs*nvlx+jbsv;
			for (kii=0;kii<blai;kii++) *plx++ += *px++;
		};
	};

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bslx[i];
			int i1 = ia[ilist]+1;
			int i2 = ia[ilist+1];
			j = ia[ilist];
			ibs = bsa[j];
			pablk = a+ibs;
			plx = vloc;
			for (kii=0;kii<blai;kii++) *plx++ = czero;
			plx = vloc;
			kiib = 0;
			for (kii=0;kii<blai;kii++) {
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk+kiib;
				for (kjj=0;kjj<blai;kjj++) {
//					vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
//					*plx += *pl++ * *px++;
					(*plx).PlEq (*pl++ * *px++);
				};
				plx++;
				kiib += blai;
			};
			ibsv = _bslx[i];
			px = vloc;
			plx = _lx.vect+irhs*nvlx+ibsv;
			for (kii=0;kii<blai;kii++) *plx++ = *px++;
			for (int j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = _sprnds[k+1]-_sprnds[k];
				jbsv = _bslx[k];
				ibs = bsa[j];
				pablk = a+ibs;
				plx = _lx.vect+irhs*nvlx+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
//						*plx -= *pl++ * *px++;
						(*plx).MnEq (*pl++ * *px++);
					};
					plx++;
					kiib += blai;
				};
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with the complex conjugate block column of L
//========================================================================================
void CSMatrixCS::SolveBlockColumnLConj (int *_sprnds, // Solve triangular system with the complex conjugate block column of L
										int *_bsx, const CSVectorC &_x, 
										int *_bslx, CSVectorC &_lx) const { 

	const char *funcname = "SolveBlockColumnLConj";

// Attention: lx on entry should be initialized (by zeroes or update values, not by x data) !!!

	int i,kii,kjj,j;

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	int ilist, blai, blaj, ibs, ibsv, jbsv, kiib, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvlx = _lx.nv;

	dcmplx *px, *plx, *pl, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			jbsv = _bslx[i];
			px = _x.vect+irhs*nvx+ibsv;
			plx = _lx.vect+irhs*nvlx+jbsv;
			for (kii=0;kii<blai;kii++) *plx++ += *px++;
		};
	};

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bslx[i];
			int i1 = ia[ilist]+1;
			int i2 = ia[ilist+1];
			j = ia[ilist];
			ibs = bsa[j];
			pablk = a+ibs;
			plx = vloc;
			for (kii=0;kii<blai;kii++) *plx++ = czero;
			plx = vloc;
			kiib = 0;
			for (kii=0;kii<blai;kii++) {
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk+kiib;
				for (kjj=0;kjj<blai;kjj++) {
//					vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
//					*plx += *pl++ * *px++;
					(*plx).PlEq (conj(*pl) * *px++);
					pl++;
				};
				plx++;
				kiib += blai;
			};
			ibsv = _bslx[i];
			px = vloc;
			plx = _lx.vect+irhs*nvlx+ibsv;
			for (kii=0;kii<blai;kii++) *plx++ = *px++;
			for (int j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = _sprnds[k+1]-_sprnds[k];
				jbsv = _bslx[k];
				ibs = bsa[j];
				pablk = a+ibs;
				plx = _lx.vect+irhs*nvlx+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
//						*plx -= *pl++ * *px++;
						(*plx).MnEq (conj(*pl) * *px++);
						pl++;
					};
					plx++;
					kiib += blai;
				};
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with L stored by block columns
//========================================================================================
void CGSMatrixCS::SolveL (const CSVectorC &_x, CSVectorC &_lx) const { // Solve triangular system with L stored by block columns

//	const char *funcname = "SolveL_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize lx by zeroes

	int nrhsx = _x.nrhs;
	int nvlx = _lx.nv;

	for (int kii=0;kii<nrhsx*nvlx;kii++) _lx.vect[kii] = czero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		mtrarr[iblk].SolveBlockColumnL (sprndr, sprndr, _x, sprndr, _lx);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with L stored by block columns
//========================================================================================
void CGSMatrixCS::SolveLConj (const CSVectorC &_x, CSVectorC &_lx) { // Solve triangular system with L stored by block columns

//	const char *funcname = "SolveLConj_00";

	dcmplx czero (0.0e0,0.0e0);

// Initialize lx by zeroes

	int nrhsx = _x.nrhs;
	int nvlx = _lx.nv;

	for (int kii=0;kii<nrhsx*nvlx;kii++) _lx.vect[kii] = czero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		ReadBlock(iblk);

		mtrarr[iblk].SolveBlockColumnLConj (sprndr, sprndr, _x, sprndr, _lx);

		FreeBlock(iblk);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with L stored by block columns in parallel mode
//========================================================================================
void CGSMatrixCS::SolveL (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with L stored by block columns in parallel mode
							const CSVectorC &_x, CSVectorC &_lx) {

//	const char *funcname = "SolveL_01";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px by zeroes

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;

// Search the tree

	int inode, inodetype;
	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevend;ilev>=ilevbeg;ilev--) {

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

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary and store data

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'A', '+', _tree.myid, childcpuid, 3*_tree.nnodes+inode,
											_mvm.px,
											nlistloc, listloc,
											blksr, sprndr, bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid2, 3*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

				ReadBlock (iblk);
				mtrarr[iblk].SolveBlockColumnL (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid, 3*_tree.nnodes+fathernode,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'R', _tree.myid, childcpuid, 3*_tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

				};

			};

		};

NextNode:;

	};

// Store the result of computations in lx

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_lx.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with complex conjugate L stored by block columns in parallel mode
//========================================================================================
void CGSMatrixCS::SolveLConj (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with complex conjugate L stored by block columns in parallel mode
								const CSVectorC &_x, CSVectorC &_lx) {

//	const char *funcname = "SolveLConj";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px by zeroes

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;

// Search the tree

	int inode, inodetype;
	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevend;ilev>=ilevbeg;ilev--) {

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

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary and store data

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'A', '+', _tree.myid, childcpuid, 3*_tree.nnodes+inode,
											_mvm.px,
											nlistloc, listloc,
											blksr, sprndr, bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid2, 3*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

				ReadBlock (iblk);
				mtrarr[iblk].SolveBlockColumnLConj (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid, 3*_tree.nnodes+fathernode,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (_tree.comm, 'R', _tree.myid, childcpuid, 3*_tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

				};

			};

		};

NextNode:;

	};

// Store the result of computations in lx

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_lx.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with U stored by rows
//========================================================================================
void CSMatrixCS::SolveU (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by rows

	const char *funcname = "SolveU_02";

	int i, kii, kjj;
//	int nloc=n;

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	int blai, blaj, ibs, ibsv, jbsv, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	dcmplx *px, *pux, *pu, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=nsupr-1;i>=0;i--) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			int i1 = ia[i]+1;
			int i2 = ia[i+1];
			pux = vloc;
			px = _x.vect+irhs*nvx+ibsv;
			for (kii=0;kii<blai;kii++) {
//				vloc[kii] = _x.vect[ibsv+kii];
				*pux++ = *px++;
			};
			int j;
			for (j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = sprndr[k+1]-sprndr[k];
				jbsv = sprndr[k];
				ibs = bsa[j];
				pablk = a+ibs;
				pux = vloc;
				for (kii=0;kii<blai;kii++) {
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
//						*pux -= *pu * *px++;
						(*pux).MnEq (*pu * *px++);
						pu += blai;
					};
					pux++;
				};
			};
			j = ia[i];
			ibs = bsa[j];
			pablk = a+ibs;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) *pux++ = czero;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) {
				px = vloc;
				pu = pablk+kii;
				for (kjj=0;kjj<blai;kjj++) {
//					_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
//					*pux += *pu * *px++;
					(*pux).PlEq (*pu * *px++);
					pu += blai;
				};
				pux++;
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with U stored by rows in out-of-core mode
//========================================================================================
void CSMatrixCS::SolveU (int _nblks, int *_blks, FILE **_mtrufiles, // Solve triangular system with U stored by rows in out-of-core mode
							const CSVectorC &_x, CSVectorC &_ux) const {

	const char *funcname = "SolveU_03";

	dcmplx czero (0.0e0,0.0e0);

	int i, kii, kjj;
//	int nloc=n;

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Allocate work arrays for reading 

	int nsupmx = 0;

	for (i=0;i<nsupc;i++) {
		int nsupi = ia[i+1]-ia[i];
		if (nsupi>nsupmx) nsupmx = nsupi;
	};

	dcmplx *aloc;

	aloc = new dcmplx [nsupmx*blamx*blamx];
	if (!aloc) MemoryFail (funcname);

// Backward solve

	int blai, blaj, ibs, ibs0, ibsv, jbsv;
	int iblk, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	dcmplx *px, *pux, *pu, *pablk;

	for (iblk=_nblks-1;iblk>=0;iblk--) {
		for (i=_blks[iblk+1]-1;i>=_blks[iblk];i--) {

			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];

// Read data of the supernode

			ibs0 = bsa[ia[i]];

			int nzloc = 0;

			int j;
			for (j=ia[i];j<ia[i+1];j++) {
				int jj = ja[j];
				blaj = sprndr[jj+1]-sprndr[jj];
				nzloc += blai*blaj;
			};

			FGet (_mtrufiles[iblk],nzloc,aloc,ibs0);

// Solve

			for (irhs=0;irhs<nrhsx;irhs++) {
				int i1 = ia[i]+1;
				int i2 = ia[i+1];
				pux = vloc;
				px = _x.vect+irhs*nvx+ibsv;
				for (kii=0;kii<blai;kii++) {
//					vloc[kii] = _x.vect[ibsv+kii];
					*pux++ = *px++;
				};
				for (j=i1; j<i2; j++) {
					int k = ja[j];
					blaj = sprndr[k+1]-sprndr[k];
					jbsv = sprndr[k];
					ibs = bsa[j];
					pablk = aloc+ibs-ibs0;
					pux = vloc;
					for (kii=0;kii<blai;kii++) {
						px = _ux.vect+irhs*nvux+jbsv;
						pu = pablk+kii;
						for (kjj=0;kjj<blaj;kjj++) {
//							vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
//							*pux -= *pu * *px++;
							(*pux).MnEq (*pu * *px++);
							pu += blai;
						};
						pux++;
					};
				};
				j = ia[i];
				ibs = bsa[j];
				pablk = aloc+ibs-ibs0;
				pux = _ux.vect+irhs*nvux+ibsv;
				for (kii=0;kii<blai;kii++) *pux++ = czero;
				pux = _ux.vect+irhs*nvux+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = vloc;
					pu = pablk+kii;
					for (kjj=0;kjj<blai;kjj++) {
//						_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
//						*pux += *pu * *px++;
						(*pux).PlEq (*pu * *px++);
						pu += blai;
					};
					pux++;
				};
			};
		};
	};

	delete [] vloc;
	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with the block row of U
//========================================================================================
void CSMatrixCS::SolveBlockRowU (int *_sprnds, // Solve triangular system with the block row of U
									int *_bsx, const CSVectorC &_x, 
									int *_bsux, CSVectorC &_ux) const {

	const char *funcname = "SolveBlockRowU";

	int i, kii, kjj;
//	int nloc=n;

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	dcmplx *px, *pux, *pu, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=nlist-1;ilist>=0;ilist--) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			int i1 = ia[ilist]+1;
			int i2 = ia[ilist+1];
			pux = vloc;
			px = _x.vect+irhs*nvx+ibsv;
			for (kii=0;kii<blai;kii++) {
//				vloc[kii] = _x.vect[ibsv+kii];
				*pux++ = *px++;
			};
			int j;
			for (j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = _sprnds[k+1]-_sprnds[k];
				jbsv = _bsux[k];
				ibs = bsa[j];
				pablk = a+ibs;
				pux = vloc;
				for (kii=0;kii<blai;kii++) {
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
//						*pux -= *pu * *px++;
						(*pux).MnEq (*pu * *px++);
						pu += blai;
					};
					pux++;
				};
			};
			j = ia[ilist];
			ibs = bsa[j];
			pablk = a+ibs;
			ibsv = _bsux[i];
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) *pux++ = czero;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) {
				px = vloc;
				pu = pablk+kii;
				for (kjj=0;kjj<blai;kjj++) {
//					_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
//					*pux += *pu * *px++;
					(*pux).PlEq (*pu * *px++);
					pu += blai;
				};
				pux++;
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve triangular system with the complex conjugate block row of U
//========================================================================================
void CSMatrixCS::SolveBlockRowUConj (int *_sprnds, // Solve triangular system with the complex conjugate block row of U
									int *_bsx, const CSVectorC &_x, 
									int *_bsux, CSVectorC &_ux) const {

	const char *funcname = "SolveBlockRowUConj";

	int i, kii, kjj;
//	int nloc=n;

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *vloc;

	vloc = new dcmplx [blamx];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	dcmplx *px, *pux, *pu, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=nlist-1;ilist>=0;ilist--) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			int i1 = ia[ilist]+1;
			int i2 = ia[ilist+1];
			pux = vloc;
			px = _x.vect+irhs*nvx+ibsv;
			for (kii=0;kii<blai;kii++) {
//				vloc[kii] = _x.vect[ibsv+kii];
				*pux++ = *px++;
			};
			int j;
			for (j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = _sprnds[k+1]-_sprnds[k];
				jbsv = _bsux[k];
				ibs = bsa[j];
				pablk = a+ibs;
				pux = vloc;
				for (kii=0;kii<blai;kii++) {
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
//						*pux -= *pu * *px++;
						(*pux).MnEq (conj(*pu) * *px++);
						pu += blai;
					};
					pux++;
				};
			};
			j = ia[ilist];
			ibs = bsa[j];
			pablk = a+ibs;
			ibsv = _bsux[i];
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) *pux++ = czero;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) {
				px = vloc;
				pu = pablk+kii;
				for (kjj=0;kjj<blai;kjj++) {
//					_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
//					*pux += *pu * *px++;
					(*pux).PlEq (conj(*pu) * *px++);
					pu += blai;
				};
				pux++;
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with U stored by block rows
//========================================================================================
void CGSMatrixCS::SolveU (const CSVectorC &_x, CSVectorC &_ux) const { // Solve triangular system with U stored by block rows

//	const char *funcname = "SolveU_00";

	dcmplx czero (0.0e0,0.0e0);

// Cycle over the block rows

	int iblk;

	for (iblk=nblksr-1;iblk>=0;iblk--) {

		mtrarr[iblk].SolveBlockRowU (sprndr, sprndr, _x, sprndr, _ux);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with U stored by block rows in parallel mode
//========================================================================================
void CGSMatrixCS::SolveU (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with U stored by block rows in parallel mode
							const CSVectorC &_x, CSVectorC &_ux) {

//	const char *funcname = "SolveU_01";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;

// Search the tree

	int inode, inodetype;

	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevbeg;ilev<=ilevend;ilev++) {

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

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid, 4*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'R', '=', _tree.myid, childcpuid, 4*_tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

				};

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkend;iblk>=iblkbeg;iblk--) {

				ReadBlock (iblk);
				mtrarr[iblk].SolveBlockRowU (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			if (_tree.nodes[inode].nchilds > 1) {

				for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

					int childnode = _tree.nodes[inode].childs[ichild];
					int childcpuid = _tree.nodes[childnode].nodecpu;

					int nlistloc = _mvm.nlistdnsndarr[ilev];
					int *listloc = _mvm.listdnsndarr[ilev];

					if (childcpuid != _tree.myid) {

						_mvm.SendVectData (_tree.comm, 'r', _tree.myid, childcpuid, 4*_tree.nnodes+childnode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

					};

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid2, 4*_tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		};

NextNode:;

	};

// Store the result of computations in ux

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ux.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Solve triangular system with complex conjugate U stored by block rows in parallel mode
//========================================================================================
void CGSMatrixCS::SolveUConj (const CTree &_tree, CMvmC &_mvm, // Solve triangular system with complex conjugate U stored by block rows in parallel mode
							const CSVectorC &_x, CSVectorC &_ux) {

//	const char *funcname = "SolveUConj";

	dcmplx czero (0.0e0,0.0e0);

// Initialize px

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = czero;

// Search the tree

	int inode, inodetype;

	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevbeg;ilev<=ilevend;ilev++) {

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

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (_tree.comm, 'A', '=', _tree.myid, fathercpuid, 4*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (_tree.comm, 'R', '=', _tree.myid, childcpuid, 4*_tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

				};

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkend;iblk>=iblkbeg;iblk--) {

				ReadBlock (iblk);
				mtrarr[iblk].SolveBlockRowUConj (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			if (_tree.nodes[inode].nchilds > 1) {

				for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

					int childnode = _tree.nodes[inode].childs[ichild];
					int childcpuid = _tree.nodes[childnode].nodecpu;

					int nlistloc = _mvm.nlistdnsndarr[ilev];
					int *listloc = _mvm.listdnsndarr[ilev];

					if (childcpuid != _tree.myid) {

						_mvm.SendVectData (_tree.comm, 'r', _tree.myid, childcpuid, 4*_tree.nnodes+childnode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

					};

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (_tree.comm, 'A', _tree.myid, fathercpuid2, 4*_tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		};

NextNode:;

	};

// Store the result of computations in ux

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ux.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Solve system with ordered preconditioner
//========================================================================================
void CSMatrixCS::SolveLUOrd (int *_order, CSMatrixCS &_mtrl, CSMatrixCS &_mtru, // Solve system with ordered preconditioner
										CSVectorC &_x, CSVectorC &_px) {

//	const char *funcname = "SolveLUOrd";

	dcmplx czero (0.0e0,0.0e0);

// Get sizes

	int nloc = _mtrl.n;
	int nrhsloc = _x.GetNrhs ();

// Create local work data

	CSVectorC xtemp (nloc,nrhsloc);
	CSVectorC axtemp (nloc,nrhsloc);

// Reorder initial data

	int irhs, i;

	for (irhs=0;irhs<nrhsloc;irhs++) {
		for (i=0;i<nloc;i++) {
			xtemp.vect[irhs*nloc+_order[i]] = _x.vect[irhs*nloc+i];
//			xtemp.vect[i] = _x.vect[_order[i]];
		};
	};

// Solve

	_mtrl.SolveL  (xtemp, axtemp);
	_mtru.SolveU (axtemp, xtemp);

// Reorder back and shift

	for (irhs=0;irhs<nrhsloc;irhs++) {
		for (i=0;i<nloc;i++) {
			_px.vect[irhs*nloc+i] = xtemp.vect[irhs*nloc+_order[i]];
//			_px.vect[_order[i]] = xtemp.vect[i];
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform block Daxpy y += x*a
//========================================================================================
void CSVectorC::BlockDaxpy (const CSVectorC &_x, int _lda, dcmplx *_a) { // Perform block Daxpy y += x*a

//	const char *funcname = "BlockDaxpy";

	dcmplx czero (0.0e0,0.0e0);

// Perform updates

	int ncolx = _x.nrhs;
	int ncoly = nrhs;
	int nrows = nv;

	int kii, kjj, kkk;

	for (kii=0;kii<nrows;kii++) {
		for (kjj=0;kjj<ncoly;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<ncolx;kkk++) {
				caux += _x.vect[kkk*nrows+kii] * _a[kjj*_lda+kkk];
			};
			vect[kjj*nrows+kii] += caux;
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform block Daxpy y += x*a
//========================================================================================
void CSVectorC::BlockDaxpy (int _ncolx, int _ldx, dcmplx *_x, int _lda, dcmplx *_a) { // Perform block Daxpy y += x*a

//	const char *funcname = "BlockDaxpy";

	dcmplx czero (0.0e0,0.0e0);

// Perform updates

	int ncolx = _ncolx;
	int ncoly = nrhs;
	int nrows = nv;

	int kii, kjj, kkk;

	for (kii=0;kii<nrows;kii++) {
		for (kjj=0;kjj<ncoly;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<ncolx;kkk++) {
				caux += _x[kkk*_ldx+kii] * _a[kjj*_lda+kkk];
			};
			vect[kjj*nrows+kii] += caux;
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform block Dot a = y^h*x
//========================================================================================
void CSVectorC::BlockDot (const CSVectorC &_x, int _lda, dcmplx *_a) const { // Perform block Dot a = y^h*x

//	const char *funcname = "BlockDot";

	dcmplx czero (0.0e0,0.0e0);

// Perform updates

	int ncolx = _x.nrhs;
	int ncoly = nrhs;
	int nrows = nv;

	int kii, kjj, kkk;

	for (kii=0;kii<ncoly;kii++) {
		for (kjj=0;kjj<ncolx;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<nrows;kkk++) {
				caux += _x.vect[kjj*nrows+kkk] * conj(vect[kii*nrows+kkk]);
			};
			_a[kjj*_lda+kii] = caux;
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform block Dot a = y^h*x
//========================================================================================
void CSVectorC::BlockDot (int _ncolx, int _ldx, dcmplx *_x, int _lda, dcmplx *_a) const { // Perform block Dot a = y^h*x

//	const char *funcname = "BlockDot";

	dcmplx czero (0.0e0,0.0e0);

// Perform updates

	int ncolx = _ncolx;
	int ncoly = nrhs;
	int nrows = nv;

	int kii, kjj, kkk;

	for (kii=0;kii<ncoly;kii++) {
		for (kjj=0;kjj<ncolx;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<nrows;kkk++) {
				caux += _x[kjj*_ldx+kkk] * conj(vect[kii*nrows+kkk]);
			};
			_a[kjj*_lda+kii] = caux;
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform block Dot a = x^h*y
//========================================================================================
void CSVectorC::BlockDotH (int _ncolx, int _ldx, dcmplx *_x, int _lda, dcmplx *_a) const { // Perform block Dot a = x^h*y

//	const char *funcname = "BlockDotH";

	dcmplx czero (0.0e0,0.0e0);

// Perform updates

	int ncolx = _ncolx;
	int ncoly = nrhs;
	int nrows = nv;

	int kii, kjj, kkk;

	for (kii=0;kii<ncolx;kii++) {
		for (kjj=0;kjj<ncoly;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<nrows;kkk++) {
				caux += vect[kjj*nrows+kkk] * conj(_x[kii*_ldx+kkk]);
			};
			_a[kjj*_lda+kii] = caux;
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform multiplication by the corrector matrix
//========================================================================================
void CSVectorC::MvCorr (const CTree &_tree, // Perform multiplication by the corrector matrix
						const CSVectorC &_ycorr, const CSVectorC &_wcorr,
						int _lddcorr, dcmplx *_dcorr, 
						int _ldcoef, dcmplx *_coef1, dcmplx *_coef2) {

//	const char *funcname = "MvCorr";

	dcmplx czero (0.0e0,0.0e0);

	int ncorr = _ycorr.nrhs;

	if (ncorr == 0) goto ExitMvCorr;

// Dot operations

	int kii, kjj, kkk;

	for (kii=0;kii<_ldcoef*nrhs;kii++) _coef1[kii] = czero;

	_ycorr.BlockDot (*this, _ldcoef, _coef1);
	_wcorr.BlockDot (*this, _ldcoef, _coef1+ncorr);

// Perform multiprocessor dot operation if necessary

	if (_tree.nproc != 1) {
		CMPIComm pcomm = _tree.GetComm ();
		CMPIExchange::ExchangeArrayMPI (pcomm,
													DOUBLEVALUE, ADD,
													_ldcoef*nrhs*2, _coef1, _coef2);
	} else {
		for (kii=0;kii<_ldcoef*nrhs;kii++) _coef2[kii] = _coef1[kii];
	};

// Multiply by the block diagonal

	for (kii=0;kii<_ldcoef*nrhs;kii++) _coef1[kii] = czero;

	for (kii=0;kii<ncorr*2;kii++) {
		for (kjj=0;kjj<nrhs;kjj++) {
			dcmplx caux = czero;
			for (kkk=0;kkk<ncorr*2;kkk++) {
				caux += _dcorr[kkk*_lddcorr+kii] * _coef2[kjj*_ldcoef+kkk];
			};
			_coef1[kjj*_ldcoef+kii] = caux;
		};
	};

// Daxpy operations

	BlockDaxpy (_ycorr, _ldcoef, _coef1);
	BlockDaxpy (_wcorr, _ldcoef, _coef1+ncorr);

ExitMvCorr:;

};

// Author: Kharchenko S.A.
// CSVectorC: Perform multiplication by the corrector matrix
//========================================================================================
void CSVectorC::MvCorr (const CTree &_tree, // Perform multiplication by the corrector matrix
						const CCorrectorC &_corr,
						int _ldcoef, dcmplx *_coef1, dcmplx *_coef2) {

//	const char *funcname = "MvCorr";

	dcmplx czero (0.0e0,0.0e0);

	int ncorr = _corr.ncorr;

	int kii, kjj, kkk;
	int icorr;

	if (ncorr == 0) goto ExitMvCorr;

// Main corrector cycle 

	for (icorr=ncorr-1;icorr>=0;icorr--) {

		int isize = _corr.corrsizes[icorr];

// Dot operations

		for (kii=0;kii<_ldcoef*nrhs;kii++) _coef1[kii] = czero;

		BlockDotH (isize, _corr.ldcorrarr[icorr], _corr.yarr[icorr], _ldcoef, _coef1);
		BlockDotH (isize, _corr.ldcorrarr[icorr], _corr.warr[icorr], _ldcoef, _coef1+isize);

// Perform multiprocessor dot operation if necessary

		if (_tree.nproc != 1) {
			CMPIComm pcomm = _tree.GetComm ();
			CMPIExchange::ExchangeArrayMPI (pcomm,
														DOUBLEVALUE, ADD,
														_ldcoef*nrhs*2, _coef1, _coef2);
		} else {
			for (kii=0;kii<_ldcoef*nrhs;kii++) _coef2[kii] = _coef1[kii];
		};

// Multiply by the block diagonal

		for (kii=0;kii<_ldcoef*nrhs;kii++) _coef1[kii] = czero;

		for (kii=0;kii<isize*2;kii++) {
			for (kjj=0;kjj<nrhs;kjj++) {
				dcmplx caux = czero;
				for (kkk=0;kkk<isize*2;kkk++) {
					caux += _corr.darr[icorr][kkk*isize*2+kii] * _coef2[kjj*_ldcoef+kkk];
				};
				_coef1[kjj*_ldcoef+kii] = caux;
			};
		};

// Daxpy operations

		BlockDaxpy (isize, _corr.ldcorrarr[icorr], _corr.yarr[icorr], _ldcoef, _coef1);
		BlockDaxpy (isize, _corr.ldcorrarr[icorr], _corr.warr[icorr], _ldcoef, _coef1+isize);

	};

ExitMvCorr:;

};

// Author: Kharchenko S.A.
// CMvmSchur: Implement abstract block Mvm with Schur complement parallelization
//========================================================================================
void CMvmSchur::AbstractBlockMvmA (char _transp, // Implement abstract block Mvm with Schur complement parallelization
												CSVectorC &_x, CSVectorC &_ax) {

//	const char *funcname = "AbstractBlockMvmA";

// Exchange necessary x data

	ExchangeX (_x);

// Init necessary Ax by zeroes

	InitAxByZeroes ();

// Multiply

	PerformLocalMultiplications (_transp);

// Exchange ax and store

	ExchangeAX (_ax);

};

// Author: Kharchenko S.A.
// CMvmSchur: Implement abstract block L solve with Schur complement parallelization
//========================================================================================
void CMvmSchur::AbstractBlockSolveL (char _lutype, // Implement abstract block L solve with Schur complement parallelization
													CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
													CSVectorC &_x, CSVectorC &_px) {

	const char *funcname = "AbstractBlockSolveL";

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Open the tree structure

	int nlevtot = _tree.GetNlev ();
//	int nnodestot = _tree.GetNnodes ();
	int *pcpuidend = _tree.GetCpuidend ();
	CNode *pnodes = _tree.GetNode ();

// Open the maximal factorization block sparsity structure

	int *piablksp = _ablksp.GetIa ();
	int *pjablksp = _ablksp.GetJa ();

// Get description of the tree for the current cpu

	int *i_lev_arr;

	i_lev_arr = new int [nlevtot];
	if (!i_lev_arr) MemoryFail (funcname);

	int ilev, ifather, inodeloc;

	inodeloc = pcpuidend[myid];

	while (true) {
		ilev = pnodes[inodeloc].GetNodelv ();
		i_lev_arr[ilev] = inodeloc;
		ifather = pnodes[inodeloc].FatherId ();
		if (ifather == inodeloc) break;
		inodeloc = ifather;
	};

// Allocate the work arrays

	int *iagroup;
	int *grp2cpu;
	int *grptype;
	int *iagrp2cpu;
	int *jagrp2cpu;
	int *imaskcpu;
	int *listcpu;
	int *iaset2grp;
	int *listfct;
	int *listschur;
	int *imaskschur;

	iagroup = new int [_nblks+1];
	if (!iagroup) MemoryFail (funcname);
	grp2cpu = new int [_nblks];
	if (!grp2cpu) MemoryFail (funcname);
	grptype = new int [_nblks];
	if (!grptype) MemoryFail (funcname);
	iagrp2cpu = new int [_nblks+1];
	if (!iagrp2cpu) MemoryFail (funcname);
	jagrp2cpu = new int [_nblks*nproc];
	if (!jagrp2cpu) MemoryFail (funcname);
	imaskcpu = new int [nproc];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);
	iaset2grp = new int [_nblks+1];
	if (!iaset2grp) MemoryFail (funcname);
	listfct = new int [_nblks];
	if (!listfct) MemoryFail (funcname);
	listschur = new int [_nblks];
	if (!listschur) MemoryFail (funcname);
	imaskschur = new int [_nblks];
	if (!imaskschur) MemoryFail (funcname);

	int icyclecpu = -1;

	int i;

	for (i=0;i<nproc;i++) imaskcpu[i] = icyclecpu;

	int icycleschur = -1;

	for (i=0;i<_nblks;i++) imaskschur[i] = icycleschur;

// Init Schur blocks

	int j, iblk, jjblk;

	icycleschur++;

	int nlistschur = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid) {
			if (imaskschur[iblk] != icycleschur) {
				listschur[nlistschur] = iblk;
				nlistschur++;
				imaskschur[iblk] = icycleschur;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskschur[jjblk] != icycleschur) {
						listschur[nlistschur] = jjblk;
						nlistschur++;
						imaskschur[jjblk] = icycleschur;
					};
				};
			};
		};
	};

// Init Schur blocks

	InitPxBlocks (_x);

// Main cycle over the nodes

	int inode, ibegblk, iendblk, ngroups;
	int iproc, jproc, jbegblk, jendblk, icheck, nz, nlistcpu, nsets, igroup;
	int jgroup, jjbegblk, jjendblk, iset, ibeggrp, iendgrp, myidloc;
	int nchildsloc, ichild0, ichild1, iproc0, iproc1, nlistfct;
	int *pchildsloc;

	for (ilev=nlevtot-1;ilev>=0;ilev--) {

		inode = i_lev_arr[ilev];

		ibegblk = pnodes[inode].GetIndbeg ();
		iendblk = pnodes[inode].GetIndend ();

// Create the whole set of continuous groups of blocks (set of consecutive blocks on the same cpu)

		iblk = ibegblk;
		iproc = _blk2cpu[iblk];
		ngroups = 0;
		iagroup[0] = 0;
		grp2cpu[ngroups] = iproc;
		iagroup[ngroups+1] = iblk+1-ibegblk;
		while (iblk < iendblk) {
			iblk++;
			jproc = _blk2cpu[iblk];
			if (jproc == iproc) {
				iagroup[ngroups+1] = iblk+1-ibegblk;
			} else {
				ngroups++;
				iproc = jproc;
				grp2cpu[ngroups] = iproc;
				iagroup[ngroups+1] = iblk+1-ibegblk;
			};
		};
		ngroups++;

// Mark all groups (not used, Schur or compute)

		for (i=0;i<ngroups;i++) {
			iproc = grp2cpu[i];
			if (iproc == myid) {
				grptype[i] = 1;
			} else {
				jbegblk = iagroup[i]+ibegblk;
				jendblk = iagroup[i+1]-1+ibegblk;
				icheck = -1;
				for (iblk=jbegblk;iblk<=jendblk;iblk++) {
					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk <= iblk) {
							if (_blk2cpu[jjblk] == myid) icheck = 1;
						};
					};
				};
				if (icheck > 0) {
					grptype[i] = 0;
				} else {
					grptype[i] = -1;
				};
			};
		};

// For each working group determine the list of cpu's

		iagrp2cpu[0] = 0;
		nz = 0;

		for (i=0;i<ngroups;i++) {
			icyclecpu++;
			jbegblk = iagroup[i]+ibegblk;
			jendblk = iagroup[i+1]-1+ibegblk;
			iproc = grp2cpu[i];
			listcpu[0] = iproc;
			imaskcpu[iproc] = icyclecpu;
			nlistcpu = 1;
			for (iblk=jbegblk;iblk<=jendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk <= iblk) {
						jproc = _blk2cpu[jjblk];
						if (imaskcpu[jproc] != icyclecpu) {
							listcpu[nlistcpu] = jproc;
							nlistcpu++;
							imaskcpu[jproc] = icyclecpu;
						};
					};
				};
			};
			if (nlistcpu > 1) qsort (listcpu+1,nlistcpu-1,sizeof(int),compint);
			for (j=0;j<nlistcpu;j++) jagrp2cpu[nz+j] = listcpu[j];
			nz += nlistcpu;
			iagrp2cpu[i+1] = nz;
		};

// Find the sets of sequential groups that should be processed in parallel

		nsets = 0;

		igroup = 0;
		iaset2grp[0] = 0;
		iaset2grp[nsets+1] = igroup+1;
		jbegblk = iagroup[igroup]+ibegblk;
		jendblk = iagroup[igroup+1]-1+ibegblk;
		while (igroup < ngroups-1) {
			jgroup = igroup+1;
			jjbegblk = iagroup[jgroup]+ibegblk;
			jjendblk = iagroup[jgroup+1]-1+ibegblk;
			icheck = 1;
			for (iblk=jjbegblk;iblk<=jjendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk >= jbegblk && jjblk <= jendblk) icheck = -1;
				};
			};
			if (icheck == 1) {
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			} else {
				nsets++;
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jbegblk = iagroup[igroup]+ibegblk;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			};
		};
		nsets++;

// Main cycle over the sets

		for (iset=0;iset<nsets;iset++) {

			ibeggrp = iaset2grp[iset];
			iendgrp = iaset2grp[iset+1]-1;

// Make binary tree type Schur exchanges and summation for all groups of blocks in the set

			for (igroup=ibeggrp;igroup<=iendgrp;igroup++) {
				if (grptype[igroup] >= 0) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistschur = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listschur[nlistschur] = i;
						nlistschur++;
					};

// Plan sends and receives

					myidloc = -1;

					nlistcpu = 0;
					for (j=iagrp2cpu[igroup];j<iagrp2cpu[igroup+1];j++) {
						jproc = jagrp2cpu[j];
						if (jproc == myid) myidloc = nlistcpu;
						listcpu[nlistcpu] = jproc;
						nlistcpu++;
					};

// Do exchanges and summation according to the plan

					inodeloc = pcpuidend[myidloc];

					while (true) {
						nchildsloc = pnodes[inodeloc].GetNchilds ();
						pchildsloc = pnodes[inodeloc].GetChilds ();
						iproc = pnodes[inodeloc].GetNodecpu ();
						if (iproc != myidloc) break;
						if (nchildsloc > 1) {
							ichild0 = pchildsloc[0];
							ichild1 = pchildsloc[1];
							iproc0 = pnodes[ichild0].GetNodecpu ();
							iproc1 = pnodes[ichild1].GetNodecpu ();
							if (iproc0 == myidloc) {
								if (iproc1 < nlistcpu) {
									jproc = listcpu[iproc1];
									ReceiveBlocks (true, jproc, nlistschur, listschur);
								};
							} else if (iproc1 == myidloc) {
								if (iproc0 < nlistcpu) {
									jproc = listcpu[iproc0];
									ReceiveBlocks (true, jproc, nlistschur, listschur);
								};
							};
						};
						ifather = pnodes[inodeloc].FatherId ();
						if (ifather == inodeloc) break;
						iproc0 = pnodes[ifather].GetNodecpu ();
						if (iproc0 != myidloc) {
							jproc = listcpu[iproc0];
							SendBlocks (1, &jproc, nlistschur, listschur);
						};
						inodeloc = ifather;
					};

				};
			};

// Wait for completion of sends

			WaitSends ();

// Perform solution for own blocks

			for (igroup=ibeggrp;igroup<=iendgrp;igroup++) {

				if (grptype[igroup] == 1) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistfct = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listfct[nlistfct] = i;
						nlistfct++;
					};

					SolveLBlocks (_lutype, nlistfct, listfct);

				};

			};

		};

	};

	StorePx (_px);

// Free work arrays

	delete [] i_lev_arr;
	delete [] iagroup;
	delete [] grp2cpu;
	delete [] grptype;
	delete [] iagrp2cpu;
	delete [] jagrp2cpu;
	delete [] imaskcpu;
	delete [] listcpu;
	delete [] iaset2grp;
	delete [] listfct;
	delete [] listschur;
	delete [] imaskschur;

};

// Author: Kharchenko S.A.
// CMvmSchur: Implement abstract block U solve with Schur complement parallelization
//========================================================================================
void CMvmSchur::AbstractBlockSolveU (char _lutype, // Implement abstract block U solve with Schur complement parallelization
													CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
													CSVectorC &_x, CSVectorC &_px) {

	const char *funcname = "AbstractBlockSolveU";

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Open the tree structure

	int nlevtot = _tree.GetNlev ();
//	int nnodestot = _tree.GetNnodes ();
	int *pcpuidend = _tree.GetCpuidend ();
	CNode *pnodes = _tree.GetNode ();

// Open the maximal factorization block sparsity structure

	int *piablksp = _ablksp.GetIa ();
	int *pjablksp = _ablksp.GetJa ();

// Get description of the tree for the current cpu

	int *i_lev_arr;

	i_lev_arr = new int [nlevtot];
	if (!i_lev_arr) MemoryFail (funcname);

	int ilev, ifather, inodeloc;

	inodeloc = pcpuidend[myid];

	while (true) {
		ilev = pnodes[inodeloc].GetNodelv ();
		i_lev_arr[ilev] = inodeloc;
		ifather = pnodes[inodeloc].FatherId ();
		if (ifather == inodeloc) break;
		inodeloc = ifather;
	};

// Allocate the work arrays

	int *iagroup;
	int *grp2cpu;
	int *grptype;
	int *iagrp2cpu;
	int *jagrp2cpu;
	int *imaskcpu;
	int *listcpu;
	int *iaset2grp;
	int *listfct;
	int *listschur;
	int *imaskschur;

	iagroup = new int [_nblks+1];
	if (!iagroup) MemoryFail (funcname);
	grp2cpu = new int [_nblks];
	if (!grp2cpu) MemoryFail (funcname);
	grptype = new int [_nblks];
	if (!grptype) MemoryFail (funcname);
	iagrp2cpu = new int [_nblks+1];
	if (!iagrp2cpu) MemoryFail (funcname);
	jagrp2cpu = new int [_nblks*nproc];
	if (!jagrp2cpu) MemoryFail (funcname);
	imaskcpu = new int [nproc];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);
	iaset2grp = new int [_nblks+1];
	if (!iaset2grp) MemoryFail (funcname);
	listfct = new int [_nblks];
	if (!listfct) MemoryFail (funcname);
	listschur = new int [_nblks];
	if (!listschur) MemoryFail (funcname);
	imaskschur = new int [_nblks];
	if (!imaskschur) MemoryFail (funcname);

	int icyclecpu = -1;

	int i;

	for (i=0;i<nproc;i++) imaskcpu[i] = icyclecpu;

	int icycleschur = -1;

	for (i=0;i<_nblks;i++) imaskschur[i] = icycleschur;

// Init Schur blocks

	int j, iblk, jjblk;

	icycleschur++;

	int nlistschur = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid) {
			if (imaskschur[iblk] != icycleschur) {
				listschur[nlistschur] = iblk;
				nlistschur++;
				imaskschur[iblk] = icycleschur;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskschur[jjblk] != icycleschur) {
						listschur[nlistschur] = jjblk;
						nlistschur++;
						imaskschur[jjblk] = icycleschur;
					};
				};
			};
		};
	};

// Init Schur blocks

	InitPxBlocks (_x);

// Main cycle over the nodes

	int inode, ibegblk, iendblk, ngroups;
	int iproc, jproc, jbegblk, jendblk, icheck, nz, nlistcpu, nsets, igroup;
	int jgroup, jjbegblk, jjendblk, iset, ibeggrp, iendgrp;
	int nlistfct;

	for (ilev=0;ilev<nlevtot;ilev++) {

		inode = i_lev_arr[ilev];

		ibegblk = pnodes[inode].GetIndbeg ();
		iendblk = pnodes[inode].GetIndend ();

// Create the whole set of continuous groups of blocks (set of consecutive blocks on the same cpu)

		iblk = ibegblk;
		iproc = _blk2cpu[iblk];
		ngroups = 0;
		iagroup[0] = 0;
		grp2cpu[ngroups] = iproc;
		iagroup[ngroups+1] = iblk+1-ibegblk;
		while (iblk < iendblk) {
			iblk++;
			jproc = _blk2cpu[iblk];
			if (jproc == iproc) {
				iagroup[ngroups+1] = iblk+1-ibegblk;
			} else {
				ngroups++;
				iproc = jproc;
				grp2cpu[ngroups] = iproc;
				iagroup[ngroups+1] = iblk+1-ibegblk;
			};
		};
		ngroups++;

// Mark all groups (not used, Schur or compute)

		for (i=0;i<ngroups;i++) {
			iproc = grp2cpu[i];
			if (iproc == myid) {
				grptype[i] = 1;
			} else {
				jbegblk = iagroup[i]+ibegblk;
				jendblk = iagroup[i+1]-1+ibegblk;
				icheck = -1;
				for (iblk=jbegblk;iblk<=jendblk;iblk++) {
					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk <= iblk) {
							if (_blk2cpu[jjblk] == myid) icheck = 1;
						};
					};
				};
				if (icheck > 0) {
					grptype[i] = 0;
				} else {
					grptype[i] = -1;
				};
			};
		};

// For each working group determine the list of cpu's

		iagrp2cpu[0] = 0;
		nz = 0;

		for (i=0;i<ngroups;i++) {
			icyclecpu++;
			jbegblk = iagroup[i]+ibegblk;
			jendblk = iagroup[i+1]-1+ibegblk;
			iproc = grp2cpu[i];
			listcpu[0] = iproc;
			imaskcpu[iproc] = icyclecpu;
			nlistcpu = 1;
			for (iblk=jbegblk;iblk<=jendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk <= iblk) {
						jproc = _blk2cpu[jjblk];
						if (imaskcpu[jproc] != icyclecpu) {
							listcpu[nlistcpu] = jproc;
							nlistcpu++;
							imaskcpu[jproc] = icyclecpu;
						};
					};
				};
			};
			if (nlistcpu > 1) qsort (listcpu+1,nlistcpu-1,sizeof(int),compint);
			for (j=0;j<nlistcpu;j++) jagrp2cpu[nz+j] = listcpu[j];
			nz += nlistcpu;
			iagrp2cpu[i+1] = nz;
		};

// Find the sets of sequential groups that should be processed in parallel

		nsets = 0;

		igroup = 0;
		iaset2grp[0] = 0;
		iaset2grp[nsets+1] = igroup+1;
		jbegblk = iagroup[igroup]+ibegblk;
		jendblk = iagroup[igroup+1]-1+ibegblk;
		while (igroup < ngroups-1) {
			jgroup = igroup+1;
			jjbegblk = iagroup[jgroup]+ibegblk;
			jjendblk = iagroup[jgroup+1]-1+ibegblk;
			icheck = 1;
			for (iblk=jjbegblk;iblk<=jjendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk >= jbegblk && jjblk <= jendblk) icheck = -1;
				};
			};
			if (icheck == 1) {
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			} else {
				nsets++;
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jbegblk = iagroup[igroup]+ibegblk;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			};
		};
		nsets++;

// Main cycle over the sets

		for (iset=nsets-1;iset>=0;iset--) {

			ibeggrp = iaset2grp[iset];
			iendgrp = iaset2grp[iset+1]-1;

// Perform solution for own blocks

			for (igroup=iendgrp;igroup>=ibeggrp;igroup--) {

				if (grptype[igroup] == 1) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistfct = 0;
					for (i=jendblk;i>=jbegblk;i--) {
						listfct[nlistfct] = i;
						nlistfct++;
					};

					SolveUBlocks (_lutype, nlistfct, listfct);

				};

			};

// Make solution exchanges for all groups of blocks in the set

			for (igroup=iendgrp;igroup>=ibeggrp;igroup--) {
				if (grptype[igroup] >= 0) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistschur = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listschur[nlistschur] = i;
						nlistschur++;
					};

// Plan sends and receives

					nlistcpu = 0;
					for (j=iagrp2cpu[igroup];j<iagrp2cpu[igroup+1];j++) {
						jproc = jagrp2cpu[j];
						listcpu[nlistcpu] = jproc;
						nlistcpu++;
					};

					if (nlistcpu > 1) {

						if (grptype[igroup] == 1) {
							SendBlocks (nlistcpu-1, listcpu+1, nlistschur, listschur);
						} else {
							ReceiveBlocks (false, listcpu[0], nlistschur, listschur);
						};

					};
				};
			};

// Wait for completion of sends

			WaitSends ();

		};

	};

	StorePx (_px);

// Free work arrays

	delete [] i_lev_arr;
	delete [] iagroup;
	delete [] grp2cpu;
	delete [] grptype;
	delete [] iagrp2cpu;
	delete [] jagrp2cpu;
	delete [] imaskcpu;
	delete [] listcpu;
	delete [] iaset2grp;
	delete [] listfct;
	delete [] listschur;
	delete [] imaskschur;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CMvmSchurC::CMvmSchurC()
//========================================================================================
CMvmSchurC::CMvmSchurC (): CMvmSchur () { // Memory allocation zero data constructor

	const char *funcname = "CMvmSchurC_00";

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [0];
	if (!ibsvect) MemoryFail (funcname);
	listblk = new int [0];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [0];
	if (!listblkext) MemoryFail (funcname);
	iasndblk = new int [0];
	if (!iasndblk) MemoryFail (funcname);
	jasndblk = new int [0];
	if (!jasndblk) MemoryFail (funcname);
	iarcvblk = new int [0];
	if (!iarcvblk) MemoryFail (funcname);
	jarcvblk = new int [0];
	if (!jarcvblk) MemoryFail (funcname);

	nlistblk = 0;
	nlistblkext = 0;

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	pvect = new CSVectorC;
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVectorC;
	if (!pvect1) MemoryFail (funcname);

	PrepareWorkData ();

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CMvmSchurC::CMvmSchurC()
//========================================================================================
CMvmSchurC::CMvmSchurC (std::ofstream *_pfout, // Memory allocation zero data constructor
								int _nblks, int *_blks, int *_bl2nd, int *_blk2cpu,
								CTree *_ptree, CSMatrix *_pablkstr, 
								void *_pgmtral, void *_pgmtrau, 
								void *_pgmtrl, void *_pgmtru,
								FUNC_MULTC *_pmvmal, FUNC_MULTC *_pmvmau,
								FUNC_SOLVEC *_psolvel, FUNC_MULTC *_psolveu,
								CSlvParam *_pparam
								): CMvmSchur () {

	const char *funcname = "CMvmSchurR_01";

	pfout = _pfout;
	nblks = _nblks;
	pblks = _blks;
	pbl2nd = _bl2nd;
	pblk2cpu = _blk2cpu;
	ptree = _ptree;
	pablkstr = _pablkstr;
	pgmtral = _pgmtral;
	pgmtrau = _pgmtrau;
	pgmtrl = _pgmtrl;
	pgmtru = _pgmtru;
	pmvmal = _pmvmal;
	pmvmau = _pmvmau;
	psolvel = _psolvel;
	psolveu = _psolveu;

	pparam = _pparam;

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [0];
	if (!ibsvect) MemoryFail (funcname);
	listblk = new int [0];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [0];
	if (!listblkext) MemoryFail (funcname);
	iasndblk = new int [0];
	if (!iasndblk) MemoryFail (funcname);
	jasndblk = new int [0];
	if (!jasndblk) MemoryFail (funcname);
	iarcvblk = new int [0];
	if (!iarcvblk) MemoryFail (funcname);
	jarcvblk = new int [0];
	if (!jarcvblk) MemoryFail (funcname);

	nlistblk = 0;
	nlistblkext = 0;

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	pvect = new CSVectorC;
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVectorC;
	if (!pvect1) MemoryFail (funcname);

	PrepareWorkData ();

};

// Author: Kharchenko S.A.
// Description: Destructor
// CMvmSchurC::~CMvmSchurC()
//========================================================================================
CMvmSchurC::~CMvmSchurC () { // Destructor
//	cout << " On entry to CMvmSchurC destructor " << endl;

	DestroyWorkData ();

	delete [] imaskflt;
	delete [] ibsvect;
	delete [] listblk;
	delete [] listblkext;
	delete [] iasndblk;
	delete [] jasndblk;
	delete [] iarcvblk;
	delete [] jarcvblk;

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

	delete pvect;
	delete pvect1;

//	cout << " On return from CMvmSchurC destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Prepare solve working data
// CMvmSchurC::PrepareWorkData()
//========================================================================================
void CMvmSchurC::PrepareWorkData () { // Prepare solve working data

	const char *funcname = "PrepareWorkData";

// Open the tree

	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Open the sparsity

	int *piablksp = pablkstr->GetIa ();
	int *pjablksp = pablkstr->GetJa ();

// Init the mask array

	delete [] imaskflt;
	delete [] ibsvect;

	imaskflt = new int [nblks];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [nblks];
	if (!ibsvect) MemoryFail (funcname);

	int i;
	icycleflt = -1;

	for (i=0;i<nblks;i++) imaskflt[i] = icycleflt;
	for (i=0;i<nblks;i++) ibsvect[i] = -1;

// Compute the maximal list of blocks to be stored on current cpu

	int *listloc;

	listloc = new int [nblks];
	if (!listloc) MemoryFail (funcname);

	int j, iblk, jjblk, k, kkblk;

	icycleflt++;

	int nlistschur = 0;
	nlistblk = 0;

	for (iblk=0;iblk<nblks;iblk++) {
		if (pblk2cpu[iblk] == myid) {
			nlistblk++;
			if (imaskflt[iblk] != icycleflt) {
				listloc[nlistschur] = iblk;
				nlistschur++;
				imaskflt[iblk] = icycleflt;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskflt[jjblk] != icycleflt) {
						listloc[nlistschur] = jjblk;
						nlistschur++;
						imaskflt[jjblk] = icycleflt;
					};
				} else if (jjblk < iblk) {
					for (k=piablksp[jjblk];k<piablksp[jjblk+1];k++) {
						kkblk = pjablksp[k];
						if (kkblk >= iblk) {
							if (imaskflt[kkblk] != icycleflt) {
								listloc[nlistschur] = kkblk;
								nlistschur++;
								imaskflt[kkblk] = icycleflt;
							};
						};
					};
				};
			};
		};
	};

// Add tree data

	int *pcpuidend = ptree->GetCpuidend ();
	CNode *pnodes = ptree->GetNode ();

	int inode = pcpuidend[myid];

	int indbeg, indend, ifather;

	while (true) {
		indbeg = pnodes[inode].GetIndbeg ();
		indend = pnodes[inode].GetIndend ();
		for (i=indbeg;i<=indend;i++) {
			if (imaskflt[i] != icycleflt) {
				listloc[nlistschur] = i;
				nlistschur++;
				imaskflt[i] = icycleflt;
			};
		};
		ifather = pnodes[inode].FatherId ();
		if (ifather == inode) break;
		inode = ifather;
	};

// Sort the list

	qsort (listloc, nlistschur, sizeof(int), compint);

	delete [] listblk;
	delete [] listblkext;

	nlistblkext = nlistschur;

	listblk = new int [nlistblk];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [nlistblkext];
	if (!listblkext) MemoryFail (funcname);

	for (i=0;i<nlistblkext;i++) listblkext[i] = listloc[i];

	nlistblk = 0;

	for (iblk=0;iblk<nblks;iblk++) {
		if (pblk2cpu[iblk] == myid) {
			listblk[nlistblk] = iblk;
			nlistblk++;
		};
	};

// Create sends/receives list

	delete [] iasndblk;
	delete [] jasndblk;

	iasndblk = new int [nproc+1];
	if (!iasndblk) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) {
		iasndblk[i] = 0;
	};

	int iproc, jproc, ni;

	for (iproc=0;iproc<nproc;iproc++) {

		if (iproc != myid) {

			icycleflt++;

			for (iblk=0;iblk<nblks;iblk++) {
				if (pblk2cpu[iblk] == iproc) {

					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
							jproc = pblk2cpu[jjblk];
							if (jproc == myid) {
								iasndblk[iproc+1]++;
							};
							imaskflt[jjblk] = icycleflt;
						};
					};

				};
			};

		};

	};

	for (i=0;i<nproc;i++) iasndblk[i+1] = iasndblk[i]+iasndblk[i+1];

	int nzjasnd = iasndblk[nproc];

	jasndblk = new int [nzjasnd];
	if (!jasndblk) MemoryFail (funcname);

	int *iptr;

	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) iptr[i] = iasndblk[i];

	for (iproc=0;iproc<nproc;iproc++) {

		if (iproc != myid) {

			icycleflt++;

			for (iblk=0;iblk<nblks;iblk++) {
				if (pblk2cpu[iblk] == iproc) {

					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
							jproc = pblk2cpu[jjblk];
							if (jproc == myid) {
								k = iptr[iproc];
								jasndblk[k] = jjblk;
								iptr[iproc]++;
							};
							imaskflt[jjblk] = icycleflt;
						};
					};

				};
			};

		};
	};

	for (i=0;i<nproc;i++) {
		ni = iasndblk[i+1]-iasndblk[i];
		if (ni != 0) qsort (jasndblk+iasndblk[i], ni, sizeof(int), compint);
	};

	delete [] iarcvblk;
	delete [] jarcvblk;

	iarcvblk = new int [nproc+1];
	if (!iarcvblk) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) {
		iarcvblk[i] = 0;
	};

	icycleflt++;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
			jjblk = pjablksp[j];
			if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
				jproc = pblk2cpu[jjblk];
				if (jproc != myid) {
					iarcvblk[jproc+1]++;
				};
				imaskflt[jjblk] = icycleflt;
			};
		};
	};

	for (i=0;i<nproc;i++) iarcvblk[i+1] = iarcvblk[i]+iarcvblk[i+1];

	int nzjarcv= iarcvblk[nproc];

	jarcvblk = new int [nzjarcv];
	if (!jarcvblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) iptr[i] = iarcvblk[i];

	icycleflt++;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
			jjblk = pjablksp[j];
			if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
				jproc = pblk2cpu[jjblk];
				if (jproc != myid) {
					k = iptr[jproc];
					jarcvblk[k] = jjblk;
					iptr[jproc]++;
				};
				imaskflt[jjblk] = icycleflt;
			};
		};
	};

	for (i=0;i<nproc;i++) {
		ni = iarcvblk[i+1]-iarcvblk[i];
		if (ni != 0) qsort (jarcvblk+iarcvblk[i], ni, sizeof(int), compint);
	};

// Create the maximal vector work array

	int nloc = 0;

	for (i=0;i<nlistschur;i++) {
		iblk = listloc[i];
		ni = pbl2nd[iblk+1]-pbl2nd[iblk];
		ibsvect[iblk] = nloc;
		nloc += ni;
	};

	delete pvect;
	delete pvect1;

	pvect = new CSVectorC (nloc);
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVectorC (nloc);
	if (!pvect1) MemoryFail (funcname);

// Allocate send receive work arrays

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

// Allocate send/receive structures

	int nblkmax = 0;

	int nprocloc = ptree->GetNproc ();
	int nnodesloc = ptree->GetNnodes ();

	int iblkbeg, iblkend;

	for (i=0;i<nnodesloc;i++) {
		iblkbeg = pnodes[i].GetIndbeg ();
		iblkend = pnodes[i].GetIndend ();
		ni = iblkend-iblkbeg+1;
		if (ni>nblkmax) nblkmax = ni;
	};

	nsendsmax = nprocloc*nblkmax+10;
	nsends = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	for (i=0;i<nsendsmax;i++) imasksend[i] = -1;

// Free work arrays

	delete [] listloc;
	delete [] iptr;

};

// Author: Kharchenko S.A.
// Description: Exchange X data
// CMvmSchurC::ExchangeX()
//========================================================================================
void CMvmSchurC::ExchangeX (CSVectorC &_x) { // Exchange X data

	const char *funcname = "ExchangeX";

//	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Init X part

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *px = _x.GetVect ();
	dcmplx *pdvect = pvect->GetVect ();

	pvect->SetSVect (czero);

	int i, iblk, ni, ibs, j;

	int nloc = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pbl2nd[iblk+1]-pbl2nd[iblk];
		for (j=0;j<ni;j++) pdvect[ibs+j] = px[nloc+j];
		nloc += ni;
	};

// Prepare sends and receives

	int *psize;
	dcmplx **psendarr;

	psize = new int [nproc];
	if (!psize) MemoryFail (funcname);
	psendarr = new dcmplx * [nproc];
	if (!psendarr) MemoryFail (funcname);

	int jjblk, niloc, k;

	for (i=0;i<nproc;i++) {
		ni = 0;
		for (j=iasndblk[i];j<iasndblk[i+1];j++) {
			jjblk = jasndblk[j];
			ni += pbl2nd[jjblk+1]-pbl2nd[jjblk];
		};
		psendarr[i] = new dcmplx [ni];
		if (!psendarr[i]) MemoryFail (funcname);
		ni = 0;
		for (j=iasndblk[i];j<iasndblk[i+1];j++) {
			jjblk = jasndblk[j];
			ibs = ibsvect[jjblk];
			niloc = pbl2nd[jjblk+1]-pbl2nd[jjblk];
			for (k=0;k<niloc;k++) psendarr[i][ni+k] = pdvect[ibs+k];
			ni += niloc;
		};
		psize[i] = ni * sizeof(dcmplx);
	};

// Send/receive all necessary data

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
		ObjSizeSend[i] = psize[i];
		ObjSend[i] = (char *)psendarr[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (ptree->GetComm (),
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CMvmSchurR::ExchangeX: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] psendarr[i];

// Store the result

	int jproc;
	dcmplx *pdarr;

	for (i=0;i<NObjRecv;i++) {
		jproc = CpuIDRecv[i];
		pdarr = (dcmplx *)ObjRecv[i];
		ni = 0;
		for (j=iarcvblk[jproc];j<iarcvblk[jproc+1];j++) {
			jjblk = jarcvblk[j];
			ibs = ibsvect[jjblk];
			niloc = pbl2nd[jjblk+1]-pbl2nd[jjblk];
			for (k=0;k<niloc;k++) pdvect[ibs+k] = pdarr[ni+k];
			ni += niloc;
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

	delete [] psize;
	delete [] psendarr;

};

// Author: Kharchenko S.A.
// Description: Init necessary Ax by zeroes
// CMvmSchurC::InitAxByZeroes()
//========================================================================================
void CMvmSchurC::InitAxByZeroes () { // Init necessary Ax by zeroes

	dcmplx czero (0.0e0,0.0e0);

	pvect1->SetSVect (czero);

};

// Author: Kharchenko S.A.
// Description: Multiply by local block rows/columns of A
// CMvmSchurC::PerformLocalMultiplications()
//========================================================================================
void CMvmSchurC::PerformLocalMultiplications (char _transp) { // Multiply by local block rows/columns of A

//	int myid = ptree->GetMyid ();

	void *pobj = GetSchurmvm ();

	int i, iblk;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		if (_transp == 'N' || _transp == 'n') {
			(*pmvmau) (pobj, pgmtrau, iblk, ibsvect, *pvect, ibsvect, *pvect1);
			(*pmvmal) (pobj, pgmtral, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else {
			(*pmvmau) (pobj, pgmtral, iblk, ibsvect, *pvect, ibsvect, *pvect1);
			(*pmvmal) (pobj, pgmtrau, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		};
	};

};

// Author: Kharchenko S.A.
// Description: Multiply by the current block row of Au
// CParSchurMvmSlvCS::MvmBlockRowAu()
//========================================================================================
void CParSchurMvmSlvCS::AbsMvmBlockRowAu (void *_obj, void *_mtr, int _iblk, // Multiply by the current block row of Au
														int *_bsx, const CSVectorC &_x, 
														int *_bsax, CSVectorC &_ax) {

	CParSchurMvmSlvCS *pobj = (CParSchurMvmSlvCS *)_obj;
	CGSMatrixCS *pgmtrauloc = (CGSMatrixCS *)_mtr;

	int *ibssuploc = pobj->GetIbssup ();

	int *psprndsloc = pgmtrauloc->GetSprndr ();

	CSMatrixCS *pmtrarralu = ((CGSMatrixCS *)pgmtrauloc)->GetMtrarr ();

	pmtrarralu[_iblk].MvmBlockRowAu (psprndsloc, ibssuploc, _x, ibssuploc, _ax);

};

// Author: Kharchenko S.A.
// Description: Multiply by the current block column of Al without diagonal supernode
// CParSchurMvmSlvCS::MvmBlockColumnAl()
//========================================================================================
void CParSchurMvmSlvCS::AbsMvmBlockColumnAl (void *_obj, void *_mtr, int _iblk, // Multiply by the current block column of Al without diagonal supernodeL
								int *_bsx, const CSVectorC &_x, 
								int *_bsax, CSVectorC &_ax) {

	CParSchurMvmSlvCS *pobj = (CParSchurMvmSlvCS *)_obj;
	CGSMatrixCS *pgmtralloc = (CGSMatrixCS *)_mtr;

	int *ibssuploc = pobj->GetIbssup ();

	int *psprndsloc = pgmtralloc->GetSprndr ();

	CSMatrixCS *pmtrarralu = ((CGSMatrixCS *)pgmtralloc)->GetMtrarr ();

	pmtrarralu[_iblk].MvmBlockColumnAl (psprndsloc, ibssuploc, _x, ibssuploc, _ax);

};

// Author: Kharchenko S.A.
// Description: Exchange X data
// CMvmSchurC::ExchangeAX()
//========================================================================================
void CMvmSchurC::ExchangeAX (CSVectorC &_ax) { // Exchange AX data

	const char *funcname = "ExchangeAX";

//	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Prepare sends and receives

	dcmplx *pdvect = pvect1->GetVect ();

	int *psize;
	dcmplx **psendarr;

	psize = new int [nproc];
	if (!psize) MemoryFail (funcname);
	psendarr = new dcmplx * [nproc];
	if (!psendarr) MemoryFail (funcname);

	int i, iblk, ni, ibs, j;
	int jjblk, niloc, k;

	for (i=0;i<nproc;i++) {
		ni = 0;
		for (j=iarcvblk[i];j<iarcvblk[i+1];j++) {
			jjblk = jarcvblk[j];
			ni += pbl2nd[jjblk+1]-pbl2nd[jjblk];
		};
		psendarr[i] = new dcmplx [ni];
		if (!psendarr[i]) MemoryFail (funcname);
		ni = 0;
		for (j=iarcvblk[i];j<iarcvblk[i+1];j++) {
			jjblk = jarcvblk[j];
			ibs = ibsvect[jjblk];
			niloc = pbl2nd[jjblk+1]-pbl2nd[jjblk];
			for (k=0;k<niloc;k++) psendarr[i][ni+k] = pdvect[ibs+k];
			ni += niloc;
		};
		psize[i] = ni * sizeof(dcmplx);
	};

// Send/receive all necessary data

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
		ObjSizeSend[i] = psize[i];
		ObjSend[i] = (char *)psendarr[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (ptree->GetComm (),
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CMvmSchurR::ExchangeX: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] psendarr[i];

// Store the result

	int jproc;
	dcmplx *pdarr;

	for (i=0;i<NObjRecv;i++) {
		jproc = CpuIDRecv[i];
		pdarr = (dcmplx *)ObjRecv[i];
		ni = 0;
		for (j=iasndblk[jproc];j<iasndblk[jproc+1];j++) {
			jjblk = jasndblk[j];
			ibs = ibsvect[jjblk];
			niloc = pbl2nd[jjblk+1]-pbl2nd[jjblk];
			for (k=0;k<niloc;k++) pdvect[ibs+k] += pdarr[ni+k];
			ni += niloc;
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

	delete [] psize;
	delete [] psendarr;

// Store AX part

	dcmplx czero (0.0e0,0.0e0);

	pvect->SetSVect (czero);

	int nloc = 0;

	dcmplx *pax = _ax.GetVect ();

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pbl2nd[iblk+1]-pbl2nd[iblk];
		for (j=0;j<ni;j++) pax[nloc+j] = pdvect[ibs+j];
		nloc += ni;
	};

};

// Author: Kharchenko S.A.
// Description: Init the set of Px blocks
// CMvmSchurC::InitPxBlocks()
//========================================================================================
void CMvmSchurC::InitPxBlocks (CSVectorC &_x) { // Init the set of Px blocks

//	const char *funcname = "InitPxBlocks";

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Init X part

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *px = _x.GetVect ();
	dcmplx *pdvect = pvect->GetVect ();

	pvect->SetSVect (czero);
	pvect1->SetSVect (czero);

	int i, iblk, ni, ibs, j;

	int nloc = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pbl2nd[iblk+1]-pbl2nd[iblk];
		for (j=0;j<ni;j++) pdvect[ibs+j] = px[nloc+j];
		nloc += ni;
	};

};

// Author: Kharchenko S.A.
// Description: Solve L with the set of blocks
// CMvmSchurC::SolveLBlocks()
//========================================================================================
void CMvmSchurC::SolveLBlocks (char _lutype, int _nlistfct, int *_listfct) { // Solve L with the set of blocks

//	int myid = ptree->GetMyid ();

	void *pobj = GetSchurmvm ();

	int i, iblk;

	for (i=0;i<_nlistfct;i++) {
		iblk = _listfct[i];
		if (_lutype == 'L' || _lutype == 'l') {
			(*psolvel) (pobj, pgmtrl, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else if (_lutype == 'A' || _lutype == 'a') {
			(*psolvel) (pobj, pgmtral, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else if (_lutype == 'U' || _lutype == 'u') {
			(*psolvel) (pobj, pgmtru, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		};
	};

};

// Author: Kharchenko S.A.
// Description: Solve triangular system with block column of L
// CParSchurMvmSlvCS::AbsSolveBlockColumnL()
//========================================================================================
void CParSchurMvmSlvCS::AbsSolveBlockColumnL (void *_obj, void *_mtr, int _iblk, // Solve triangular system with block column of L
								int *_bsx, const CSVectorC &_x, 
								int *_bslx, CSVectorC &_lx) {

	CParSchurMvmSlvCS *pobj = (CParSchurMvmSlvCS *)_obj;
	CGSMatrixCS *pgmtrlloc = (CGSMatrixCS *)_mtr;

	int *ibssuploc = pobj->GetIbssup ();

	int *psprndsloc = pgmtrlloc->GetSprndr ();

	CSMatrixCS *pmtrarrlu = ((CGSMatrixCS *)pgmtrlloc)->GetMtrarr ();

	pmtrarrlu[_iblk].SolveBlockColumnL (psprndsloc, ibssuploc, _x, ibssuploc, _lx);

};

// Author: Kharchenko S.A.
// Description: Solve U with the set of blocks
// CMvmSchurC::SolveUBlocks()
//========================================================================================
void CMvmSchurC::SolveUBlocks (char _lutype, int _nlistfct, int *_listfct) { // Solve U with the set of blocks

//	int myid = ptree->GetMyid ();

	void *pobj = GetSchurmvm ();

	int i, iblk;

	for (i=0;i<_nlistfct;i++) {
		iblk = _listfct[i];
		if (_lutype == 'L' || _lutype == 'l') {
			(*psolveu) (pobj, pgmtrl, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else {
			(*psolveu) (pobj, pgmtru, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		};
	};

};

// Author: Kharchenko S.A.
// Description: Solve triangular system with the block row of U
// CParSchurMvmSlvCS::SolveBlockRowU()
//========================================================================================
void CParSchurMvmSlvCS::AbsSolveBlockRowU (void *_obj, void *_mtr, int _iblk, // Solve triangular system with the block row of U
								int *_bsx, const CSVectorC &_x, 
								int *_bsux, CSVectorC &_ux) {

	CParSchurMvmSlvCS *pobj = (CParSchurMvmSlvCS *)_obj;
	CGSMatrixCS *pgmtruloc = (CGSMatrixCS *)_mtr;

	int *ibssuploc = pobj->GetIbssup ();

	int *psprndsloc = pgmtruloc->GetSprndr ();

	CSMatrixCS *pmtrarrlu = ((CGSMatrixCS *)pgmtruloc)->GetMtrarr ();

	pmtrarrlu[_iblk].SolveBlockRowU (psprndsloc, ibssuploc, _x, ibssuploc, _ux);

};

// Author: Kharchenko S.A.
// Description: Send the set of Px blocks
// CMvmSchurC::SendBlocks()
//========================================================================================
void CMvmSchurC::SendBlocks (int _nlistcpu, int *_listcpu, int _nlistschur, int *_listschur) { // Send the set of Px blocks

	const char *funcname = "SendBlocks";

	if (_nlistschur == 0) return;

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Prepare the send buffer

	int i, iblk, niloc;

	int ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		niloc = pbl2nd[iblk+1]-pbl2nd[iblk];
		ni += niloc;
	};

	dcmplx *darr;

	darr = new dcmplx [ni];
	if (!darr) MemoryFail (funcname);

	dcmplx *pdvect = pvect1->GetVect ();

	int ibs, j;

	ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		ibs = ibsvect[iblk];
		niloc = pbl2nd[iblk+1]-pbl2nd[iblk];
		for (j=0;j<niloc;j++) darr[ni+j] = pdvect[ibs+j];
		ni += niloc;
	};
//	OutArr (*pfout," Send to CPUs ",_nlistcpu,_listcpu);
//	OutArr (*pfout," Send list ",_nlistschur,_listschur);
//	OutArr (*pfout," IbsVect ",nblks,ibsvect);
//	OutArr (*pfout," Blks ",nblks+1,pblks);
//	OutArr (*pfout," Send arr ",ni,darr);

// Send the set of blocks to the list of cpu's

	int jproc;

	for (i=0;i<_nlistcpu;i++) {
		jproc = _listcpu[i];

		if (nsends >= nsendsmax) throw " CMvmSchurR::SendBlocks: insufficient number of sends ";

		CMPIExchange::ISend (ptree->GetComm (), jproc, 1,
									ni*sizeof(dcmplx), (char *)darr, sendreqvarr[nsends]);

		imasksend[nsends] = 1;
		if (i==0) {
			psenddata[nsends] = (char *)darr;
		} else {
			psenddata[nsends] = 0;
		};
		nsends++;
	};

};

// Author: Kharchenko S.A.
// Description: Receive the set of Px blocks
// CMvmSchurC::ReceiveBlocks()
//========================================================================================
void CMvmSchurC::ReceiveBlocks (bool _add, int _jproc, int _nlistschur, int *_listschur) { // Receive the set of Px blocks

	const char *funcname = "ReceiveBlocks";

	if (_nlistschur == 0) return;

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Receive the set of blocks

	int i, iblk, niloc;

	int ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		niloc = pbl2nd[iblk+1]-pbl2nd[iblk];
		ni += niloc;
	};

	dcmplx *darr;

	darr = new dcmplx [ni];
	if (!darr) MemoryFail (funcname);

	CMPIStatus status;

	CMPIExchange::Recv (ptree->GetComm (), _jproc, 1,
								ni*sizeof(dcmplx), (char *)darr, status);

//	OutArr (*pfout," Recv from ",1,&_jproc);
//	OutArr (*pfout," Recv list ",_nlistschur,_listschur);
//	OutArr (*pfout," Recv arr ",ni,darr);

// Store the data

	dcmplx *pdvect = pvect1->GetVect ();

	int ibs, j;

	ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		ibs = ibsvect[iblk];
		niloc = pbl2nd[iblk+1]-pbl2nd[iblk];
		if (_add) {
			for (j=0;j<niloc;j++) pdvect[ibs+j] += darr[ni+j];
		} else {
			for (j=0;j<niloc;j++) pdvect[ibs+j] = darr[ni+j];
		};
		ni += niloc;
	};

// Free work data

	delete [] darr;

};

// Author: Kharchenko S.A.
// CMvmSchurC: Wait for completion of sends
//========================================================================================
void CMvmSchurC::WaitSends () { // Wait for completion of sends

	const char *funcname = "WaitSends";

	int nloc = nsends;

	CMPIRequest *reqvarrloc;
	CMPIStatus *statarrloc;

	reqvarrloc = new CMPIRequest [nloc];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new CMPIStatus [nloc];
	if (!statarrloc) MemoryFail (funcname);

	nloc = 0;

	int i;

	for (i=0;i<nsends;i++) {
		if (imasksend[i] == 1) {
			reqvarrloc[nloc] = sendreqvarr[i];
			nloc++;
		};
	};

	CMPIExchange::WaitAll (nloc, reqvarrloc, statarrloc);

	for (i=0;i<nsends;i++) {
		if (psenddata[i] != 0) {
			delete [] psenddata[i];
		};
	};

	delete [] reqvarrloc;
	delete [] statarrloc;

	nsends = 0;

};

// Author: Kharchenko S.A.
// Description: Store the result
// CMvmSchurC::StorePx()
//========================================================================================
void CMvmSchurC::StorePx (CSVectorC &_px) { // Store the result

//	const char *funcname = "StorePx";

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Init X part

	dcmplx *px = _px.GetVect ();
	dcmplx *pdvect = pvect1->GetVect ();

	int i, iblk, ni, ibs, j;

	int nloc = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pbl2nd[iblk+1]-pbl2nd[iblk];
		for (j=0;j<ni;j++) px[nloc+j] = pdvect[ibs+j];
		nloc += ni;
	};

};

// Author: Kharchenko S.A.
// Description: Destroy solve working data
// CMvmSchurC::DestroyWorkData()
//========================================================================================
void CMvmSchurC::DestroyWorkData () { // Destroy solve working data

	const char *funcname = "DestroyWorkData";

	delete [] imaskflt;
	delete [] ibsvect;
	delete [] listblk;
	delete [] listblkext;
	delete [] iasndblk;
	delete [] jasndblk;
	delete [] iarcvblk;
	delete [] jarcvblk;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [0];
	if (!ibsvect) MemoryFail (funcname);
	listblk = new int [0];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [0];
	if (!listblkext) MemoryFail (funcname);
	iasndblk = new int [0];
	if (!iasndblk) MemoryFail (funcname);
	jasndblk = new int [0];
	if (!jasndblk) MemoryFail (funcname);
	iarcvblk = new int [0];
	if (!iarcvblk) MemoryFail (funcname);
	jarcvblk = new int [0];
	if (!jarcvblk) MemoryFail (funcname);

	nlistblk = 0;
	nlistblkext = 0;

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	delete pvect;
	delete pvect1;

	pvect = new CSVectorC;
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVectorC;
	if (!pvect1) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destroy solve working data
// CParSchurMvmSlvCS::CreateIbssup()
//========================================================================================
void CParSchurMvmSlvCS::CreateIbssup () { // Create ibssup data for local operations

	const char *funcname = "CreateIbssup";

	if (ibssup != 0) {
		delete [] ibssup;
	};

	int *pibsvect = pmvm->GetIbsvect ();
	int nblksloc = ptree->GetNblkstree ();
	int *pblksloc = ptree->GetBlkstree ();
	int nsuploc = pgmtral->GetNsupr ();
	int *sprndsloc = pgmtral->GetSprndr ();

	ibssup = new int [nsuploc];
	if (!ibssup) MemoryFail (funcname);

	int i, iblk, ibs, j;

	for (i=0;i<nsuploc;i++) ibssup[i] = -1;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (pibsvect[iblk]>= 0) {
			ibs = pibsvect[iblk];
			for (j=pblksloc[iblk];j<pblksloc[iblk+1];j++) {
				ibssup[j] = ibs;
				ibs += sprndsloc[j+1]-sprndsloc[j];
			};
		};
	};

};
