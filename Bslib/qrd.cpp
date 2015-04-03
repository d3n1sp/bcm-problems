//------------------------------------------------------------------------------------------------
// File: qrd.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>

#include "tree.h"
#include "globals.h"
#include "qrd.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CQrd: Constructor
//========================================================================================
CQrd::CQrd (const CTree &_tree, int _nupdmax, int _blksiz) { // Constructor

	const char *funcname = "CQrdC_00";

// Allocate the data

	nlev    = _tree.nlev;
	nupdmax = _nupdmax;
	blksiz  = _blksiz;
	nfiles = 0;

	nchildsmax = 0;
	int inode;
	for (inode=0;inode<_tree.nnodes;inode++) {
		int nchilds = _tree.nodes[inode].nchilds;
		if (nchilds > nchildsmax) nchildsmax = nchilds;
	};

	mlocarr = new int * [nlev];
	if (!mlocarr) MemoryFail (funcname);
	bsqarr = new int * [nlev];
	if (!bsqarr) MemoryFail (funcname);

	int ilev;
	for (ilev=0;ilev<nlev;ilev++) {
		int *mblk, *bsqloc;

		mblk = new int [_nupdmax];
		if (!mblk) MemoryFail (funcname);

		mlocarr[ilev] = mblk;

		bsqloc = new int [_nupdmax];
		if (!bsqloc) MemoryFail (funcname);

		bsqarr[ilev] = bsqloc;

	};

	qarr = new double * [nlev];
	if (!qarr) MemoryFail (funcname);
	tauarr = new double * [nlev];
	if (!tauarr) MemoryFail (funcname);
	qrecvarr = new double [nupdmax*blksiz*blksiz];
	if (!qrecvarr) MemoryFail (funcname);
	qsendarr = new double [nupdmax*blksiz*blksiz];
	if (!qsendarr) MemoryFail (funcname);
	bl2file = new int [nfiles];
	if (!bl2file) MemoryFail (funcname);
	bsblk = new int [nfiles];
	if (!bsblk) MemoryFail (funcname);
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);

	int i;
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

// Search the tree

	int nodeend = _tree.cpuidend[_tree.myid];

	inode = nodeend;

	while (inode >= 0) {

		ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

		int fathernode = _tree.nodes[inode].fatherid;

// Allocate and init all necessary arrays 

		int *mblk = mlocarr[ilev];
		int *bsqloc = bsqarr[ilev];

		int nzq = 0;

		if (_tree.nodes[inode].nodecpu == _tree.myid) {
			for (int iupd=0;iupd<nupdmax;iupd++) {
				int mloc = (iupd+1) * nchilds * blksiz;
				mblk[iupd] = mloc;
				bsqloc[iupd] = nzq;
				nzq += mloc*blksiz;
			};
		};

		double *qloc;
		double *tauloc;

		qloc = new double [nzq];
		if (!qloc) MemoryFail (funcname);
		tauloc = new double [_nupdmax*_blksiz];
		if (!tauloc) MemoryFail (funcname);

		qarr[ilev] = qloc;
		tauarr[ilev] = tauloc;

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			inode = fathernode;
		};

	};

};

// Author: Kharchenko S.A.
// Update QR decomposition of the complex block
//========================================================================================
void CQrd::UpdateQrdBlk (int _ilev, int _iblk, // Update QR decomposition of the block
							int _ldq, double *_qblk) {

//	const char *funcname = "UpdateQrdBlk";

	int i, j, k, kk;
	double scprod;

	double *pa, *ph;

// Copy initial data

	int *bsqloc = bsqarr[_ilev];
	int *mblkloc = mlocarr[_ilev];
	double *qloc = qarr[_ilev];
	double *tauloc = tauarr[_ilev];

	int ibs = bsqloc[_iblk];
	int mloc = mblkloc[_iblk];

	for (i=0;i<blksiz;i++) {
		for (j=0;j<mloc;j++) {
			qloc[ibs+i*mloc+j] = _qblk[i*_ldq+j];
		};
	};

// Perform updates of the current block by the previous ones

	for (int iblkrd=0;iblkrd<_iblk;iblkrd++) {

		int jbs = bsqloc[iblkrd];
		int mlocj = mblkloc[iblkrd];

		for (i=0;i<blksiz;i++) {

			for (j=0;j<blksiz;j++) {

				int jnew = j+iblkrd*blksiz;

				scprod = qloc[ibs+i*mloc+jnew];

				pa = qloc+ibs+i*mloc+jnew+1;
				ph = qloc+jbs+j*mlocj+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
//					scprod += _a[i*_m+k] * _a[j*_m+k];
					scprod += *pa++ * *ph++;
//					scprod.PlEq (*pa++ * conj(*ph++));
				};

				scprod *= tauloc[jnew];
//				scprod.MlEq (conj(tauloc[jnew]));

				qloc[ibs+i*mloc+jnew] -= scprod;

				pa = qloc+ibs+i*mloc+jnew+1;
				ph = qloc+jbs+j*mlocj+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
//					_a[i*_m+k] -= scprod*_a[j*_m+k];
//					*pa++ -= scprod * conj(*ph++);
					*pa++ -= scprod * *ph++;
//					(*pa++).MnEq (scprod * *ph++);
				};
			};

		};

	};

	for (i=0;i<blksiz;i++) {

		int inew = i+_iblk*blksiz;

// Apply previous columns stored in-core to the current column

		for (j=0;j<i;j++) {

			int jnew = j+_iblk*blksiz;

			scprod = qloc[ibs+i*mloc+jnew];

			pa = qloc+ibs+i*mloc+jnew+1;
			ph = qloc+ibs+j*mloc+jnew+1;

			for (k=jnew+1;k<mloc;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
				scprod += *pa++ * *ph++;
//				scprod.PlEq (*pa++ * conj(*ph++));
			};

			scprod *= tauloc[jnew];
//			scprod.MlEq (conj(tauloc[jnew]));

			qloc[ibs+i*mloc+jnew] -= scprod;

			pa = qloc+ibs+i*mloc+jnew+1;
			ph = qloc+ibs+j*mloc+jnew+1;

			for (k=jnew+1;k<mloc;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
				*pa++ -= scprod * *ph++;
//				(*pa++).MnEq (scprod * *ph++);
			};
		};

// Compute new transformation

		j = mloc-1;
		if (inew+1 < mloc-1) j = inew+1;

		k = mloc-inew;
		kk = 1;

		dlarfg_ (&k, qloc+ibs+i*mloc+inew, qloc+ibs+i*mloc+j, &kk, tauloc+inew);

	};

};

// Author: Kharchenko S.A.
// Multiply Q factor by the current block
//========================================================================================
void CQrd::MvmQBlk (int _ilev, int _iblk, int _nrhs, // Multiply Q factor by the current block
						double *_x, int _ldx , double *_qx, int _ldqx) { 

//	const char *funcname = "MvmQBlk";

	double dzero = 0.0e0;

	int irhs, j, k;
	double scprod;

	double *pq, *pqx;

	int *bsqloc = bsqarr[_ilev];
	int *mblkloc = mlocarr[_ilev];
	double *qloc = qarr[_ilev];
	double *tauloc = tauarr[_ilev];

	int mloc = mblkloc[_iblk];

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<mloc;j++) *pqx++ = dzero;
		for (j=0;j<(_iblk+1)*blksiz;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};
	};

// Main cycle over columns

	int iblkrd;

	for (iblkrd=_iblk;iblkrd>=0;iblkrd--) {

		int jbs = bsqloc[iblkrd];
		int mlocj = mblkloc[iblkrd];

// Apply current part of Q to the set of columns

		for (j=blksiz-1;j>=0;j--) {

			int jnew = j+iblkrd*blksiz;

			for (irhs=0;irhs<_nrhs;irhs++) {

				scprod = _qx[irhs*_ldqx+jnew];

				pq = qloc+jbs+j*mlocj+jnew+1;
				pqx = _qx+irhs*_ldqx+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
					scprod += *pq++ * *pqx++;
//					scprod.PlEq (*pqx++ * conj(*pq++));
				};

				scprod = -scprod * tauloc[jnew];

				_qx[irhs*_ldqx+jnew] += scprod;

				pq = qloc+jbs+j*mlocj+jnew+1;
				pqx = _qx+irhs*_ldqx+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
					*pqx++ += scprod * *pq++;
//					(*pqx++).PlEq (scprod * *pq++);
				};

			};

		};

	};

};

// Author: Kharchenko S.A.
// Get R part of the QR decomposition
//========================================================================================
void CQrd::GetRPartQrd (const CTree &_tree, int _ilev, int _ibeg, int _iend, // Get R part of the QR decomposition
							double *_r, int _ldr) {

	double dzero = 0.0e0;

	int kii, kjj;

	for (kii=0;kii<_iend-_ibeg;kii++) {
		for (kjj=0;kjj<_iend;kjj++) {
			_r[kii*_ldr+kjj] = dzero;
		};
	};

	if (_ilev != 0 || _tree.myid == _tree.rootcpu) {
		int *mblk = mlocarr[_ilev];
		int *bsqloc = bsqarr[_ilev];
		double *qloc;
		qloc = qarr[_ilev];
		if (_ibeg%blksiz != 0) throw " Incorrect ibeg parameter in GetRPartQrd routine";
		if (_iend%blksiz != 0) throw " Incorrect iend parameter in GetRPartQrd routine";
		int iupdbeg = _ibeg / blksiz;
		int iupdend = _iend / blksiz;
		for (int iupd = iupdbeg;iupd<iupdend;iupd++) {
			int mloc = mblk[iupd];
			int ibs = bsqloc[iupd];
			for (kii=iupd*blksiz;kii<(iupd+1)*blksiz;kii++) {
				for (kjj=0;kjj<=kii;kjj++) {
					_r[(kii-_ibeg)*_ldr+kjj] = qloc[ibs+(kii-iupd*blksiz)*mloc+kjj];
				};
			};
		};

	};

};

// Author: Kharchenko S.A.
// Qrd: Setup files
//========================================================================================
void CQrd::SetupFiles (int _nfiles, char *_name, int _ldq) { // Setup files

	const char *funcname = "SetupFiles";

	int i;

	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	delete [] bl2file;
	delete [] bsblk;
	delete [] files;
	delete [] fnames;

	nfiles = _nfiles;

	bl2file = new int [nupdmax+1+nfiles];
	if (!bl2file) MemoryFail (funcname);
	bsblk = new int [nupdmax+1+nfiles];
	if (!bsblk) MemoryFail (funcname);
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

	OpenDIOFiles (_nfiles, files, fnames, "", _name, ".bin");

	int ni = (nupdmax+1+_nfiles-1)/_nfiles;

	if (ni<1) ni = 1;

	int nz;

	int ip=0;

	for (int iblk=0;iblk<_nfiles-1;iblk++) {
		nz = 0;
		for (int j=0;j<ni;j++) {
			bl2file[ip] = iblk;
			bsblk[ip] = j*_ldq*blksiz;
			ip++;
		};
	};
	nz = 0;
	while (ip<nupdmax+1) {
		bl2file[ip] = _nfiles-1;
		bsblk[ip] = nz;
		ip++;
		nz += _ldq*blksiz;
	};

};

// Author: Kharchenko S.A.
// Qrd: Close files
//========================================================================================
void CQrd::CloseFiles () { // Close files

//	const char *funcname = "CloseFiles";

//	CloseDIOFiles (nfiles, files);
	CloseDIOFiles (nfiles, files, fnames);

};

// Author: Kharchenko S.A.
// CQrdC: Constructor
//========================================================================================
CQrdC::CQrdC (const CTree &_tree, int _nupdmax, int _blksiz) { // Constructor

	const char *funcname = "CQrdC_00";

// Allocate the data

	nlev    = _tree.nlev;
	nupdmax = _nupdmax;
	blksiz  = _blksiz;
	nfiles = 0;

	nchildsmax = 0;
	int inode;
	for (inode=0;inode<_tree.nnodes;inode++) {
		int nchilds = _tree.nodes[inode].nchilds;
		if (nchilds > nchildsmax) nchildsmax = nchilds;
	};

	mlocarr = new int * [nlev];
	if (!mlocarr) MemoryFail (funcname);
	bsqarr = new int * [nlev];
	if (!bsqarr) MemoryFail (funcname);

	int ilev;
	for (ilev=0;ilev<nlev;ilev++) {
		int *mblk, *bsqloc;

		mblk = new int [_nupdmax];
		if (!mblk) MemoryFail (funcname);

		mlocarr[ilev] = mblk;

		bsqloc = new int [_nupdmax];
		if (!bsqloc) MemoryFail (funcname);

		bsqarr[ilev] = bsqloc;

	};

	qarr = new dcmplx * [nlev];
	if (!qarr) MemoryFail (funcname);
	tauarr = new dcmplx * [nlev];
	if (!tauarr) MemoryFail (funcname);
	qrecvarr = new dcmplx [nupdmax*blksiz*blksiz];
	if (!qrecvarr) MemoryFail (funcname);
	qsendarr = new dcmplx [nupdmax*blksiz*blksiz];
	if (!qsendarr) MemoryFail (funcname);
	bl2file = new int [nfiles];
	if (!bl2file) MemoryFail (funcname);
	bsblk = new int [nfiles];
	if (!bsblk) MemoryFail (funcname);
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);

	int i;
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

// Search the tree

	int nodeend = _tree.cpuidend[_tree.myid];

	inode = nodeend;

	while (inode >= 0) {

		ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

		int fathernode = _tree.nodes[inode].fatherid;

// Allocate and init all necessary arrays 

		int *mblk = mlocarr[ilev];
		int *bsqloc = bsqarr[ilev];

		int nzq = 0;

		if (_tree.nodes[inode].nodecpu == _tree.myid) {
			for (int iupd=0;iupd<nupdmax;iupd++) {
				int mloc = (iupd+1) * nchilds * blksiz;
				mblk[iupd] = mloc;
				bsqloc[iupd] = nzq;
				nzq += mloc*blksiz;
			};
		};

		dcmplx *qloc;
		dcmplx *tauloc;

		qloc = new dcmplx [nzq];
		if (!qloc) MemoryFail (funcname);
		tauloc = new dcmplx [_nupdmax*_blksiz];
		if (!tauloc) MemoryFail (funcname);

		qarr[ilev] = qloc;
		tauarr[ilev] = tauloc;

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			inode = fathernode;
		};

	};

};

// Author: Kharchenko S.A.
// Update QR decomposition of the complex block
//========================================================================================
void CQrdC::UpdateQrdBlk (int _ilev, int _iblk, // Update QR decomposition of the block
							int _ldq, dcmplx *_qblk) {

//	const char *funcname = "UpdateQrdBlk";

	int i, j, k, kk;
	dcmplx scprod;

	dcmplx *pa, *ph;

// Copy initial data

	int *bsqloc = bsqarr[_ilev];
	int *mblkloc = mlocarr[_ilev];
	dcmplx *qloc = qarr[_ilev];
	dcmplx *tauloc = tauarr[_ilev];

	int ibs = bsqloc[_iblk];
	int mloc = mblkloc[_iblk];

	for (i=0;i<blksiz;i++) {
		for (j=0;j<mloc;j++) {
			qloc[ibs+i*mloc+j] = _qblk[i*_ldq+j];
		};
	};

// Perform updates of the current block by the previous ones

	for (int iblkrd=0;iblkrd<_iblk;iblkrd++) {

		int jbs = bsqloc[iblkrd];
		int mlocj = mblkloc[iblkrd];

		for (i=0;i<blksiz;i++) {

			for (j=0;j<blksiz;j++) {

				int jnew = j+iblkrd*blksiz;

				scprod = qloc[ibs+i*mloc+jnew];

				pa = qloc+ibs+i*mloc+jnew+1;
				ph = qloc+jbs+j*mlocj+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
//					scprod += _a[i*_m+k] * _a[j*_m+k];
//					scprod += *pa++ * conj(*ph++);
					scprod.PlEq (*pa++ * conj(*ph++));
				};

//				scprod *= conj(tauloc[jnew]);
				scprod.MlEq (conj(tauloc[jnew]));

				qloc[ibs+i*mloc+jnew] -= scprod;

				pa = qloc+ibs+i*mloc+jnew+1;
				ph = qloc+jbs+j*mlocj+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
//					_a[i*_m+k] -= scprod*_a[j*_m+k];
//					*pa++ -= scprod * conj(*ph++);
//					*pa++ -= scprod * *ph++;
					(*pa++).MnEq (scprod * *ph++);
				};
			};

		};

	};

	for (i=0;i<blksiz;i++) {

		int inew = i+_iblk*blksiz;

// Apply previous columns stored in-core to the current column

		for (j=0;j<i;j++) {

			int jnew = j+_iblk*blksiz;

			scprod = qloc[ibs+i*mloc+jnew];

			pa = qloc+ibs+i*mloc+jnew+1;
			ph = qloc+ibs+j*mloc+jnew+1;

			for (k=jnew+1;k<mloc;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
//				scprod += *pa++ * conj(*ph++);
				scprod.PlEq (*pa++ * conj(*ph++));
			};

//			scprod *= conj(tauloc[jnew]);
			scprod.MlEq (conj(tauloc[jnew]));

			qloc[ibs+i*mloc+jnew] -= scprod;

			pa = qloc+ibs+i*mloc+jnew+1;
			ph = qloc+ibs+j*mloc+jnew+1;

			for (k=jnew+1;k<mloc;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
//				*pa++ -= scprod * *ph++;
				(*pa++).MnEq (scprod * *ph++);
			};
		};

// Compute new transformation

		j = mloc-1;
		if (inew+1 < mloc-1) j = inew+1;

		k = mloc-inew;
		kk = 1;

		zlarfg_ (&k, qloc+ibs+i*mloc+inew, qloc+ibs+i*mloc+j, &kk, tauloc+inew);

	};

};

// Author: Kharchenko S.A.
// Multiply Q factor by the current block
//========================================================================================
void CQrdC::MvmQBlk (int _ilev, int _iblk, int _nrhs, // Multiply Q factor by the current block
						dcmplx *_x, int _ldx , dcmplx *_qx, int _ldqx) { 

//	const char *funcname = "MvmQBlk";

	dcmplx czero (0.0e0,0.0e0);

	int irhs, j, k;
	dcmplx scprod;

	dcmplx *pq, *pqx;

	int *bsqloc = bsqarr[_ilev];
	int *mblkloc = mlocarr[_ilev];
	dcmplx *qloc = qarr[_ilev];
	dcmplx *tauloc = tauarr[_ilev];

	int mloc = mblkloc[_iblk];

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<mloc;j++) *pqx++ = czero;
		for (j=0;j<(_iblk+1)*blksiz;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};
	};

// Main cycle over columns

	int iblkrd;

	for (iblkrd=_iblk;iblkrd>=0;iblkrd--) {

		int jbs = bsqloc[iblkrd];
		int mlocj = mblkloc[iblkrd];

// Apply current part of Q to the set of columns

		for (j=blksiz-1;j>=0;j--) {

			int jnew = j+iblkrd*blksiz;

			for (irhs=0;irhs<_nrhs;irhs++) {

				scprod = _qx[irhs*_ldqx+jnew];

				pq = qloc+jbs+j*mlocj+jnew+1;
				pqx = _qx+irhs*_ldqx+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
//					scprod += conj(*pq++) * *pqx++;
					scprod.PlEq (*pqx++ * conj(*pq++));
				};

				scprod = -scprod * tauloc[jnew];

				_qx[irhs*_ldqx+jnew] += scprod;

				pq = qloc+jbs+j*mlocj+jnew+1;
				pqx = _qx+irhs*_ldqx+jnew+1;

				for (k=jnew+1;k<mlocj;k++) {
//					*pqx++ += scprod * *pq++;
					(*pqx++).PlEq (scprod * *pq++);
				};

			};

		};

	};

};

// Author: Kharchenko S.A.
// Get R part of the QR decomposition
//========================================================================================
void CQrdC::GetRPartQrd (const CTree &_tree, int _ilev, int _ibeg, int _iend, // Get R part of the QR decomposition
							dcmplx *_r, int _ldr) {

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj;

	for (kii=0;kii<_iend-_ibeg;kii++) {
		for (kjj=0;kjj<_iend;kjj++) {
			_r[kii*_ldr+kjj] = czero;
		};
	};

	if (_ilev != 0 || _tree.myid == _tree.rootcpu) {
		int *mblk = mlocarr[_ilev];
		int *bsqloc = bsqarr[_ilev];
		dcmplx *qloc;
		qloc = qarr[_ilev];
		if (_ibeg%blksiz != 0) throw " Incorrect ibeg parameter in GetRPartQrd routine";
		if (_iend%blksiz != 0) throw " Incorrect iend parameter in GetRPartQrd routine";
		int iupdbeg = _ibeg / blksiz;
		int iupdend = _iend / blksiz;
		for (int iupd = iupdbeg;iupd<iupdend;iupd++) {
			int mloc = mblk[iupd];
			int ibs = bsqloc[iupd];
			for (kii=iupd*blksiz;kii<(iupd+1)*blksiz;kii++) {
				for (kjj=0;kjj<=kii;kjj++) {
					_r[(kii-_ibeg)*_ldr+kjj] = qloc[ibs+(kii-iupd*blksiz)*mloc+kjj];
				};
			};
		};

	};

};

// Author: Kharchenko S.A.
// QrdC: Setup files
//========================================================================================
void CQrdC::SetupFiles (int _nfiles, char *_name, int _ldq) { // Setup files

	const char *funcname = "SetupFiles";

	int i;

	for (i=0;i<nfiles;i++) {
		delete [] fnames[i];
	};
	delete [] bl2file;
	delete [] bsblk;
	delete [] files;
	delete [] fnames;

	nfiles = _nfiles;

	bl2file = new int [nupdmax+1+nfiles];
	if (!bl2file) MemoryFail (funcname);
	bsblk = new int [nupdmax+1+nfiles];
	if (!bsblk) MemoryFail (funcname);
	files = new FILE * [nfiles];
	if (!files) MemoryFail (funcname);
	fnames = new char * [nfiles];
	if (!fnames) MemoryFail (funcname);
	for (i=0;i<nfiles;i++) {
		fnames[i] = new char [512];
		if (!fnames[i]) MemoryFail (funcname);
	};

	OpenDIOFiles (_nfiles, files, fnames, "", _name, ".bin");

	int ni = (nupdmax+1+_nfiles-1)/_nfiles;

	if (ni<1) ni = 1;

	int nz;

	int ip=0;

	for (int iblk=0;iblk<_nfiles-1;iblk++) {
		nz = 0;
		for (int j=0;j<ni;j++) {
			bl2file[ip] = iblk;
			bsblk[ip] = j*_ldq*blksiz;
			ip++;
		};
	};
	nz = 0;
	while (ip<nupdmax+1) {
		bl2file[ip] = _nfiles-1;
		bsblk[ip] = nz;
		ip++;
		nz += _ldq*blksiz;
	};

};

// Author: Kharchenko S.A.
// QrdC: Close files
//========================================================================================
void CQrdC::CloseFiles () { // Close files

//	const char *funcname = "CloseFiles";

//	CloseDIOFiles (nfiles, files);
	CloseDIOFiles (nfiles, files, fnames);

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
