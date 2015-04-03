//------------------------------------------------------------------------------------------------
// File: ilu2r.cpp
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
#include <iomanip>
#include <ctime>

#include "gsmatrix.h"
#include "smatrix.h"
#include "slvparam.h"
#include "globals.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSMatrixR: Allocate double data and init scaling for Ilu2
//========================================================================================
void CSMatrixR::Ilu2AllocScale (double _sclmin, // Allocate double data and init scaling for Ilu2
								double *&_fmaskl, double *&_fmasku, double *&_dpiv,
								double *&_diagg, double *&_scll, double *&_sclu) const {

	const char *funcname = "Ilu2AllocScale";

	_fmaskl = new double [n];
	if (!_fmaskl) MemoryFail (funcname);
	_fmasku = new double [n];
	if (!_fmasku) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_diagg = new double [n];
	if (!_diagg) MemoryFail (funcname);
	_scll = new double [n];
	if (!_scll) MemoryFail (funcname);
	_sclu = new double [n];
	if (!_sclu) MemoryFail (funcname);

// Compute scaling of the matrix

	int i;
	for (i=0;i<n;i++) {
		double aux = a[ia[i]];
		double auxloc = aux;
		if (auxloc < 0.0e0) auxloc = -auxloc;
		if (auxloc < _sclmin) auxloc = _sclmin;
		auxloc = sqrt(auxloc);
		_dpiv[i] = auxloc;
		auxloc = 1.0e0 / auxloc;
		_scll[i] = auxloc;
		if (aux < 0.0e0) {
			_sclu[i] = auxloc;
		} else {
			_sclu[i] = -auxloc;
		};
	};
//	for (i=0;i<n;i++) scll[i] = 1.0e0;
//	for (i=0;i<n;i++) sclu[i] = 1.0e0;

// Init working arrays

	for (i=0;i<n;i++) {
		_fmaskl[i] = 0.0e0;
		_fmasku[i] = 0.0e0;
		_diagg[i] = a[ia[i]]*_scll[i]*_sclu[i];
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate G data for Ilu2
//========================================================================================
void CSMatrixR::Ilu2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ilu2
							int *&_ig, int *&_jg, int *&_addrg, double *&_gl, double *&_gu) const {

	const char *funcname = "Ilu2AllocG";

	int nz = ia[nsupr];

	_nzjgmx = int (_memory * double(nz));
	_nzgmax = _nzjgmx;

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int [_nzjgmx];
	if (!_jg) MemoryFail (funcname);
	_addrg = new int [1];
	if (!_addrg) MemoryFail (funcname);
	_gl = new double [_nzgmax];
	if (!_gl) MemoryFail (funcname);
	_gu = new double [_nzgmax];
	if (!_gu) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate G data for Ilu2
//========================================================================================
void CSMatrixR::Ilu2AllocGDynamic ( // Allocate G data for Ilu2
									int *&_ig, int **&_jg, double **&_gl, double **&_gu) const {

	const char *funcname = "Ilu2AllocG";

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int * [nsupr+1];
	if (!_jg) MemoryFail (funcname);
	_gl = new double * [nsupr+1];
	if (!_gl) MemoryFail (funcname);
	_gu = new double * [nsupr+1];
	if (!_gu) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate G data for Ilu2
//========================================================================================
void CSMatrixR::Ilu2AllocGDynamicByBlocks ( // Allocate G data for Ilu2
											int *&_ig, int **&_jg, double **&_gl, double **&_gu,
											int *&_blks, int *&_nd2blk, int **&_jgblk,
											double **&_glblk, double **&_gublk) const {

	const char *funcname = "Ilu2AllocGDynamicByBlocks";

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int * [nsupr+1];
	if (!_jg) MemoryFail (funcname);
	_gl = new double * [nsupr+1];
	if (!_gl) MemoryFail (funcname);
	_gu = new double * [nsupr+1];
	if (!_gu) MemoryFail (funcname);
	_blks = new int [nsupr+1];
	if (!_blks) MemoryFail (funcname);
	_nd2blk = new int [nsupr];
	if (!_nd2blk) MemoryFail (funcname);
	_jgblk = new int * [nsupr+1];
	if (!_jgblk) MemoryFail (funcname);
	_glblk = new double * [nsupr+1];
	if (!_glblk) MemoryFail (funcname);
	_gublk = new double * [nsupr+1];
	if (!_gublk) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixR: Init current scaled row
//========================================================================================
void CSMatrixR::Ilu2InitRow (const CSMatrixR &_mtral, const CSMatrixR &_mtrau, // Init current scaled row
								int _irow, int &_nlist1,
								int *_lstloc, int _icycle, int *_imask, 
								double *_diagg, double *_fmaskl, double *_fmasku, 
								double *_scll, double *_sclu, int &_iops) const {

//	const char *funcname = "Ilu2InitRow";

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	_fmaskl[_irow] = _diagg[_irow];
	_fmasku[_irow] = _diagg[_irow];
	_nlist1++;

	for (int j=_mtral.ia[_irow]+1;j<_mtral.ia[_irow+1];j++) {
		int jj = _mtral.ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;
		_fmasku[jj] = _mtrau.a[j]*_scll[_irow]*_sclu[jj];
		_fmaskl[jj] = _mtral.a[j]*_sclu[_irow]*_scll[jj];
		_nlist1++;
	};

	_iops += _nlist1*2;

};

// Author: Kharchenko S.A.
// CSMatrixR: Init mask and update arrays
//========================================================================================
void CSMatrixR::Ilu2InitMask (bool _elmisr, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const {

//	const char *funcname = "Ilu2InitMask";

	_nupd = 0;

	for (int k=_iv[_irwprv];k<_ig[_irwprv+1];k++) {

		int jcolmn = _jg[k];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (_imask[jcolmn] != _icycle) {
			_imask[jcolmn] = _icycle;
			_lstloc[_nlist1] = jcolmn;
			_fmaskl[jcolmn] = 0.0e0;
			_fmasku[jcolmn] = 0.0e0;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = k;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Init mask and update arrays
//========================================================================================
void CSMatrixR::Ilu2InitMaskDynamic (bool _elmisr, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int **_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const {

//	const char *funcname = "Ilu2InitMaskDynamic";

	_nupd = 0;

	int *jgloc = _jg[_irwprv];

	for (int k=_iv[_irwprv];k<_ig[_irwprv+1];k++) {

		int kloc = k-_ig[_irwprv];

		int jcolmn = jgloc[kloc];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (_imask[jcolmn] != _icycle) {
			_imask[jcolmn] = _icycle;
			_lstloc[_nlist1] = jcolmn;
			_fmaskl[jcolmn] = 0.0e0;
			_fmasku[jcolmn] = 0.0e0;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = kloc;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Perform filtering of the current row
//========================================================================================
void CSMatrixR::Ilu2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
							int &_nlist1, int *_lstloc, int *_lstupd,
							double *_fmaskl, double *_fmasku, double *_diagg, int &_iops) const {

//	const char *funcname = "Ilu2FiltrRow";


	int nlist2 = 0;

	for (int i=1;i<_nlist1;i++) {

		int j = _lstloc[i];
		double auxloc = _fmaskl[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		double aux = auxloc;
		auxloc = _fmasku[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		if (auxloc > aux) aux = auxloc;

		if (aux >= _tau2) {
			_lstupd[nlist2] = j;
			nlist2++;
		} else {
			aux *= _theta;
			_fmaskl[_irow] += aux;
			_fmasku[_irow] += aux;
			_diagg[j] += aux;
		};

	};

	_iops += 3*(_nlist1-nlist2);

	_nlist1 = nlist2;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute diagonal pivot and scale current row
//========================================================================================
void CSMatrixR::Ilu2PivScaleRow (double _pivmin, // Compute diagonal pivot and scale current row
								int _irow, double &_deleml, double &_delemu,
								int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, 
								double *_dpiv, int &_iops) const {

//	const char *funcname = "Ilu2PivScaleRow";

	double aux = _fmaskl[_irow];

	double auxloc = aux;
	if (auxloc < 0.0e0) auxloc = -auxloc;

	if (auxloc < _pivmin) auxloc = _pivmin;

	auxloc = sqrt(auxloc);
	_dpiv[_irow] = auxloc;

	auxloc = 1.0e0 / auxloc;

	_deleml = auxloc;

	if (aux < 0.0e0) {
		_delemu = -auxloc;
	} else {
		_delemu = auxloc;
	};

	_iops += 2;

// Update the data in the list

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		_fmaskl[j] *= _delemu;
		_fmasku[j] *= _deleml;
	};

	_iops += 2*_nlist1;

};

// Author: Kharchenko S.A.
// CSMatrix: Update arrays that support transposed structures search
//========================================================================================
void CSMatrix::Ich2UpdateTranspDynamic (int _irow, int *_ibegm, int *_madj, // Update arrays that support transposed structures search
										int *_iv, int *_ig, int **_jg) const {

//	const char *funcname = "Ich2UpdateTransp";

	int irwprv = _ibegm[_irow];

	while (irwprv != -1) {

		int irwpr1 = _madj[irwprv];

		if (_iv[irwprv] >= _ig[irwprv+1]-1) {

			_madj[irwprv] = -1;

		} else {

			int j = _iv[irwprv]+1;
			int *jgloc = _jg[irwprv];
			int jloc = j-_ig[irwprv];
			int icolmn = jgloc[jloc];
			if (icolmn < 0) icolmn = -icolmn;
			_madj[irwprv] = _ibegm[icolmn];
			_ibegm[icolmn] = irwprv;

		};

		_iv[irwprv]++;

		irwprv = irwpr1;

	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Free R part of the data
//========================================================================================
void CSMatrixR::Ilu2FreeR (int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, double *_gl, double *_gu) const {

//	const char *funcname = "Ilu2FreeR";

	int iendfr = _irow;

	cout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;

	int ibgfrn = _ibegfr;

	_nzjg = _ig[_ibegfr];

	int ibeg = _ig[_ibegfr];

	int nzjro = _nzjrl;
	_nzjrl = 0;

	for (int i=_ibegfr;i<iendfr;i++) {

// Scan current row

		for (int j=ibeg;j<_ig[i+1];j++) {

			int icolmn = _jg[j];

			if (icolmn < 0) {
				icolmn = -icolmn;
				if (icolmn < iendfr) goto label2;
				_nzjrl++;
			};
	
			if (_nzjg != j) {
				_jg[_nzjg] = _jg[j];
				_gl[_nzjg] = _gl[j];
				_gu[_nzjg] = _gu[j];
			};
			_nzjg++;

label2:;

		};

		ibeg = _ig[i+1];

		_ig[i+1] = _nzjg;

		if (_nzjrl == 0 && ibgfrn == i) ibgfrn = i+1;

// Update iv

		_iv[i] += _nzjg-ibeg;

	};

	_ibegfr = ibgfrn;

	cout << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;

};

// Author: Kharchenko S.A.
// CSMatrixR: Free R part of the data
//========================================================================================
void CSMatrixR::Ilu2FreeRDynamic (int _msglev, int _irow, // Free R part of the data
								int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
								int *_iv, int *_ig, int **_jg, double **_gl, double **_gu,
								int *_list, double *_fmaskl, double *_fmasku) const {

	const char *funcname = "Ilu2FreeRDynamic";

	int iendfr = _irow;

	if (_msglev > 1) {
		cout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;
	};

	int ibgfrn = _ibegfr;

	_nzjg = _ig[_ibegfr];

	int ibeg = _ig[_ibegfr];

	int nzjro = _nzjrl;
	_nzjrl = 0;

	int i, j, jloc, nlistloc;
	int *jgloc;
	double *glloc;
	double *guloc;

	for (i=_ibegfr;i<iendfr;i++) {

// Initial scan of the current row

		jgloc = _jg[i];
		glloc = _gl[i];
		guloc = _gu[i];

		nlistloc = 0;

		for (j=ibeg;j<_ig[i+1];j++) {

			jloc = j-ibeg;

			int icolmn = jgloc[jloc];

			if (icolmn < 0) {
				icolmn = -icolmn;
				if (icolmn < iendfr) goto label;
				_nzjrl++;
			};


			_list  [nlistloc] = jgloc[jloc];
			_fmaskl[nlistloc] = glloc[jloc];
			_fmasku[nlistloc] = guloc[jloc];
			nlistloc++;
			_nzjg++;

label:;

		};

// Free/allocate the memory

		delete [] _jg[i];
		delete [] _gl[i];
		delete [] _gu[i];

		jgloc = new int [nlistloc];
		if (!jgloc) MemoryFail (funcname);
		glloc = new double [nlistloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [nlistloc];
		if (!guloc) MemoryFail (funcname);

// Final scan of the current row

		for (j=0;j<nlistloc;j++) {
			jgloc[j] = _list[j];
			glloc[j] = _fmaskl[j];
			guloc[j] = _fmasku[j];
		};

		_jg[i] = jgloc;
		_gl[i] = glloc;
		_gu[i] = guloc;

		ibeg = _ig[i+1];

		_ig[i+1] = _nzjg;

		if (_nzjrl == 0 && ibgfrn == i) ibgfrn = i+1;

// Update iv

		_iv[i] += _nzjg-ibeg;

	};

	_ibegfr = ibgfrn;

	if (_msglev > 1) {
		cout << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Free R part of the data
//========================================================================================
void CSMatrixR::Ilu2FreeRDynamicByBlocks (int _msglev, int _irow, // Free R part of the data
										int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
										int *_iv, int *_ig, int **_jg, double **_gl, double **_gu,
										int _iblk, int _nzalloc, int &_nzused,
										int *_blks, int *_nd2blk, int **_jgblk, double **_glblk, double **_gublk,
										int *_list, double *_fmaskl, double *_fmasku) const {

	const char *funcname = "Ilu2FreeRDynamicByBlocks";

	int iendfr = _irow;
	int ibegblk = _nd2blk[_ibegfr];

	if (_msglev > 1) {
		cout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << " Ibegblk = " << ibegblk << " Iendblk = " << _iblk << endl;
	};

	_ibegfr = _blks[ibegblk];

	int ibgfrn = _ibegfr;

	_nzjg = _ig[_ibegfr];

	int ibeg = _ig[_ibegfr];

	int nzjro = _nzjrl;
	_nzjrl = 0;

	int i, j, jloc, nlistloc, iblkloc, ibegnd, iendnd, nzblk, ibegloc;
	int *jgloc;
	double *glloc;
	double *guloc;
	int *jgarr;
	double *glarr;
	double *guarr;

	for (iblkloc=ibegblk;iblkloc<_iblk;iblkloc++) {

		ibegnd = _blks[iblkloc];
		iendnd = _blks[iblkloc+1];

// Initial scan of the current block

		nzblk = 0;

		for (i=ibegnd;i<iendnd;i++) {
			jgloc = _jg[i];
			if (i==ibegnd) {
				ibegloc = ibeg;
			} else {
				ibegloc = _ig[i];
			};
			for (j=ibegloc;j<_ig[i+1];j++) {
				jloc = j-ibegloc;
				int icolmn = jgloc[jloc];
				if (icolmn < 0) {
					icolmn = -icolmn;
					if (icolmn < iendfr) goto label;
				};
				nzblk++;
label:;
			};

		};

// Allocate new memory and store data of the block row by row

		jgarr = new int [nzblk];
		if (!jgarr) MemoryFail (funcname);
		glarr = new double [nzblk];
		if (!glarr) MemoryFail (funcname);
		guarr = new double [nzblk];
		if (!guarr) MemoryFail (funcname);

		nzblk = 0;

		for (i=ibegnd;i<iendnd;i++) {

// Initial scan of the current row

			jgloc = _jg[i];
			glloc = _gl[i];
			guloc = _gu[i];

			nlistloc = 0;

			for (j=ibeg;j<_ig[i+1];j++) {

				jloc = j-ibeg;

				int icolmn = jgloc[jloc];

				if (icolmn < 0) {
					icolmn = -icolmn;
					if (icolmn < iendfr) goto label1;
					_nzjrl++;
				};

				_list  [nlistloc] = jgloc[jloc];
				_fmaskl[nlistloc] = glloc[jloc];
				_fmasku[nlistloc] = guloc[jloc];
				nlistloc++;
				_nzjg++;

label1:;

			};

// Final scan of the current row

			jgloc = jgarr+nzblk;
			glloc = glarr+nzblk;
			guloc = guarr+nzblk;

			for (j=0;j<nlistloc;j++) {
				jgloc[j] = _list[j];
				glloc[j] = _fmaskl[j];
				guloc[j] = _fmasku[j];
			};

			_jg[i] = jgloc;
			_gl[i] = glloc;
			_gu[i] = guloc;

			ibeg = _ig[i+1];

			_ig[i+1] = _nzjg;

			if (_nzjrl == 0 && ibgfrn == i) ibgfrn = i+1;

// Update iv

			_iv[i] += _nzjg-ibeg;

			nzblk += nlistloc;

		};

// Free previous memory and store the new one

		delete [] _jgblk[iblkloc];
		delete [] _glblk[iblkloc];
		delete [] _gublk[iblkloc];

		_jgblk[iblkloc] = jgarr;
		_glblk[iblkloc] = glarr;
		_gublk[iblkloc] = guarr;

	};

// Scan the last block for completeness

	ibegnd = _blks[_iblk];
	iendnd = _blks[_iblk+1];

	for (i=ibegnd;i<iendnd;i++) {

// Initial scan of the current row

		jgloc = _jg[i];

		for (j=ibeg;j<_ig[i+1];j++) {

			jloc = j-ibeg;

			int icolmn = jgloc[jloc];

			if (icolmn < 0) {
				icolmn = -icolmn;
				_nzjrl++;
			};

			_nzjg++;

		};

		ibeg = _ig[i+1];

		_ig[i+1] = _nzjg;

// Update iv

		_iv[i] += _nzjg-ibeg;

	};

	_ibegfr = ibgfrn;

	if (_msglev > 1) {
		cout << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Store current row
//========================================================================================
void CSMatrixR::Ilu2StoreRow (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double _deleml, double _delemu, 
								int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku,
								int *_ig, int *_jg, double *_gl, double *_gu) const {

//	const char *funcname = "Ilu2StoreRow";

	_jg[_nzjg] = _irow;
	_gl[_nzjg] = _deleml;
	_gu[_nzjg] = _delemu;
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double auxloc = _fmaskl[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		double aux = auxloc;
		auxloc = _fmasku[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		if (auxloc > aux) aux = auxloc;

		if (aux >= _tau1) {
			_jg[_nzjg] = j;
			_gl[_nzjg] = _fmaskl[j];
			_gu[_nzjg] = _fmasku[j];
			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
			_gl[_nzjg] = _fmaskl[j];
			_gu[_nzjg] = _fmasku[j];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

};

// Author: Kharchenko S.A.
// CSMatrixR: Store current row
//========================================================================================
void CSMatrixR::Ilu2StoreRowDynamic (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double _deleml, double _delemu, 
								int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku,
								int *_ig, int **_jg, double **_gl, double **_gu) const {

	const char *funcname = "Ilu2StoreRowDynamic";

	int *jgloc = new int [_nlist1+1];
	if (!jgloc) MemoryFail (funcname);
	double *glloc = new double [_nlist1+1];
	if (!glloc) MemoryFail (funcname);
	double *guloc = new double [_nlist1+1];
	if (!guloc) MemoryFail (funcname);

	jgloc[0] = _irow;
	glloc[0] = _deleml;
	guloc[0] = _delemu;
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double auxloc = _fmaskl[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		double aux = auxloc;
		auxloc = _fmasku[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		if (auxloc > aux) aux = auxloc;

		if (aux >= _tau1) {
			jgloc[i+1] = j;
			glloc[i+1] = _fmaskl[j];
			guloc[i+1] = _fmasku[j];
			_nzjg++;
			_nzju++;
		} else {
			jgloc[i+1] = -j;
			glloc[i+1] = _fmaskl[j];
			guloc[i+1] = _fmasku[j];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

	_jg[_irow] = jgloc;
	_gl[_irow] = glloc;
	_gu[_irow] = guloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Store current row
//========================================================================================
void CSMatrixR::Ilu2StoreRowDynamicByBlocks (int _irow, double _tau1, // Store current row
												int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
												double _deleml, double _delemu, 
												int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku,
												int *_ig, int **_jg, double **_gl, double **_gu,
												int &_iblk, int _nzallocmax, int &_nzalloc, int &_nzused,
												int *_blks, int *_nd2blk, int **_jgblk, 
												double **_glblk, double **_gublk) const {

	const char *funcname = "Ilu2StoreRowDynamicByBlocks";

	int *jgloc;
	double *glloc, *guloc;

	if (_nlist1+1 > _nzalloc-_nzused) {
		_nzalloc = _nzallocmax;
		if (_nzalloc < _nlist1+1) _nzalloc = _nlist1+1;
		jgloc = new int [_nzalloc];
		if (!jgloc) MemoryFail (funcname);
		glloc = new double [_nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [_nzalloc];
		if (!guloc) MemoryFail (funcname);
		_iblk++;
		_jgblk[_iblk] = jgloc;
		_glblk[_iblk] = glloc;
		_gublk[_iblk] = guloc;
		_nzused = 0;
	} else {
		jgloc = _jgblk[_iblk]+_nzused;
		glloc = _glblk[_iblk]+_nzused;
		guloc = _gublk[_iblk]+_nzused;
	};

	jgloc[0] = _irow;
	glloc[0] = _deleml;
	guloc[0] = _delemu;
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double auxloc = _fmaskl[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		double aux = auxloc;
		auxloc = _fmasku[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		if (auxloc > aux) aux = auxloc;

		if (aux >= _tau1) {
			jgloc[i+1] = j;
			glloc[i+1] = _fmaskl[j];
			guloc[i+1] = _fmasku[j];
			_nzjg++;
			_nzju++;
		} else {
			jgloc[i+1] = -j;
			glloc[i+1] = _fmaskl[j];
			guloc[i+1] = _fmasku[j];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

	_jg[_irow] = jgloc;
	_gl[_irow] = glloc;
	_gu[_irow] = guloc;

	_blks[_iblk+1] = _irow+1;
	_nd2blk[_irow] = _iblk;
	_nzused += (_nlist1+1);

};

// Author: Kharchenko S.A.
// CSMatrixR: Store the resulting U and scale it
//========================================================================================
void CSMatrixR::Ich2StoreScaleUDynamic (double *_scl, // Store the resulting U and scale it
								int *_ig, int **_jg, double **_g,
								CSMatrixR &_temp, double &_ops) const {

//	const char *funcname = "Ich2StoreScaleUDynamic";

	_temp.m = n;
	_temp.n = n;
	_temp.nsupr = nsupr;
	_temp.nsupc = nsupc;
	_temp.nlist = nsupr;

	int i;

	for (i=0;i<nsupr;i++) _temp.list[i] = i;

	_temp.ia[0] = 0;

	int nz = 0;

	int jloc;
	int *jgloc;
	double *gloc;

	for (i=0;i<nsupr;i++) {
		jgloc = _jg[i];
		gloc = _g[i];
		for (int j=_ig[i];j<_ig[i+1];j++) {
			jloc = j-_ig[i];
			int jj = jgloc[jloc];
			if (jj >= 0) {
				_temp.ja[nz] = jj;
				if (jj == i) {
					_temp.a [nz] = gloc[jloc] * _scl[i];
				} else {
					_temp.a [nz] = gloc[jloc] / _scl[jj];
				};
				nz++;
			};
		};
		_temp.ia[i+1] = nz;
	};

	_temp.nzja = nz;
	_temp.nza = nz;

	_ops += nz;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Store the resulting U and scale it
//========================================================================================
void CGSMatrixR::Ich2StoreScaleUDynamic (int _nblks, int *_blks, double *_scl, // Store the resulting U and scale it
								int *_ig, int **_jgblk, double **_gblk,
								CGSMatrixR &_gich) {

	const char *funcname = "Ich2StoreScaleUDynamic";

	int nloc = _blks[_nblks];

	_gich.m = nloc;
	_gich.n = nloc;
	_gich.nblksr = _nblks;
	_gich.nblksc = _nblks;
	_gich.nsupr = nloc;
	_gich.nsupc = nloc;

	int iblk, isup, isupbeg, isupend;

	for (iblk=0;iblk<=_nblks;iblk++) _gich.blksr[iblk] = _blks[iblk];
	for (iblk=0;iblk<=_nblks;iblk++) _gich.blksc[iblk] = _blks[iblk];
	for (isup=0;isup<=nloc;isup++) _gich.sprndr[isup] = isup;
	for (isup=0;isup<=nloc;isup++) _gich.sprndc[isup] = isup;

	_gich.bl2ndr[0] = 0;
	for (iblk=0;iblk<_nblks;iblk++) {
		isupbeg=_gich.blksr[iblk];
		isupend=_gich.blksr[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		_gich.bl2ndr[iblk+1] = _gich.bl2ndr[iblk] + ni;
	};

	_gich.bl2ndc[0] = 0;
	for (iblk=0;iblk<_nblks;iblk++) {
		isupbeg=_gich.blksc[iblk];
		isupend=_gich.blksc[iblk+1]-1;
		int ni=isupend-isupbeg+1;
		_gich.bl2ndc[iblk+1] = _gich.bl2ndc[iblk] + ni;
	};

	_gich.nlist=0;

	for (iblk=0;iblk<_nblks;iblk++) {
		_gich.listb[_gich.nlist++] = iblk;
	};

// Create the set of CSMatrixR block rows

	CSMatrixR mtrdummy;

	int nsupl, ibeg;
	int nzjal, nzal, i, j, jj;

	for (iblk=0;iblk<_nblks;iblk++) {

		isupbeg = _blks[iblk];
		isupend = _blks[iblk+1]-1;

		nsupl = isupend-isupbeg+1;

// Allocate the data

		nzjal = 0;
		nzal = _ig[isupend+1]-_ig[isupbeg];

		for (i=0;i<nzal;i++) {
			jj = _jgblk[iblk][i];
			if (jj >= 0) nzjal++;
		};

		nzal = nzjal;

		CSMatrixR temp (nsupl, nzjal, nzal);

// Init the matrix

		temp.m = nloc;
		temp.n = nloc;
		temp.nsupr = nsupl;
		temp.nsupc = nsupl;
		temp.nlist = nsupl;
		temp.nzja = nzjal;
		temp.nza = nzal;
		temp.nzatot = nzal;

		nzjal = 0;

		temp.ia[0] = 0;

		ibeg = _ig[isupbeg];
		for (isup=isupbeg;isup<=isupend;isup++) {
			temp.list[isup-isupbeg] = isup;
			for (j=_ig[isup];j<_ig[isup+1];j++) {
				jj = _jgblk[iblk][j-ibeg];
				if (jj >= 0) {
					temp.ja[nzjal] = jj;
					if (jj == isup) {
						temp.a [nzjal] = _gblk[iblk][j-ibeg] * _scl[isup];
					} else {
						temp.a [nzjal] = _gblk[iblk][j-ibeg] / _scl[jj];
					};
					nzjal++;
				};
			};
			temp.ia[isup-isupbeg+1] = nzjal;
		};

// Store the matrix

		_gich.mtrarr[iblk] = temp;

		temp = mtrdummy;

// Free preconditioner data

		delete [] _jgblk[iblk];
		delete [] _gblk[iblk];

		_jgblk[iblk] = new int [0];
		if (!_jgblk[iblk]) MemoryFail (funcname);
		_gblk[iblk] = new double [0];
		if (!_gblk[iblk]) MemoryFail (funcname);

	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Second order Incomplete LU decomposition
//========================================================================================
void CSMatrixR::Ilu2 (ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition
						const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
						CSMatrixR &_mtrl, CSMatrixR &_mtru) { 

//	const char *funcname = "Ilu2";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzgmax, nzjgmx, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	double lelem, uelem, deleml, delemu;
	int *ig, *jg, *addrg;
	double *gl, *gu;
	bool elmisr;

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;
	int iops;

	icheck = _param.ichkfct;
	ops = 0.0e0;

// Init time measurement

	time0 = clock ();

// Define working arrays

	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	double *fmaskl, *fmasku, *dpiv, *diagg, *scll, *sclu;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	_mtral.Ilu2AllocScale (_param.sclmin, fmaskl, fmasku, dpiv, diagg, scll, sclu);

// Output hystogram of the scaling values

	int itype = 1;
	int nloc = _mtral.n;

	if (_param.msglev > 1) {
		csvhyst ("SclDia", nloc, dpiv);
	};

// Allocate factorization arrays

	_mtral.Ilu2AllocG (_param.memory, nzjgmx, nzgmax, ig, jg, addrg, gl, gu);

// Create dummy matrix

	CSMatrixR mtrdummy;

// Main cycle over the rows

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;

	ig[0] = 0;

	for (irow=0;irow<_mtral.n;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {
			if (_param.msglev > 1) {
				cout << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			};
		};

// Init current scaled row

		mtrdummy.Ilu2InitRow (_mtral, _mtrau,
							irow, nlist1, lstloc, icycle, imask, 
							diagg, fmaskl, fmasku, scll, sclu, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];
			lelem = gl[j];
			uelem = gu[j];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			_mtral.Ilu2InitMask (elmisr, irwprv, iv, ig, jg, 
									nupd, lstupd, indlst,
									nlist1, lstloc, icycle, imask, 
									fmaskl, fmasku);

// Perform updates

			_mtral.Ich2UpdateRow (nupd, lstupd, indlst, 
									uelem, fmaskl, gl, iops);
			_mtral.Ich2UpdateRow (nupd, lstupd, indlst, 
									lelem, fmasku, gu, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		_mtral.Ich2UpdateTransp (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		_mtral.Ilu2FiltrRow (irow, _param.theta, _param.tau2, 
								nlist1, lstloc, lstupd, 
								fmaskl, fmasku, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin, irow, deleml, delemu, 
								nlist1, lstupd, fmaskl, fmasku, dpiv, iops);

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx) {

// Free the memory containing R part that is not in use any more

			_mtral.Ilu2FreeR (irow, 
								ibegfr, nzjg, nzju, nzjrl, 
								iv, ig, jg, gl, gu);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx) {
			if (_param.msglev > 1) {
				cout << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
			};
			throw " Insufficient workspace in Ilu2";
		};

// Sort the list of indices

		if (nlist1 != 0) qsort (lstupd, nlist1, sizeof(int), compint);

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist1 != 0) {
			iv[irow] = nzjg+1;
			madj[irow] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = irow;
		};

// Store the data in the global arrays

		_mtral.Ilu2StoreRow (irow, _param.tau1,
							nzjg, nzju, nzjr, nzjrl,
							deleml, delemu, 
							nlist1, lstupd, fmaskl, fmasku,
							ig, jg, gl, gu);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

	itype = 2;
	nloc = _mtral.n;

	if (_param.msglev > 1) {
		csvhyst ("Ilu2Dia", nloc, dpiv);
	};

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] fmaskl;
	delete [] fmasku;
	delete [] dpiv;
	delete [] diagg;

// Perform the final filtering and scaling and compute matrix on return

	CSMatrixR temp (_mtral.n,nzju);

	_mtral.Ich2StoreScaleU (scll, ig, jg, gl, 
							temp, ops);

	_mtrl = temp;

	_mtral.Ich2StoreScaleU (sclu, ig, jg, gu, 
							temp, ops);

	_mtru = temp;

	temp = mtrdummy;

// Free factorization arrays 

	delete [] scll;
	delete [] sclu;
	delete [] ig;
	delete [] jg;
	delete [] addrg;
	delete [] gl;
	delete [] gu;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ich2 factorization statistics

	nz = _mtral.ia[_mtral.n];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	_param.density  = (double) nzju / (double) nz;
	_param.density2 = (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = 2*_mtral.ia[_mtral.n]-_mtral.n;

	ops = ops / (double) nz;

	if (_param.msglev > 1) {
		cout << " Ilu2 preconditioner generation statistics: " << endl;
		cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout << "     Costs = " << ops << " MvmA flops. " << endl;
		cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};
	if (_param.msglev >= 1) {
		_fout << " Ilu2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Second order Incomplete LU decomposition with dynamic memory allocation
//========================================================================================
void CSMatrixR::Ilu2Dynamic (ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition with dynamic memory allocation
								const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
								CSMatrixR &_mtrl, CSMatrixR &_mtru) { 

//	const char *funcname = "Ilu2Dynamic";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	double lelem, uelem, deleml, delemu;
	int *ig, **jg;
	double **gl, **gu;
	bool elmisr;

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;
	int iops;

	icheck = _param.ichkfct;
	ops = 0.0e0;

// Init time measurement

	time0 = clock ();

// Define working arrays

	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	double *fmaskl, *fmasku, *dpiv, *diagg, *scll, *sclu;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	_mtral.Ilu2AllocScale (_param.sclmin, fmaskl, fmasku, dpiv, diagg, scll, sclu);

// Output hystogram of the scaling values

	int itype = 1;
	int nloc = _mtral.n;

	if (_param.msglev > 1) {
		csvhyst ("SclDia", nloc, dpiv);
	};

// Allocate factorization arrays

	_mtral.Ilu2AllocGDynamic (ig, jg, gl, gu);

// Create dummy matrix

	CSMatrixR mtrdummy;

// Main cycle over the rows

	int jloc;
	int *jgloc;
	double *glloc, *guloc;

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;

	ig[0] = 0;

// Set free R parameters

	int nfreer = (int) _param.nfreer;

	int nrowsfree;

	if (nfreer <= 0) {
		nrowsfree = 0;
	} else {
		nrowsfree = _mtral.n / nfreer;
	};

	for (irow=0;irow<_mtral.n;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {
			if (_param.msglev > 1) {
				cout << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << endl;
			};
		};

// Init current scaled row

		mtrdummy.Ilu2InitRow (_mtral, _mtrau,
							irow, nlist1, lstloc, icycle, imask, 
							diagg, fmaskl, fmasku, scll, sclu, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			jgloc = jg[irwprv];
			glloc = gl[irwprv];
			guloc = gu[irwprv];

			j = iv[irwprv];

			jloc = j-ig[irwprv];

			icolmn = jgloc[jloc];
			lelem = glloc[jloc];
			uelem = guloc[jloc];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			_mtral.Ilu2InitMaskDynamic (elmisr, irwprv, iv, ig, jg,
										nupd, lstupd, indlst,
										nlist1, lstloc, icycle, imask, 
										fmaskl, fmasku);

// Perform updates

			_mtral.Ich2UpdateRow (nupd, lstupd, indlst, 
									uelem, fmaskl, glloc, iops);
			_mtral.Ich2UpdateRow (nupd, lstupd, indlst, 
									lelem, fmasku, guloc, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		_mtral.Ich2UpdateTranspDynamic (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		_mtral.Ilu2FiltrRow (irow, _param.theta, _param.tau2, 
								nlist1, lstloc, lstupd, 
								fmaskl, fmasku, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin, irow, deleml, delemu, 
								nlist1, lstupd, fmaskl, fmasku, dpiv, iops);

// Sort the list of indices

		if (nlist1 != 0) qsort (lstupd, nlist1, sizeof(int), compint);

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist1 != 0) {
			iv[irow] = nzjg+1;
			madj[irow] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = irow;
		};

// Store the data in the global arrays

		_mtral.Ilu2StoreRowDynamic (irow, _param.tau1,
							nzjg, nzju, nzjr, nzjrl,
							deleml, delemu, 
							nlist1, lstupd, fmaskl, fmasku,
							ig, jg, gl, gu);

// Free R part if necessary

		if (nrowsfree > 0 && irow > 0 && irow % nrowsfree == 0 && irow != _mtral.n-1) {

			int msglev = (int) _param.msglev;

			Ilu2FreeRDynamic (msglev, irow+1,
								ibegfr, nzjg, nzju, nzjrl,
								iv, ig, jg, gl, gu,
								lstupd, fmaskl, fmasku);

		};

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

	itype = 2;
	nloc = _mtral.n;

	if (_param.msglev > 1) {
		csvhyst ("Ilu2Dia", nloc, dpiv);
	};

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] fmaskl;
	delete [] fmasku;
	delete [] dpiv;
	delete [] diagg;

// Perform the final filtering and scaling and compute matrix on return

	CSMatrixR temp (_mtral.n,nzju);

	_mtral.Ich2StoreScaleUDynamic (scll, ig, jg, gl, 
									temp, ops);

	_mtrl = temp;

	_mtral.Ich2StoreScaleUDynamic (sclu, ig, jg, gu, 
									temp, ops);

	_mtru = temp;

	temp = mtrdummy;

// Free factorization arrays 

	delete [] scll;
	delete [] sclu;
	delete [] ig;
	for (irow=0;irow<_mtral.n;irow++) {
		delete [] jg[irow];
		delete [] gl[irow];
		delete [] gu[irow];
	};
	delete [] jg;
	delete [] gl;
	delete [] gu;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ilu2 factorization statistics

	nz = _mtral.ia[_mtral.n];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	_param.density  = (double) nzju / (double) nz;
	_param.density2 = (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = 2*_mtral.ia[_mtral.n]-_mtral.n;

	ops = ops / (double) nz;

	if (_param.msglev > 1) {
		cout << " Ilu2 preconditioner generation statistics: " << endl;
		cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout << "     Costs = " << ops << " MvmA flops. " << endl;
		cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev >= 1) {
		_fout << " Ilu2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Second order Incomplete LU decomposition with dynamic memory allocation by blocks
//========================================================================================
void CSMatrixR::Ilu2DynamicByBlocks (ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition with dynamic memory allocation by blocks
										const CSMatrixR &_mtral, const CSMatrixR &_mtrau,
										CSMatrixR &_mtrl, CSMatrixR &_mtru) { 

//	const char *funcname = "Ilu2DynamicByBlocks";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	int nblks, iblk, nzallocmax, nzalloc, nzused;
	double lelem, uelem, deleml, delemu;

// Local arrays

	int *ig, **jg;
	double **gl, **gu;
	bool elmisr;
	int *blks, *nd2blk;
	int **jgblk;
	double **glblk, **gublk;

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	int iops;

	icheck = _param.ichkfct;
	ops = 0.0e0;

// Init time measurement

	time0 = clock ();

// Define working arrays

	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	double *fmaskl, *fmasku, *dpiv, *diagg, *scll, *sclu;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	_mtral.Ilu2AllocScale (_param.sclmin, fmaskl, fmasku, dpiv, diagg, scll, sclu);

// Output hystogram of the scaling values

	int itype = 1;
	int nloc = _mtral.n;

	if (_param.msglev > 1) {
		csvhyst ("SclDia", nloc, dpiv);
	};

// Allocate factorization arrays

	_mtral.Ilu2AllocGDynamicByBlocks (ig, jg, gl, gu, 
										blks, nd2blk, jgblk, glblk, gublk);

// Create dummy matrix

	CSMatrixR mtrdummy;

// Main cycle over the rows

//	int nrowstime = 10000;

	int jloc;
	int *jgloc;
	double *glloc, *guloc;

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;

	iblk = -1;
	nzalloc = 0;
	nzallocmax = (int) ( (double)_mtral.nzja * (-_param.memory) );
	if (nzallocmax < 0) nzallocmax = 0;
	nzused = 0;

	ig[0] = 0;
	blks[0] = 0;

// Set free R parameters

	int nfreer = (int) _param.nfreer;

	int nrowsfree;

	if (nfreer <= 0) {
		nrowsfree = 0;
	} else {
		nrowsfree = _mtral.n / nfreer;
	};

	time2 = clock ();
	time3 = clock ();

	for (irow=0;irow<_mtral.n;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {
			if (_param.msglev > 1) {
				cout << " Ilu2: Irow = " << irow << " Iblk = " << iblk << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << endl;
			};
		};

#ifdef __FlowVisionOut
		if (irow%nrowstime == 0) {

			time3 = clock ();

			tottim = (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			if (tottim > _param.timefctout) {

				nz = _mtral.ia[_mtral.n];

				densu = (double) nzju / (double) nz;
				densr = (double) nzjr / (double) nz;

				int prcent = 100*irow / _mtral.n;

				_param.meth->Err = densu;
				_param.meth->Residual = densr;
				_param.meth->N = prcent;
				if (_param.numpar) {
					for (int j=0; j<_param.numpar; j++) {
						_param.methj[j]->Err = _param.meth->Err;
						_param.methj[j]->Residual = _param.meth->Residual;
						_param.methj[j]->N   = _param.meth->N;
						ErrorOutput (_param.pregion, _param.methj[j]);
					};
				} else {
					ErrorOutput (_param.pregion, _param.meth);
				};

//				cout << " Ilu2 Fv: prcent = " << prcent << endl;
//				cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;

				time2 = clock ();

			};

		};
#endif

// Init current scaled row

		mtrdummy.Ilu2InitRow (_mtral, _mtrau,
							irow, nlist1, lstloc, icycle, imask, 
							diagg, fmaskl, fmasku, scll, sclu, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			jgloc = jg[irwprv];
			glloc = gl[irwprv];
			guloc = gu[irwprv];

			j = iv[irwprv];

			jloc = j-ig[irwprv];

			icolmn = jgloc[jloc];
			lelem = glloc[jloc];
			uelem = guloc[jloc];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			_mtral.Ilu2InitMaskDynamic (elmisr, irwprv, iv, ig, jg,
										nupd, lstupd, indlst,
										nlist1, lstloc, icycle, imask, 
										fmaskl, fmasku);

// Perform updates

			_mtral.Ich2UpdateRow (nupd, lstupd, indlst, 
									uelem, fmaskl, glloc, iops);
			_mtral.Ich2UpdateRow (nupd, lstupd, indlst, 
									lelem, fmasku, guloc, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		_mtral.Ich2UpdateTranspDynamic (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		_mtral.Ilu2FiltrRow (irow, _param.theta, _param.tau2, 
								nlist1, lstloc, lstupd, 
								fmaskl, fmasku, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin, irow, deleml, delemu, 
								nlist1, lstupd, fmaskl, fmasku, dpiv, iops);

// Sort the list of indices

		if (nlist1 != 0) qsort (lstupd, nlist1, sizeof(int), compint);

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist1 != 0) {
			iv[irow] = nzjg+1;
			madj[irow] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = irow;
		};

// Store the data in the global arrays

		_mtral.Ilu2StoreRowDynamicByBlocks (irow, _param.tau1,
											nzjg, nzju, nzjr, nzjrl,
											deleml, delemu, 
											nlist1, lstupd, fmaskl, fmasku,
											ig, jg, gl, gu,
											iblk, nzallocmax, nzalloc, nzused,
											blks, nd2blk, jgblk, glblk, gublk);

// Free R part if necessary

		if (nrowsfree > 0 && irow > 0 && irow % nrowsfree == 0 && irow != _mtral.n-1) {

			int msglev = (int) _param.msglev;

			Ilu2FreeRDynamicByBlocks (msglev, irow+1,
										ibegfr, nzjg, nzju, nzjrl,
										iv, ig, jg, gl, gu,
										iblk, nzalloc, nzused,
										blks, nd2blk, jgblk, glblk, gublk,
										lstupd, fmaskl, fmasku);

		};

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

	itype = 2;
	nloc = _mtral.n;

	if (_param.msglev > 1) {
		csvhyst ("Ilu2Dia", nloc, dpiv);
	};

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] fmaskl;
	delete [] fmasku;
	delete [] dpiv;
	delete [] diagg;

// Perform the final filtering and scaling and compute matrix on return

	nblks = iblk+1;

	CSMatrixR temp (_mtral.n,nzju);

	_mtral.Ich2StoreScaleUDynamic (scll, ig, jg, gl, 
									temp, ops);

	for (iblk=0;iblk<nblks;iblk++) {
		delete [] glblk[iblk];
	};

	_mtrl = temp;

	_mtral.Ich2StoreScaleUDynamic (sclu, ig, jg, gu, 
									temp, ops);

	for (iblk=0;iblk<nblks;iblk++) {
		delete [] gublk[iblk];
	};

	_mtru = temp;

	temp = mtrdummy;

// Free factorization arrays 

	delete [] scll;
	delete [] sclu;
	delete [] ig;
	delete [] jg;
	delete [] gl;
	delete [] gu;
	delete [] blks;
	delete [] nd2blk;
	for (iblk=0;iblk<nblks;iblk++) {
		delete [] jgblk[iblk];
	};
	delete [] jgblk;
	delete [] glblk;
	delete [] gublk;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ilu2 factorization statistics

	nz = _mtral.ia[_mtral.n];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	_param.density  = (double) nzju / (double) nz;
	_param.density2 = (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = 2*_mtral.ia[_mtral.n]-_mtral.n;

	ops = ops / (double) nz;

	if (_param.msglev > 1) {
		cout << " Ilu2 preconditioner generation statistics: " << endl;
		cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout << "     Costs = " << ops << " MvmA flops. " << endl;
		cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev >= 1) {
		_fout << " Ilu2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
