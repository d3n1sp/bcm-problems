//------------------------------------------------------------------------------------------------
// File: ich2r.cpp
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
// CSMatrix: Allocate and init integer data for Ich2/Ilu2
//========================================================================================
void CSMatrix::Ich2AllocInt (int &_icycle,int *&_imask, int *&_lstloc, // Allocate and init integer data for Ich2/Ilu2
							int *&_iv, int *&_madj, int *&_ibegm, int *&_lstupd, int *&_indlst) const { 

	const char *funcname = "Ich2AllocInt";

	_imask = new int [nsupr];
	if (!_imask) MemoryFail (funcname);
	_lstloc = new int [nsupr];
	if (!_lstloc) MemoryFail (funcname);
	_iv = new int [nsupr];
	if (!_iv) MemoryFail (funcname);
	_madj = new int [nsupr];
	if (!_madj) MemoryFail (funcname);
	_ibegm = new int [nsupr];
	if (!_ibegm) MemoryFail (funcname);
	_lstupd = new int [nsupr];
	if (!_lstupd) MemoryFail (funcname);
	_indlst = new int [nsupr];
	if (!_indlst) MemoryFail (funcname);

	_icycle = -1;

	for (int i=0;i<nsupr;i++) {
		_imask[i] = _icycle;
		_iv   [i] = -1;
		_madj [i] = -1;
		_ibegm[i] = -1;
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate double data and init scaling for Ich2
//========================================================================================
void CSMatrixR::Ich2AllocScale (int _msglev, double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
								double *&_diagg, double *&_scl) const {

	const char *funcname = "Ich2AllocScale";

	_fmask = new double [n];
	if (!_fmask) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_diagg = new double [n];
	if (!_diagg) MemoryFail (funcname);
	_scl = new double [n];
	if (!_scl) MemoryFail (funcname);

// Compute scaling of the matrix

	int i;

	for (i=0;i<n;i++) {
		double aux = sqrt(a[ia[i]]);
		_scl[i] = 1.0e0 / aux;
		_dpiv[i] = aux;
	};
//	for (i=0;i<n;i++) scl[i] = 1.0e0;

// Init working arrays

	for (i=0;i<n;i++) {
		_fmask[i] = 0.0e0;
		_diagg[i] = a[ia[i]]*_scl[i]*_scl[i];
	};

// Output hystogram of the scaling values

//	int itype = 1;
	int nloc = n;

	if (_msglev > 1) {
		csvhyst ("SclDia", nloc, _dpiv);
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate G data for Ich2
//========================================================================================
void CSMatrixR::Ich2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ich2
							int *&_ig, int *&_jg, int *&_addrg, double *&_g) const {

	const char *funcname = "Ich2AllocG";

	int nz = ia[nsupr];

	_nzjgmx = int (_memory * double(nz));
	_nzgmax = _nzjgmx;

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int [_nzjgmx];
	if (!_jg) MemoryFail (funcname);
	_addrg = new int [1];
	if (!_addrg) MemoryFail (funcname);
	_g = new double [_nzgmax];
	if (!_g) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate G data for Ich2
//========================================================================================
void CSMatrixR::Ich2AllocGDynamic (int *&_ig, int **&_jg, double **&_g) const { // Allocate G data for Ich2

	const char *funcname = "Ich2AllocGDynamic";

//	int nz = ia[nsupr];

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int * [nsupr+1];
	if (!_jg) MemoryFail (funcname);
	_g = new double * [nsupr+1];
	if (!_g) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixR: Allocate G data for Ich2
//========================================================================================
void CSMatrixR::Ich2AllocGDynamicByBlocks (int *&_ig, int **&_jg, double **&_g, // Allocate G data for Ich2
											int *&_blks, int *&_nd2blk, 
											int **&_jgblk, double **&_gblk) const {

	const char *funcname = "Ich2AllocGDynamic";

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int * [nsupr+1];
	if (!_jg) MemoryFail (funcname);
	_g = new double * [nsupr+1];
	if (!_g) MemoryFail (funcname);
	_blks = new int [nsupr+1];
	if (!_blks) MemoryFail (funcname);
	_nd2blk = new int [nsupr];
	if (!_nd2blk) MemoryFail (funcname);
	_jgblk = new int * [nsupr+1];
	if (!_jgblk) MemoryFail (funcname);
	_gblk = new double * [nsupr+1];
	if (!_gblk) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixR: Init current scaled row
//========================================================================================
void CSMatrixR::Ich2InitRow (int _irow, int &_nlist1, // Init current scaled row
							int *_lstloc, int _icycle, int *_imask, 
							double *_diagg, double *_fmask, double *_scl, int &_iops) const {

//	const char *funcname = "Ich2InitRow";

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	_fmask[_irow] = _diagg[_irow];
	_nlist1++;

	for (int j=ia[_irow]+1;j<ia[_irow+1];j++) {
		int jj = ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;
		_fmask[jj] = a[j]*_scl[_irow]*_scl[jj];
		_nlist1++;
	};

	_iops += _nlist1;

};

// Author: Kharchenko S.A.
// CSMatrixR: Init mask and update arrays
//========================================================================================
void CSMatrixR::Ich2InitMask (bool _elmisr, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmask) const {

//	const char *funcname = "Ich2InitMask";

	_nupd = 0;

	for (int k=_iv[_irwprv];k<_ig[_irwprv+1];k++) {

		int jcolmn = _jg[k];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (_imask[jcolmn] != _icycle) {
			_imask[jcolmn] = _icycle;
			_lstloc[_nlist1] = jcolmn;
			_fmask[jcolmn] = 0.0e0;
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
void CSMatrixR::Ich2InitMaskDynamic (bool _elmisr, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int **_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmask) const {

//	const char *funcname = "Ich2InitMaskDynamic";

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
			_fmask[jcolmn] = 0.0e0;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = kloc;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Update current row
//========================================================================================
void CSMatrixR::Ich2UpdateRow (int _nupd, int *_lstupd, int *_indlst, // Update current row
							double _lelem, double *_fmask, double *_g, int &_iops) const {

//	const char *funcname = "Ich2UpdateRow";

	for (int k=0;k<_nupd;k++) {
		int icol = _lstupd[k];
		int indg = _indlst[k];
		_fmask[icol] = _fmask[icol] - _lelem * _g[indg];
	};

	_iops += _nupd;

};

// Author: Kharchenko S.A.
// CSMatrix: Update arrays that support transposed structures search
//========================================================================================
void CSMatrix::Ich2UpdateTransp (int _irow, int *_ibegm, int *_madj, // Update arrays that support transposed structures search
								int *_iv, int *_ig, int *_jg) const {

//	const char *funcname = "Ich2UpdateTransp";

	int irwprv = _ibegm[_irow];

	while (irwprv != -1) {

		int irwpr1 = _madj[irwprv];

		if (_iv[irwprv] >= _ig[irwprv+1]-1) {

			_madj[irwprv] = -1;

		} else {

			int j = _iv[irwprv]+1;
			int icolmn = _jg[j];
			if (icolmn < 0) icolmn = -icolmn;
			_madj[irwprv] = _ibegm[icolmn];
			_ibegm[icolmn] = irwprv;

		};

		_iv[irwprv]++;

		irwprv = irwpr1;

	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Perform filtering of the current row
//========================================================================
void CSMatrixR::Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
							int &_nlist1, int *_lstloc, int *_lstupd,
							double *_fmask, double *_diagg, int &_iops) const {

//	const char *funcname = "Ich2FiltrRow";


	int nlist2 = 0;

	for (int i=1;i<_nlist1;i++) {

		int j = _lstloc[i];
		double aux = _fmask[j];
		if (aux < 0.0e0) aux = -aux;

		if (aux >= _tau2) {
			_lstupd[nlist2] = j;
			nlist2++;
		} else {
			aux *= _theta;
			_fmask[_irow] += aux;
			_diagg[j] += aux;
		};

	};

	_iops += 2*(_nlist1-nlist2);

	_nlist1 = nlist2;

};

// Author: Kharchenko S.A.
// CSMatrixR: Compute diagonal pivot and scale current row
//========================================================================
void CSMatrixR::Ich2PivScaleRow (int _msglev, int _irow, double &_delem, // Compute diagonal pivot and scale current row
								int _nlist1, int *_lstupd, double *_fmask, double *_dpiv, int &_iops) const {

//	const char *funcname = "Ich2PivScaleRow";

	double aux = _fmask[_irow];

	if (aux <= 0.0e0) {
		if (_msglev > 1) {
			cout << " Ich2: Irow = " << _irow << " Diagonal pivot is nonpositive";
		};
		throw " Ich2: Diagonal pivot is nonpositive";
	};

	aux = sqrt(aux);
	_dpiv[_irow] = aux;

	aux = 1.0e0 / aux;

	_delem = aux;

	_iops += 2;

// Update the data in the list

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		_fmask[j] *= aux;
	};

	_iops += _nlist1;

};

// Author: Kharchenko S.A.
// CSMatrixR: Free R part of the data
//========================================================================
void CSMatrixR::Ich2FreeR (int _msglev, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, double *_g) const {

//	const char *funcname = "Ich2FreeR";

	int iendfr = _irow;

	if (_msglev > 1) {
		cout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;
	};

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
				_g [_nzjg] = _g [j];
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

	if (_msglev > 1) {
		cout << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Free R part of the data
//========================================================================
void CSMatrixR::Ich2FreeRDynamic (int _msglev, int _irow, // Free R part of the data
								int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
								int *_iv, int *_ig, int **_jg, double **_g,
								int *_list, double *_fmask) const {

	const char *funcname = "Ich2FreeRDynamic";

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
	double *gloc;

	for (i=_ibegfr;i<iendfr;i++) {

// Initial scan of the current row

		jgloc = _jg[i];
		gloc = _g[i];

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
			_fmask[nlistloc] = gloc[jloc];
			nlistloc++;
			_nzjg++;

label:;

		};

// Free/allocate the memory

		delete [] _jg[i];
		delete [] _g[i];

		int *jgloc = new int [nlistloc];
		if (!jgloc) MemoryFail (funcname);
		double *gloc = new double [nlistloc];
		if (!gloc) MemoryFail (funcname);

// Final scan of the current row

		for (j=0;j<nlistloc;j++) {
			jgloc[j] = _list[j];
			gloc[j] = _fmask[j];
		};

		_jg[i] = jgloc;
		_g[i] = gloc;

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
//========================================================================
void CSMatrixR::Ich2FreeRDynamicByBlocks (int _msglev, int _irow, // Free R part of the data
										int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
										int *_iv, int *_ig, int **_jg, double **_g,
										int _iblk, int _nzalloc, int &_nzused,
										int *_blks, int *_nd2blk, int **_jgblk, double **_gblk,
										int *_list, double *_fmask) const {

	const char *funcname = "Ich2FreeRDynamicByBlocks";

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
	double *gloc;
	int *jgarr;
	double *garr;

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
		garr = new double [nzblk];
		if (!garr) MemoryFail (funcname);

		nzblk = 0;

		for (i=ibegnd;i<iendnd;i++) {

// Initial scan of the current row

			jgloc = _jg[i];
			gloc = _g[i];

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
				_fmask[nlistloc] = gloc[jloc];
				nlistloc++;
				_nzjg++;

label1:;

			};

// Final scan of the current row

			jgloc = jgarr+nzblk;
			gloc = garr+nzblk;

			for (j=0;j<nlistloc;j++) {
				jgloc[j] = _list[j];
				gloc[j] = _fmask[j];
			};

			_jg[i] = jgloc;
			_g[i] = gloc;

			ibeg = _ig[i+1];

			_ig[i+1] = _nzjg;

			if (_nzjrl == 0 && ibgfrn == i) ibgfrn = i+1;

// Update iv

			_iv[i] += _nzjg-ibeg;

			nzblk += nlistloc;

		};

// Free previous memory and store the new one

		delete [] _jgblk[iblkloc];
		delete [] _gblk[iblkloc];

		_jgblk[iblkloc] = jgarr;
		_gblk[iblkloc] = garr;

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
//========================================================================
void CSMatrixR::Ich2StoreRow (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double _delem, int _nlist1, int *_lstupd, double *_fmask, 
								int *_ig, int *_jg, double *_g) const {

//	const char *funcname = "Ich2StoreRow";

	_jg[_nzjg] = _irow;
	_g [_nzjg] = _delem;
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double aux1 = _fmask[j];
		if (aux1 < 0.0e0) aux1 = -aux1;
		if (aux1 >= _tau1) {
			_jg[_nzjg] = j;
			_g [_nzjg] = _fmask[j];
			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
			_g [_nzjg] = _fmask[j];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

};

// Author: Kharchenko S.A.
// CSMatrixR: Store current row
//========================================================================
void CSMatrixR::Ich2StoreRowDynamic (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double _delem, 
								int _nlist1, int *_lstupd, double *_fmask,
								int *_ig, int **_jg, double **_g) const {

	const char *funcname = "Ich2StoreRowDynamic";

	int *jgloc = new int [_nlist1+1];
	if (!jgloc) MemoryFail (funcname);
	double *gloc = new double [_nlist1+1];
	if (!gloc) MemoryFail (funcname);

	jgloc[0] = _irow;
	gloc[0] = _delem;
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double auxloc = _fmask[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		double aux = auxloc;

		if (aux >= _tau1) {
			jgloc[i+1] = j;
			gloc[i+1] = _fmask[j];
			_nzjg++;
			_nzju++;
		} else {
			jgloc[i+1] = -j;
			gloc[i+1] = _fmask[j];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

	_jg[_irow] = jgloc;
	_g[_irow] = gloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Store current row
//========================================================================
void CSMatrixR::Ich2StoreRowDynamicByBlocks (int _irow, double _tau1, // Store current row
												int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
												double _delem,
												int _nlist1, int *_lstupd, double *_fmask,
												int *_ig, int **_jg, double **_g,
												int &_iblk, int _nzallocmax, int &_nzalloc, int &_nzused,
												int *_blks, int *_nd2blk, int **_jgblk, 
												double **_gblk) const {

	const char *funcname = "Ich2StoreRowDynamicByBlocks";

	int *jgloc;
	double *gloc;

	if (_nlist1+1 > _nzalloc-_nzused) {
		_nzalloc = _nzallocmax;
		if (_nzalloc < _nlist1+1) _nzalloc = _nlist1+1;
		jgloc = new int [_nzalloc];
		if (!jgloc) MemoryFail (funcname);
		gloc = new double [_nzalloc];
		if (!gloc) MemoryFail (funcname);
		_iblk++;
		_jgblk[_iblk] = jgloc;
		_gblk[_iblk] = gloc;
		_nzused = 0;
	} else {
		jgloc = _jgblk[_iblk]+_nzused;
		gloc = _gblk[_iblk]+_nzused;
	};

	jgloc[0] = _irow;
	gloc[0] = _delem;
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double auxloc = _fmask[j];
		if (auxloc < 0.0e0) auxloc = -auxloc;
		double aux = auxloc;

		if (aux >= _tau1) {
			jgloc[i+1] = j;
			gloc[i+1] = _fmask[j];
			_nzjg++;
			_nzju++;
		} else {
			jgloc[i+1] = -j;
			gloc[i+1] = _fmask[j];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

	_jg[_irow] = jgloc;
	_g[_irow] = gloc;

	_blks[_iblk+1] = _irow+1;
	_nd2blk[_irow] = _iblk;
	_nzused += (_nlist1+1);

};

// Author: Kharchenko S.A.
// CSMatrixR: Store the resulting U and scale it
//========================================================================
void CSMatrixR::Ich2StoreScaleU (double *_scl, // Store the resulting U and scale it
								int *_ig, int *_jg, double *_g,
								CSMatrixR &_temp, double &_ops) const {

//	const char *funcname = "Ich2StoreScaleU";

	_temp.m = n;
	_temp.n = n;
	_temp.nsupr = nsupr;
	_temp.nsupc = nsupc;
	_temp.nlist = nsupr;

	int i;

	for (i=0;i<nsupr;i++) _temp.list[i] = i;

	_temp.ia[0] = 0;

	int nz = 0;

	for (i=0;i<nsupr;i++) {
		for (int j=_ig[i];j<_ig[i+1];j++) {
			int jj = _jg[j];
			if (jj >= 0) {
				_temp.ja[nz] = jj;
				if (jj == i) {
					_temp.a [nz] = _g[j] * _scl[i];
				} else {
					_temp.a [nz] = _g[j] / _scl[jj];
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
// CSMatrixR: Second order Incomplete Cholessky decomposition
//========================================================================
CSMatrixR CSMatrixR::Ich2 (CSlvParam &_param) const { // Second order Incomplete Cholessky decomposition

//	const char *funcname = "Ich2";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzgmax, nzjgmx, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	double lelem, delem;
	int *ig, *jg, *addrg;
	double *g;
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
	double *fmask, *dpiv, *diagg, *scl;

// Allocate and init working arrays including the scaling

// Integer arrays

	Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	Ich2AllocScale (_param.msglev, fmask, dpiv, diagg, scl);

// Allocate factorization arrays

	Ich2AllocG (_param.memory, nzjgmx, nzgmax, ig, jg, addrg, g);

// Main cycle over the rows

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;

	ig[0] = 0;

	for (irow=0;irow<n;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (_param.msglev > 1) {
			if (irow > 0 && irow%icheck==0) {
				cout << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			};
		};

// Init current scaled row

		Ich2InitRow (irow, nlist1, lstloc, icycle, imask, 
					diagg, fmask, scl, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];
			lelem = g[j];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMask (elmisr, irwprv, iv, ig, jg, 
							nupd, lstupd, indlst,
							nlist1, lstloc, icycle, imask, fmask);

// Perform updates

			Ich2UpdateRow (nupd, lstupd, indlst, 
							lelem, fmask, g, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTransp (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow (irow, _param.theta, _param.tau2, 
						nlist1, lstloc, lstupd, fmask, diagg, iops);

// Compute and store the pivot value and scale current row

		Ich2PivScaleRow (_param.msglev, irow, delem, nlist1, lstupd, fmask, dpiv, iops);

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx) {

// Free the memory containing R part that is not in use any more

			Ich2FreeR (_param.msglev, irow, 
						ibegfr, nzjg, nzju, nzjrl, 
						iv, ig, jg, g);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx) {
			if (_param.msglev > 1) {
				cout << " Insufficient workspace in Ich2, Irow = " << irow << endl;
			};
			throw " Insufficient workspace in Ich2";
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

		Ich2StoreRow (irow, _param.tau1,
						nzjg, nzju, nzjr, nzjrl,
						delem, nlist1, lstupd, fmask, 
						ig, jg, g);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = n;

	if (_param.msglev > 1) {
		csvhyst ("Ich2Dia", nloc, dpiv);
	};

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] fmask;
	delete [] dpiv;
	delete [] diagg;

// Perform the final filtering and scaling and compute matrix on return

	CSMatrixR temp (n,nzju);

	Ich2StoreScaleU (scl, ig, jg, g, 
					temp, ops);

// Free factorization arrays 

	delete [] scl;
	delete [] ig;
	delete [] jg;
	delete [] addrg;
	delete [] g;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ich2 factorization statistics

	nz = ia[n];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	_param.density  = (double) nzju / (double) nz;
	_param.density2 = (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = 2*ia[n]-n;

	ops = ops / (double) nz;

	if (_param.msglev > 1) {
		cout << " Ich2 preconditioner generation statistics: " << endl;
		cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout << "     Costs = " << ops << " MvmA flops. " << endl;
		cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

// Return the result

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Second order Incomplete Cholessky decomposition with dynamic memory allocation
//========================================================================
CSMatrixR CSMatrixR::Ich2Dynamic (ofstream &_fout, CSlvParam &_param) const { // Second order Incomplete Cholessky decomposition with dynamic memory allocation

//	const char *funcname = "Ich2Dynamic";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	double lelem, delem;
	int *ig, **jg;
	double **g;
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
	double *fmask, *dpiv, *diagg, *scl;

// Allocate and init working arrays including the scaling

// Integer arrays

	Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	Ich2AllocScale (_param.msglev, fmask, dpiv, diagg, scl);

// Allocate factorization arrays

	Ich2AllocGDynamic (ig, jg, g);

// Set free R parameters

	int nfreer = (int) _param.nfreer;

	int nrowsfree;

	if (nfreer <= 0) {
		nrowsfree = 0;
	} else {
		nrowsfree = n / nfreer;
	};

// Main cycle over the rows

	int *jgloc;
	double *gloc;

	int jloc;

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;

	ig[0] = 0;

	for (irow=0;irow<n;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (_param.msglev > 1) {
			if (irow > 0 && irow%icheck==0) {
				cout << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << endl;
			};
		};

// Init current scaled row

		Ich2InitRow (irow, nlist1, lstloc, icycle, imask, 
					diagg, fmask, scl, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			jgloc = jg[irwprv];
			gloc = g[irwprv];

			j = iv[irwprv];

			jloc = j-ig[irwprv];

			icolmn = jgloc[jloc];
			lelem = gloc[jloc];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMaskDynamic (elmisr, irwprv, iv, ig, jg, 
								nupd, lstupd, indlst,
								nlist1, lstloc, icycle, imask, fmask);

// Perform updates

			Ich2UpdateRow (nupd, lstupd, indlst, 
							lelem, fmask, gloc, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTranspDynamic (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow (irow, _param.theta, _param.tau2, 
						nlist1, lstloc, lstupd, fmask, diagg, iops);

// Compute and store the pivot value and scale current row

		Ich2PivScaleRow (_param.msglev, irow, delem, nlist1, lstupd, fmask, dpiv, iops);

// Sort the list of indices

		if (nlist1 != 0) qsort (lstupd, nlist1, sizeof(int), compint);

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist1 != 0) {
			iv[irow] = nzjg+1;
			madj[irow] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = irow;
		};

// Store the data in the global arrays

		Ich2StoreRowDynamic (irow, _param.tau1,
								nzjg, nzju, nzjr, nzjrl,
								delem, nlist1, lstupd, fmask, 
								ig, jg, g);

// Free R part if necessary

		if (nrowsfree > 0 && irow > 0 && irow % nrowsfree == 0 && irow != n-1) {

			int msglev = (int) _param.msglev;

			Ich2FreeRDynamic (msglev, irow+1,
								ibegfr, nzjg, nzju, nzjrl,
								iv, ig, jg, g,
								lstupd, fmask);

		};

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = n;

	if (_param.msglev > 1) {
		csvhyst ("Ich2Dia", nloc, dpiv);
	};

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] fmask;
	delete [] dpiv;
	delete [] diagg;

// Perform the final filtering and scaling and compute matrix on return

	CSMatrixR temp (n,nzju);

	Ich2StoreScaleUDynamic (scl, ig, jg, g, 
							temp, ops);

// Free factorization arrays 

	delete [] scl;
	delete [] ig;
	for (irow=0;irow<n;irow++) {
		delete [] jg[irow];
		delete [] g[irow];
	};
	delete [] jg;
	delete [] g;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ich2 factorization statistics

	nz = ia[n];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	_param.density  = (double) nzju / (double) nz;
	_param.density2 = (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = 2*ia[n]-n;

	ops = ops / (double) nz;

	if (_param.msglev > 1) {
		cout << " Ich2 preconditioner generation statistics: " << endl;
		cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout << "     Costs = " << ops << " MvmA flops. " << endl;
		cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev >= 1) {
		_fout << " Ich2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

// Return the result

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixR: Second order Incomplete Cholessky decomposition with dynamic memory allocation by blocks
//========================================================================
void CSMatrixR::Ich2DynamicByBlocks (ofstream &_fout, CSlvParam &_param,
													int _istore, CGSMatrixR &_gich, CSMatrixR &_ich) const { // Second order Incomplete Cholessky decomposition with dynamic memory allocation by blocks

//	const char *funcname = "Ich2DynamicByBlocks";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	int nblks, iblk, nzallocmax, nzalloc, nzused;
	double lelem, delem;
	int *ig, **jg;
	double **g;
	bool elmisr;
	int *blks, *nd2blk;
	int **jgblk;
	double **gblk;

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
	double *fmask, *dpiv, *diagg, *scl;

// Allocate and init working arrays including the scaling

// Integer arrays

	Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	Ich2AllocScale (_param.msglev, fmask, dpiv, diagg, scl);

// Allocate factorization arrays

	Ich2AllocGDynamicByBlocks (ig, jg, g,
								blks, nd2blk, jgblk, gblk);

// Set free R parameters

	int nfreer = (int) _param.nfreer;

	int nrowsfree;

	if (nfreer <= 0) {
		nrowsfree = 0;
	} else {
		nrowsfree = n / nfreer;
	};

// Main cycle over the rows

//	int nrowstime = 10000;

	int *jgloc;
	double *gloc;

	int jloc;

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;

	iblk = -1;
	nzalloc = 0;
	nzallocmax = (int) ( (double) nzja * (-_param.memory) );
	if (nzallocmax < 0) nzallocmax = 0;
	nzused = 0;

	ig[0] = 0;
	blks[0] = 0;

	time2 = clock ();
	time3 = clock ();

	for (irow=0;irow<n;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (_param.msglev > 1) {
			if (irow > 0 && irow%icheck==0) {
				cout << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << endl;
			};
		};

#ifdef __FlowVisionOut
		if (irow%nrowstime == 0) {

			time3 = clock ();

			tottim = (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			if (tottim > _param.timefctout) {

				nz = ia[n];

				densu = (double) nzju / (double) nz;
				densr = (double) nzjr / (double) nz;

				int prcent = 100*irow / n;

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

//				cout << " Ich2 Fv: prcent = " << prcent << endl;
//				cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;

				time2 = clock ();

			};

		};
#endif

// Init current scaled row

		Ich2InitRow (irow, nlist1, lstloc, icycle, imask, 
					diagg, fmask, scl, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			jgloc = jg[irwprv];
			gloc = g[irwprv];

			j = iv[irwprv];

			jloc = j-ig[irwprv];

			icolmn = jgloc[jloc];
			lelem = gloc[jloc];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMaskDynamic (elmisr, irwprv, iv, ig, jg, 
								nupd, lstupd, indlst,
								nlist1, lstloc, icycle, imask, fmask);

// Perform updates

			Ich2UpdateRow (nupd, lstupd, indlst, 
							lelem, fmask, gloc, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTranspDynamic (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow (irow, _param.theta, _param.tau2, 
						nlist1, lstloc, lstupd, fmask, diagg, iops);

// Compute and store the pivot value and scale current row

		Ich2PivScaleRow (_param.msglev, irow, delem, nlist1, lstupd, fmask, dpiv, iops);

// Sort the list of indices

		if (nlist1 != 0) qsort (lstupd, nlist1, sizeof(int), compint);

// Add data of the current row into the arrays iv, madj and ibegm

		if (nlist1 != 0) {
			iv[irow] = nzjg+1;
			madj[irow] = ibegm[lstupd[0]];
			ibegm[lstupd[0]] = irow;
		};

// Store the data in the global arrays

		Ich2StoreRowDynamicByBlocks (irow, _param.tau1,
										nzjg, nzju, nzjr, nzjrl,
										delem, nlist1, lstupd, fmask,
										ig, jg, g,
										iblk, nzallocmax, nzalloc, nzused,
										blks, nd2blk, jgblk, gblk);

// Free R part if necessary

		if (nrowsfree > 0 && irow > 0 && irow % nrowsfree == 0 && irow != n-1) {

			int msglev = (int) _param.msglev;

			Ich2FreeRDynamicByBlocks (msglev, irow+1,
										ibegfr, nzjg, nzju, nzjrl,
										iv, ig, jg, g,
										iblk, nzalloc, nzused,
										blks, nd2blk, jgblk, gblk,
										lstupd, fmask);

		};

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = n;

	if (_param.msglev > 1) {
		csvhyst ("Ich2Dia", nloc, dpiv);
	};

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] fmask;
	delete [] dpiv;
	delete [] diagg;

// Perform the final filtering and scaling and compute matrix on return

	nblks = iblk+1;

	if (_istore == 0) {

		CSMatrixR temp (n,nzju);

		Ich2StoreScaleUDynamic (scl, ig, jg, g, 
								temp, ops);

		for (iblk=0;iblk<nblks;iblk++) {
			delete [] jgblk[iblk];
			delete [] gblk[iblk];
		};

		_ich = temp;

		CSMatrixR mtrdummy;

		temp = mtrdummy;

	} else {

		CGSMatrixR gtemp (nblks,n,nblks);

		_gich = gtemp;

		CGSMatrixR::Ich2StoreScaleUDynamic (nblks, blks, scl, ig, jgblk, gblk, 
									_gich);

		for (iblk=0;iblk<nblks;iblk++) {
			delete [] jgblk[iblk];
			delete [] gblk[iblk];
		};

	};

// Free factorization arrays 

	delete [] scl;
	delete [] ig;
//	for (irow=0;irow<n;irow++) {
//		delete [] jg[irow];
//		delete [] g[irow];
//	};
	delete [] jg;
	delete [] g;
	delete [] blks;
	delete [] nd2blk;
	delete [] jgblk;
	delete [] gblk;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ich2 factorization statistics

	nz = ia[n];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	_param.density  = (double) nzju / (double) nz;
	_param.density2 = (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = 2*ia[n]-n;

	ops = ops / (double) nz;

	if (_param.msglev > 1) {
		cout << " Ich2 preconditioner generation statistics: " << endl;
		cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout << "     Costs = " << ops << " MvmA flops. " << endl;
		cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (_param.msglev >= 1) {
		_fout << " Ich2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Incomplete Cholessky decomposition
//========================================================================
CSMatrixR CSMatrixR::Ich () const { // Incomplete Cholessky decomposition

//	const char *funcname = "Ich";

	cout << " Compute incomplete Cholessky, N =" << n << endl;

	clock_t time0, time1;
	double tottim;

// Init time measurement

	time0 = clock ();

// Compute transposed matrix

	CSMatrixR at;

	at = this->TranspMtr ();

// Compute Ich

	CSMatrixR ich(n,nzja);

// Copy sparsity structure

	int i;

	for (i=0;i<=n;i++) ich.ia[i] = at.ia[i];
	for (i=0;i<nzja;i++) ich.ja[i] = at.ja[i];

	double r;

// Formation of subdiagonal elements:

	for (i=0; i<n; i++)
	{
		int i1   = at.ia[i];
		int i2   = at.ia[i+1]-1;
		int j;
		for (j=i1; j<i2; j++)
		{
			r       = at.a[j];
			int k   = at.ja[j];
			int ik  = at.ia[k];
			int ik2 = at.ia[k+1]-1;
			int it  = i1;
			while (ik<ik2 && it<j)
			{
				if (at.ja[ik] < at.ja[it])
				{
					ik++;
				}
				else if (at.ja[ik] > at.ja[it])
				{
					it++;
				}
				else
				{
	             r -= ich.a[ik] * ich.a[it];
		          it++;
			       ik++;
				}
			}
			ich.a[j] = r * ich.a[ik2];
		}
// Diagonal 
		r = at.a[i2];
		for (j=i1; j<i2; j++)
		{
			r -= ich.a[j] * ich.a[j];
		}
		if (r < 1.0e-20) r = 1.0e-20;
		ich.a[i2] = 1.0e0 / sqrt (r);
	};

// Transpose the matrix back

	CSMatrixR icht;

	icht = ich.TranspMtr ();

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	cout << "     Fct Time  = " << tottim << " sec." << endl;

	return icht;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
