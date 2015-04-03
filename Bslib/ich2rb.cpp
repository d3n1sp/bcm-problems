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

#include "smatrix.h"
#include "slvparam.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixRB: Allocate double data and init scaling for Ich2
//========================================================================================
void CSMatrixRB::Ich2AllocScale (double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
								double *&_diagg, double *&_scl, double *&_sclinv,
								double *&_lelem, double *&_delem,
								double *&_aloc, double *&_vloc, 
								double *&_eig, int &_lwork, double *&_work,
								double &_ops) const {

	const char *funcname = "Ich2AllocScale";

	_fmask = new double [n*bla];
	if (!_fmask) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_diagg = new double [n*bla];
	if (!_diagg) MemoryFail (funcname);
	_scl = new double [n*bla];
	if (!_scl) MemoryFail (funcname);
	_sclinv = new double [n*bla];
	if (!_sclinv) MemoryFail (funcname);
	_lelem = new double [bla_2];
	if (!_lelem) MemoryFail (funcname);
	_delem = new double [bla_2];
	if (!_delem) MemoryFail (funcname);

// Allocate working arrays

	int info=0;

	_lwork = 10 * bla;

	_aloc = new double [bla_2];
	if (!_aloc) MemoryFail (funcname);
	_vloc = new double [bla_2];
	if (!_vloc) MemoryFail (funcname);
	_eig = new double [bla];
	if (!_eig) MemoryFail (funcname);
	_work = new double [_lwork];
	if (!_work) MemoryFail (funcname);

// Compute scaling and its inverse

	int j, isup, kii, kjj, kkk, ibs, ibssc;
	double aux;

	for (isup=0;isup<nsupr;isup++) {

// Assign diagonal supernode

		j = ia[isup];

		for (kii=0;kii<bla_2;kii++) {
			_aloc[kii] = a[j*bla_2+kii];
		};

// Compute its symmetric spectral decomposition

		int blal = bla;

		dsyev_ ("V", "U", &blal, _aloc, &blal, _eig, _work, &_lwork, &info);

		if (info != 0) throw " Error in the Lapack routine DSYEV";

		_ops += 8 * bla_3;

// Check for negative pivots

		for (kii=0;kii<bla;kii++) {
			if (_eig[kii] < 0.0e0) throw " Negative eigenvalue when scaling ";
		};

// Store scaling factors

		for (kii=0;kii<bla;kii++) {
			aux = _eig[kii];
			aux = sqrt(aux);
			_eig[kii] = aux;
			_dpiv[isup*bla+kii] = _eig[kii];
		};

		ibs = isup*bla_2;

		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				_sclinv[ibs+kjj*bla+kii] = _aloc[kjj*bla+kii] * _eig[kjj];
			};
		};

		for (kii=0;kii<bla;kii++) {
			aux = _eig[kii];
			aux = 1.0e0 / aux;
			_eig[kii] = aux;
		};

		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				_scl[ibs+kjj*bla+kii] = _aloc[kii*bla+kjj] * _eig[kii];
			};
		};

		_ops += 2*bla_2 + 2*bla;

	};

// Reassign scaling factors if necessary

	if (false) {

		for (kii=0;kii<n*bla;kii++) _scl[kii] = 0.0e0;
		for (kii=0;kii<n*bla;kii++) _sclinv[kii] = 0.0e0;

		for (isup=0;isup<nsupr;isup++) {
			for (kii=0;kii<bla;kii++) {
				_scl[isup*bla_2+kii*bla+kii] = 1.0e0;
				_sclinv[isup*bla_2+kii*bla+kii] = 1.0e0;
			};
		};

	};

// Init working arrays

	for (kii=0;kii<n*bla;kii++) _fmask[kii] = 0.0e0;

	for (isup=0;isup<nsupr;isup++) {
		j = ia[isup];
		ibs = j*bla_2;
		ibssc = isup*bla_2;
		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				aux = 0.0e0;
				for (kkk=0;kkk<bla;kkk++) {
					aux += a[ibs+kjj*bla+kkk] * _scl[ibssc+kkk*bla+kii];
				};
				_aloc[kjj*bla+kii] = aux;
			};
		};
		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				aux = 0.0e0;
				for (kkk=0;kkk<bla;kkk++) {
					aux += _aloc[kkk*bla+kii] * _scl[ibssc+kkk*bla+kjj];
				};
				_vloc[kjj*bla+kii] = aux;
			};
		};
		for (kii=0;kii<bla_2;kii++) _diagg[ibssc+kii] = _vloc[kii];
	};

	_ops += 2 * nsupr * bla_3;

// Output hystogram of the scaling values

	int itype = 1;
	int nloc = n;

	csvhyst (cout, "ScDia",nloc, _dpiv);
//	csvhyst (_fout,"ScDia",nloc, _dpiv);

};

// Author: Kharchenko S.A.
// CSMatrixRB: Allocate G data for Ich2
//========================================================================================
void CSMatrixRB::Ich2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ich2
							int *&_ig, int *&_jg, int *&_addrg, double *&_g) const {

	const char *funcname = "Ich2AllocG";

	int nz = ia[nsupr];

	_nzjgmx = int (_memory * double(nz));
	_nzgmax = _nzjgmx*bla_2;

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
// CSMatrixRB: Init current scaled row
//========================================================================================
void CSMatrixRB::Ich2InitRow (int _irow, int &_nlist1, // Init current scaled row
							int *_lstloc, int _icycle, int *_imask, 
							double *_diagg, double *_fmask, double *_scl, 
							double *_aloc, double *_vloc,
							int &_iops) const {

	const char *funcname = "Ich2InitRow";

	int kii, kjj, kkk, ibs, ibssc, jbssc;
	int kjjb;
	double aux;

	double *paloc;
	double *pscloc;

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	for (kii=0;kii<bla_2;kii++) {
		_fmask[_irow*bla_2+kii] = _diagg[_irow*bla_2+kii];
	};
	_nlist1++;

	for (int j=ia[_irow]+1;j<ia[_irow+1];j++) {
		int jj = ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;

// Local scaling

		ibs = j*bla_2;
		ibssc = _irow*bla_2;
		jbssc = jj*bla_2;
		for (kii=0;kii<bla;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<bla;kjj++) {
				paloc = a+ibs+kjjb;
				pscloc = _scl+ibssc+kii;
				aux = 0.0e0;
				for (kkk=0;kkk<bla;kkk++) {
//					aux += a[ibs+kjjb+kkk] * _scl[ibssc+kkk*bla+kii];
					aux += *paloc * *pscloc;
					paloc++;
					pscloc += bla;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += bla;
			};
		};
		for (kii=0;kii<bla;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<bla;kjj++) {
				aux = 0.0e0;
				paloc = _aloc+kii;
				pscloc = _scl+jbssc+kjj;
				for (kkk=0;kkk<bla;kkk++) {
//					aux += _aloc[kkk*bla+kii] * _scl[jbssc+kkk*bla+kjj];
					aux += *paloc * *pscloc;
					pscloc += bla;
					paloc += bla;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += bla;
			};
		};
		for (kii=0;kii<bla_2;kii++) _fmask[jbssc+kii] = _vloc[kii];

		_nlist1++;

	};

	_iops += 2*_nlist1*bla_3;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Init mask and update arrays
//========================================================================================
void CSMatrixRB::Ich2InitMask (bool _elmisr, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmask) const {

	const char *funcname = "Ich2InitMask";

	int kii;

	_nupd = 0;

	for (int k=_iv[_irwprv];k<_ig[_irwprv+1];k++) {

		int jcolmn = _jg[k];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (_imask[jcolmn] != _icycle) {
			_imask[jcolmn] = _icycle;
			_lstloc[_nlist1] = jcolmn;
			for (kii=0;kii<bla_2;kii++) _fmask[jcolmn*bla_2+kii] = 0.0e0;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = k;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixRB: Update current row
//========================================================================================
void CSMatrixRB::Ich2UpdateRow (int _nupd, int *_lstupd, int *_indlst, // Update current row
							double *_lelem, double *_fmask, double *_g, int &_iops) const {

	const char *funcname = "Ich2UpdateRow";

	int kii, kjj, ibs, jbs;
	int kiib, kjjb;

	double *pfmloc;
	double *plloc;
	double *pgloc;

	for (int k=0;k<_nupd;k++) {
		int icol = _lstupd[k];
		int indg = _indlst[k];
		ibs = icol*bla_2;
		jbs = indg*bla_2;
		kiib = 0;
		for (kii=0;kii<bla;kii++) {
			pfmloc = _fmask+ibs+kii;
			kjjb = 0;
			for (kjj=0;kjj<bla;kjj++) {
				plloc = _lelem+kiib;
				pgloc = _g+jbs+kjjb;
				double *end = plloc+bla;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*bla+kii] -= _lelem[kii*bla+kkk] * _g[jbs+kjj*bla+kkk];
					*pfmloc -= *plloc * *pgloc;
				};
				pfmloc += bla;
				kjjb += bla;
			};
			kiib += bla;
		};
	};

	_iops += _nupd*bla_3;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Perform filtering of the current row
//========================================================================================
void CSMatrixRB::Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
								int &_nlist1, int *_lstloc, int *_lstupd,
								double *_fmask, double *_diagg, int &_iops) const {

	const char *funcname = "Ich2FiltrRow";

	int nlist2 = 0;
	int kii, kjj;

	for (int i=1;i<_nlist1;i++) {

		int j = _lstloc[i];

		double aux = 0.0e0;

		for (kii=0;kii<bla_2;kii++) {
			double auxloc = _fmask[j*bla_2+kii];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};

		if (aux >= _tau2) {
			_lstupd[nlist2] = j;
			nlist2++;
		} else {
			for (kii=0;kii<bla;kii++) {
				for (kjj=0;kjj<bla;kjj++) {
					double auxloc = _fmask[j*bla_2+kjj*bla+kii];
					if (auxloc < 0.0e0) auxloc = -auxloc;
					auxloc *= _theta;
					_fmask[_irow*bla_2+kii*bla+kii] += auxloc;
					_diagg[j*bla_2+kjj*bla+kjj] += auxloc;
				};
			};
		};

	};

	_nlist1 = nlist2;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Compute diagonal pivot and scale current row
//========================================================================================
void CSMatrixRB::Ich2PivScaleRow (int _irow, double *_delem, // Compute diagonal pivot and scale current row
								int _nlist1, int *_lstupd, double *_fmask, double *_dpiv, 
								double *_aloc, double *_eig, int _lwork, double *_work,
								int &_iops) const {

	const char *funcname = "Ich2PivScaleRow";

	int kii, kjj, kkk, ibs;
	double aux;

// Copy diagonal supernode

	ibs = _irow * bla_2;

	for (kii=0;kii<bla_2;kii++) _aloc[kii] = _fmask[ibs+kii];

// Compute spectral decomposition of the diagonal supernode

	int blal = bla;
	int info = 0;

	dsyev_ ("V", "U", &blal, _aloc, &blal, _eig, _work, &_lwork, &info);

	if (info != 0) throw " Error in the Lapack routine DSYEV";

	_iops += 8 * bla_3;

// Check for negative pivots

	for (kii=0;kii<bla;kii++) {
		if (_eig[kii] < 0.0e0) throw " Negative pivot in ICH2 ";
		_eig[kii] = sqrt(_eig[kii]);
		_dpiv[_irow*bla+kii] = _eig[kii];
		_eig[kii] = 1.0e0/_eig[kii];
	};

// Compute diagonal supernode

	for (kii=0;kii<bla;kii++) {
		for (kjj=0;kjj<bla;kjj++) {
			_delem[kjj*bla+kii] = _aloc[kii*bla+kjj] * _eig[kii];
		};
	};

	_iops += bla_2;

// Perform scaling of the data in the list

	double *pdeloc;
	double *pfmloc;
	double *palloc;

	int kjjb;

	for (int i=0;i<_nlist1;i++) {

		int j = _lstupd[i];
		ibs = j*bla_2;

		for (kii=0;kii<bla;kii++) {
			palloc = _aloc+kii;
			for (kjj=0;kjj<bla;kjj++) {
				kjjb = kjj*bla;
				pfmloc = _fmask+ibs+kjjb;
				pdeloc = _delem+kii;
				aux = 0.0e0;
				for (kkk=0;kkk<bla;kkk++) {
//					aux += _delem[kkk*bla+kii] * _fmask[ibs+kjj*bla+kkk];
					aux += *pdeloc * *pfmloc;
					pfmloc++;
					pdeloc += bla;
				};
//				_aloc[kjj*bla+kii] = aux;
				*palloc = aux;
				palloc += bla;
			};
		};

		for (kii=0;kii<bla_2;kii++) _fmask[ibs+kii] = _aloc[kii];

	};

	_iops += _nlist1*bla_3;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Free R part of the data
void CSMatrixRB::Ich2FreeR (int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, double *_g) const {

	const char *funcname = "Ich2FreeR";

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
				for (int kii=0;kii<bla_2;kii++) _g [_nzjg*bla_2+kii] = _g [j*bla_2+kii];
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
// CSMatrixRB: Store current row
//========================================================================================
void CSMatrixRB::Ich2StoreRow (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl,
								double *_delem, int _nlist1, int *_lstupd, double *_fmask, 
								int *_ig, int *_jg, double *_g) const {

	const char *funcname = "Ich2StoreRow";

	int kii, kjj;

	_jg[_nzjg] = _irow;
	for (kii=0;kii<bla;kii++) {
		for (kjj=0;kjj<bla;kjj++) {
			_g [_nzjg*bla_2+kjj*bla+kii] = _delem[kii*bla+kjj];
		};
	};
	_nzjg++;
	_nzju++;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		double aux = 0.0e0;
		for (kii=0;kii<bla_2;kii++) {
			double auxloc = _fmask[j*bla_2+kii];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		if (aux >= _tau1) {
			_jg[_nzjg] = j;
			for (kii=0;kii<bla_2;kii++) _g [_nzjg*bla_2+kii] = _fmask[j*bla_2+kii];
			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
			for (kii=0;kii<bla_2;kii++) _g [_nzjg*bla_2+kii] = _fmask[j*bla_2+kii];
			_nzjg++;
			_nzjr++;
			_nzjrl++;
		};
	};

	_ig[_irow+1] = _nzjg;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Store the resulting U and scale it
//========================================================================================
void CSMatrixRB::Ich2StoreScaleU (double *_scl, double *_sclinv, // Store the resulting U and scale it
								int *_ig, int *_jg, double *_g,
								CSMatrixRB &_temp, 
								double *_aloc,
								double &_ops) const {

	const char *funcname = "Ich2StoreScaleU";

	_temp.m = n;
	_temp.n = n;
	_temp.nsupr = nsupr;
	_temp.nsupc = nsupc;
	_temp.nlist = nsupr;

	int i;

	for (i=0;i<nsupr;i++) _temp.list[i] = i;

	_temp.ia[0] = 0;

	int kii, kjj, kkk, kiib, kjjb;
	int ibs, jbs, jjbs;
	double aux;

	double *pscloc;
	double *pgloc;
	double *paloc;

	int nz = 0;

	for (i=0;i<nsupr;i++) {
		ibs = i*bla_2;
		for (int j=_ig[i];j<_ig[i+1];j++) {
			jbs = j*bla_2;
			int jj = _jg[j];
			jjbs = jj*bla_2;
			if (jj >= 0) {
				_temp.ja[nz] = jj;
				if (jj == i) {
					kiib = 0;
					for (kii=0;kii<bla;kii++) {
						kjjb = 0;
						paloc = _aloc+kii;
						for (kjj=0;kjj<bla;kjj++) {
							pgloc = _g+jbs+kjjb;
							pscloc = _scl+ibs+kiib;
							aux = 0.0e0;
							for (kkk=0;kkk<bla;kkk++) {
//								aux += _scl[ibs+kii*bla+kkk] * _g[jbs+kjj*bla+kkk];
								aux += *pscloc * *pgloc;
								pgloc++;
								pscloc++;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = aux;
							paloc += bla;
							kjjb += bla;
						};
						kiib += bla;
					};
					for (kii=0;kii<bla_2;kii++) _temp.a[nz*bla_2+kii] = _aloc[kii];
				} else {
					for (kii=0;kii<bla;kii++) {
						kjjb = 0;
						paloc = _aloc+kii;
						for (kjj=0;kjj<bla;kjj++) {
							aux = 0.0e0;
							pgloc = _g+jbs+kii;
							pscloc = _sclinv+jjbs+kjj;
							for (kkk=0;kkk<bla;kkk++) {
//								aux += _g[jbs+kkk*bla+kii] * _sclinv[jjbs+kkk*bla+kjj];
								aux += *pgloc * *pscloc;
								pgloc += bla;
								pscloc += bla;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = aux;
							paloc += bla;
							kjjb += bla;
						};
					};
					for (kii=0;kii<bla_2;kii++) _temp.a[nz*bla_2+kii] = _aloc[kii];
				};
				nz++;
			};
		};
		_temp.ia[i+1] = nz;
	};

	_temp.nzja = nz;
	_temp.nza = nz*bla_2;

	_ops += nz*bla_3;

};

// Author: Kharchenko S.A.
// CSMatrixRB: Second order Incomplete Cholessky decomposition
//========================================================================================
CSMatrixRB CSMatrixRB::Ich2 (const CSlvParam _param) const { // Second order Incomplete Cholessky decomposition

	const char *funcname = "Ich2";

// Local variables

	int irow, j, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzgmax, nzjgmx, nz, nzju, nzjr, nzjrl, icheck, ibegfr;
	int kii, kjj;
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

	int lwork;
	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	double *fmask, *dpiv, *diagg, *scl, *sclinv;
	double *lelem, *delem;
	double *aloc, *vloc, *eig, *work;

// Allocate and init working arrays including the scaling

// Integer arrays

	Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	Ich2AllocScale (fmask, dpiv, diagg, scl, sclinv,
					lelem, delem,
					aloc, vloc, eig, lwork, work,
					ops);

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

	for (irow=0;irow<nsupr;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			cout << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;

// Output hystogram of the scaling values

			int itype = 2;
			int nloc = irow*bla;

			csvhyst (cout, "Ich2Dia", nloc, dpiv);
//			csvhyst (_fout,"Ich2Dia", nloc, dpiv);

		};

// Init current scaled row

		Ich2InitRow (irow, nlist1, lstloc, icycle, imask, 
					diagg, fmask, scl, 
					aloc, vloc,
					iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];

			for (kii=0;kii<bla;kii++) {
				for (kjj=0;kjj<bla;kjj++) {
					lelem[kjj*bla+kii] = g[j*bla_2+kjj*bla+kii];
				};
			};

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

		Ich2PivScaleRow (irow, delem, nlist1, lstupd, fmask, dpiv, 
							aloc, eig, lwork, work,
							iops);

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx) {

// Free the memory containing R part that is not in use any more

			Ich2FreeR (irow, 
						ibegfr, nzjg, nzju, nzjrl, 
						iv, ig, jg, g);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx) {
			cout << " Insufficient workspace in Ich2, Irow = " << irow << endl;
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

	int itype = 2;
	int nloc = n;

	csvhyst (cout, "Ich2Dia", nloc, dpiv);
//	csvhyst (_fout,"Ich2Dia", nloc, dpiv);

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

	CSMatrixRB temp (nsupr,nzju,bla);

	Ich2StoreScaleU (scl, sclinv, ig, jg, g, 
					temp, aloc,
					ops);

// Free factorization arrays 

	delete [] scl;
	delete [] sclinv;
	delete [] lelem;
	delete [] delem;
	delete [] aloc;
	delete [] vloc;
	delete [] eig;
	delete [] work;

	delete [] ig;
	delete [] jg;
	delete [] addrg;
	delete [] g;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Ich2 factorization statistics

	nz = ia[nsupr];

	densu = 1.0e2 * (double) nzju / (double) nz;
	densr = 1.0e2 * (double) nzjr / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	nz = (2*ia[nsupr]-nsupr) * bla_2;

	ops = ops / (double) nz;

	cout << " Ich2 preconditioner generation statistics: " << endl;
	cout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	cout << "     Costs = " << ops << " MvmA flops. " << endl;
	cout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Return the result

	return temp;

};
