//------------------------------------------------------------------------------------------------
// File: ich2rs.cpp
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

#include "smatrix.h"
#include "slvparam.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CSMatrixRS: Allocate double data and init scaling for Ich2
//========================================================================================
void CSMatrixRS::Ich2AllocScale (ofstream &_fout, double *&_fmask, double *&_dpiv, // Allocate double data and init scaling for Ich2
								int *&_bsdia, double *&_diagg, double *&_scl, double *&_sclinv,
								double *&_lelem, double *&_delem,
								double *&_aloc, double *&_vloc, 
								double *&_eig, int &_lwork, double *&_work,
								double &_ops) const {

	const char *funcname = "Ich2AllocScale";

	_fmask = new double [n*blamx];
	if (!_fmask) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_bsdia = new int [nsupr];
	if (!_bsdia) MemoryFail (funcname);
	_diagg = new double [n*blamx];
	if (!_diagg) MemoryFail (funcname);
	_scl = new double [n*blamx];
	if (!_scl) MemoryFail (funcname);
	_sclinv = new double [n*blamx];
	if (!_sclinv) MemoryFail (funcname);
	_lelem = new double [blamx*blamx];
	if (!_lelem) MemoryFail (funcname);
	_delem = new double [blamx*blamx];
	if (!_delem) MemoryFail (funcname);

// Allocate working arrays

	int info=0;

	_lwork = 10 * blamx;

	_aloc = new double [blamx*blamx];
	if (!_aloc) MemoryFail (funcname);
	_vloc = new double [blamx*blamx];
	if (!_vloc) MemoryFail (funcname);
	_eig = new double [blamx];
	if (!_eig) MemoryFail (funcname);
	_work = new double [_lwork];
	if (!_work) MemoryFail (funcname);

// Compute array bsdia

	int isup;
	int j, kii, kjj, kkk, ibssc;
	int blai, ibs, ibsv, ibsd, blai_2, kjjb;

	int nz=0;

	for (isup=0;isup<nsupr;isup++) {
		_bsdia[isup] = nz;
		blai = sprndr[isup+1]-sprndr[isup];
		nz += blai * blai;
	};

// Compute scaling and its inverse

	double *paloc, *pa;
	double *peig, *pscl;
	double *pfmask;

	double aux;

	for (isup=0;isup<nsupr;isup++) {

		blai = sprndr[isup+1]-sprndr[isup];
		blai_2 = blai*blai;
		ibsv = sprndr[isup];
		ibsd = _bsdia[isup];

// Assign diagonal supernode

		j = ia[isup];
		ibs = bsa[j];

		paloc = _aloc;
		pa = a+ibs;
//		for (kii=0;kii<blai_2;kii++) _aloc[kii] = a[ibs+kii];
		for (kii=0;kii<blai_2;kii++) *paloc++ = *pa++;

// Compute its symmetric spectral decomposition

		int blal = blai;

		dsyev_ ("V", "U", &blal, _aloc, &blal, _eig, _work, &_lwork, &info);

		if (info != 0) throw " Error in the Lapack routine DSYEV";

		_ops += 8 * blai_2*blai;

// Check for negative pivots

		for (kii=0;kii<blai;kii++) {
			if (_eig[kii] < 0.0e0) {
				cout << " Negative eigenvalue when scaling Isup = " << isup << endl;
				throw " Negative eigenvalue when scaling ";
			};
		};

// Store scaling factors

		for (kii=0;kii<blai;kii++) {
			aux = _eig[kii];
			aux = sqrt(aux);
			_eig[kii] = aux;
			_dpiv[ibsv+kii] = _eig[kii];
		};

		for (kii=0;kii<blai;kii++) {
			peig = _eig;
			paloc = _aloc+kii;
			pscl = _sclinv+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinv[ibsd+kjj*blai+kii] = _aloc[kjj*blai+kii] * _eig[kjj];
				*pscl = *paloc * *peig++;
				paloc += blai;
				pscl += blai;
			};
		};

		for (kii=0;kii<blai;kii++) {
			aux = _eig[kii];
			aux = 1.0e0 / aux;
			_eig[kii] = aux;
		};

		peig = _eig;
		for (kii=0;kii<blai;kii++) {
			paloc = _aloc+kii*blai;
			pscl = _scl+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_scl[ibsd+kjj*blai+kii] = _aloc[kii*blai+kjj] * _eig[kii];
				*pscl = *paloc++ * *peig;
				pscl += blai;
			};
			peig++;
		};

		_ops += 2*blai_2 + 2*blai;

	};

// Reassign scaling factors if necessary

	if (false) {

		for (kii=0;kii<n*blamx;kii++) _scl[kii] = 0.0e0;
		for (kii=0;kii<n*blamx;kii++) _sclinv[kii] = 0.0e0;

		for (isup=0;isup<nsupr;isup++) {
			blai = sprndr[isup+1]-sprndr[isup];
			ibsd = _bsdia[isup];
			for (kii=0;kii<blai;kii++) {
				_scl[ibsd+kii*blai+kii] = 1.0e0;
				_sclinv[ibsd+kii*blai+kii] = 1.0e0;
			};
		};

	};

// Init working arrays

	pfmask = _fmask;
	for (kii=0;kii<n*blamx;kii++) *pfmask++ = 0.0e0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		blai_2 = blai*blai;
		j = ia[isup];
		ibs = bsa[j];
		ibssc = _bsdia[isup];
		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			paloc = _aloc+kii;
			for (kjj=0;kjj<blai;kjj++) {
				aux = 0.0e0;
				pa = a+ibs+kjjb;
				pscl = _scl+ibssc+kii;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scl[ibssc+kkk*blai+kii];
					aux += *pa++ * *pscl;
					pscl += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*paloc = aux;
				paloc += blai;
				kjjb += blai;
			};
		};
		for (kii=0;kii<blai;kii++) {
			pa = _vloc+kii;
			for (kjj=0;kjj<blai;kjj++) {
				aux = 0.0e0;
				paloc = _aloc+kii;
				pscl = _scl+ibssc+kjj;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _scl[ibssc+kkk*blai+kjj];
					aux += *paloc * *pscl;
					paloc += blai;
					pscl += blai;
				};
//				_vloc[kjj*blai+kii] = aux;
				*pa = aux;
				pa += blai;
			};
		};
		pa = _diagg+ibssc;
		paloc = _vloc;
//		for (kii=0;kii<blai*blai;kii++) _diagg[ibssc+kii] = _vloc[kii];
		for (kii=0;kii<blai_2;kii++) *pa++ = *paloc++;
		_ops += 2 * blai_2 * blai;
	};

// Output hystogram of the scaling values

//	int itype = 1;
	int nloc = n;

	csvhyst (cout, "ScDia",nloc, _dpiv);
	csvhyst (_fout,"ScDia",nloc, _dpiv);

};

// Author: Kharchenko S.A.
// CSMatrixRS: Allocate G data for Ich2
//========================================================================================
void CSMatrixRS::Ich2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ich2
							int *&_ig, int *&_jg, int *&_addrg, double *&_g) const {

	const char *funcname = "Ich2AllocG";

	int nz = ia[nsupr];

	_nzjgmx = int (_memory * double(nz));
	_nzgmax = int (_memory * double(nza));

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int [_nzjgmx];
	if (!_jg) MemoryFail (funcname);
	_addrg = new int [_nzjgmx];
	if (!_addrg) MemoryFail (funcname);
	_g = new double [_nzgmax];
	if (!_g) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixRS: Init current scaled row
//========================================================================================
void CSMatrixRS::Ich2InitRow (int _irow, int &_nlist1, // Init current scaled row
							int *_lstloc, int _icycle, int *_imask, 
							int *_bsdia, double *_diagg, double *_fmask, double *_scl, 
							double *_aloc, double *_vloc,
							int &_iops) const {

//	const char *funcname = "Ich2InitRow";

	int kii, kjj, kkk, ibs, ibssc, jbssc;
	int kjjb, blai, blaj, ibsd, ibsv, jbsv, blaij, blai_2, blaj_2;
	double aux;

	double *paloc;
	double *pscloc;
	double *pfmask;
	double *pvloc;

	blai = sprndr[_irow+1]-sprndr[_irow];
	blai_2 = blai*blai;
	ibsd = _bsdia[_irow];
	ibsv = sprndr[_irow];

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	pfmask = _fmask+ibsv*blai;
	pvloc = _diagg+ibsd;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	for (kii=0;kii<blai_2;kii++) *pfmask++ = *pvloc++;
	_nlist1++;

	for (int j=ia[_irow]+1;j<ia[_irow+1];j++) {
		int jj = ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;

// Local scaling

		ibs = bsa[j];
		ibssc = _bsdia[_irow];
		jbssc = _bsdia[jj];
		blaj = sprndr[jj+1]-sprndr[jj];
		jbsv = sprndr[jj];
		blaij = blai * blaj;
		blaj_2 = blaj*blaj;
		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				paloc = a+ibs+kjjb;
				pscloc = _scl+ibssc+kii;
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scl[ibssc+kkk*blai+kii];
					aux += *paloc++ * *pscloc;
					pscloc += blai;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				aux = 0.0e0;
				paloc = _aloc+kii;
				pscloc = _scl+jbssc+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _scl[jbssc+kkk*blaj+kjj];
					aux += *paloc * *pscloc;
					pscloc += blaj;
					paloc += blai;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		pfmask = _fmask+jbsv*blai;
		pvloc = _vloc;
//		for (kii=0;kii<blai*blaj;kii++) _fmask[jbsv*blai+kii] = _vloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmask++ = *pvloc++;

		_nlist1++;

		_iops += blai*(blaij+blaj_2);

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Init mask and update arrays
//========================================================================================
void CSMatrixRS::Ich2InitMask (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmask) const {

//	const char *funcname = "Ich2InitMask";

	int kii;
	int blai, blaj, jbsv, blaij;
	double *pfmask;

	_nupd = 0;

	blai = sprndr[_irow+1]-sprndr[_irow];

	for (int k=_iv[_irwprv];k<_ig[_irwprv+1];k++) {

		int jcolmn = _jg[k];

		if (_elmisr && jcolmn < 0) goto label1;

		if (jcolmn < 0) jcolmn = -jcolmn;
		if (_imask[jcolmn] != _icycle) {
			_imask[jcolmn] = _icycle;
			_lstloc[_nlist1] = jcolmn;
			blaj = sprndr[jcolmn+1]-sprndr[jcolmn];
			jbsv = sprndr[jcolmn];
			blaij = blai*blaj;
			pfmask = _fmask+jbsv*blai;
//			for (kii=0;kii<blai*blaj;kii++) _fmask[jbsv*blai+kii] = 0.0e0;
			for (kii=0;kii<blaij;kii++) *pfmask++ = 0.0e0;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = k;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Update current row
//========================================================================================
void CSMatrixRS::Ich2UpdateRow (int _irow, int _irwprv, int _nupd, int *_lstupd, int *_indlst, // Update current row
							double *_lelem, double *_fmask, int *_addrg, double *_g, int &_iops) const {

//	const char *funcname = "Ich2UpdateRow";

	int kii, kjj, ibs, jbs;
	int kikb, kjkb;
	int blai, blaj, blak, jbsv, blaik;

	double *pfmloc;
	double *plloc;
	double *pgloc;

	blai = sprndr[_irow+1]-sprndr[_irow];
	blak = sprndr[_irwprv+1]-sprndr[_irwprv];
	blaik = blai*blak;

	for (int k=0;k<_nupd;k++) {
		int icol = _lstupd[k];
		int indg = _indlst[k];
		blaj = sprndr[icol+1]-sprndr[icol];
		jbsv = sprndr[icol];
		ibs = jbsv*blai;
		jbs = _addrg[indg];
		kikb = 0;
		for (kii=0;kii<blai;kii++) {
			pfmloc = _fmask+ibs+kii;
			kjkb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				plloc = _lelem+kikb;
				pgloc = _g+jbs+kjkb;
				double *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
					*pfmloc -= *plloc * *pgloc;
				};
				pfmloc += blai;
				kjkb += blak;
			};
			kikb += blak;
		};

		_iops += blaik * blaj;

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Perform filtering of the current row
//========================================================================================
void CSMatrixRS::Ich2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
								int &_nlist1, int *_lstloc, int *_lstupd,
								double *_fmask, int *_bsdia, double *_diagg, int &_iops) const {

//	const char *funcname = "Ich2FiltrRow";

	int nlist2 = 0;
	int kii, kjj, blai, blaj, ibsv, jbsv, jbsd, blaij;
	double *pfmask, *pfmask1, *pdiag;

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsv = sprndr[_irow];

	for (int i=1;i<_nlist1;i++) {

		int j = _lstloc[i];

		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		jbsd = _bsdia[j];
		blaij = blai*blaj;

		double aux = 0.0e0;

		pfmask = _fmask+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmask[jbsv*blai+kii];
			double auxloc = *pfmask++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};

		if (aux >= _tau2) {
			_lstupd[nlist2] = j;
			nlist2++;
		} else {
			pfmask1 = _fmask+ibsv*blai;
			for (kii=0;kii<blai;kii++) {
				pfmask = _fmask+jbsv*blai+kii;
				pdiag = _diagg+jbsd;
				for (kjj=0;kjj<blaj;kjj++) {
//					double auxloc = _fmask[jbsv*blai+kjj*blai+kii];
					double auxloc = *pfmask;
					if (auxloc < 0.0e0) auxloc = -auxloc;
					auxloc *= _theta;
//					_fmask[ibsv*blai+kii*blai+kii] += auxloc;
					*pfmask1 += auxloc;
//					_diagg[jbsd+kjj*blaj+kjj] += auxloc;
					*pdiag += auxloc;
					pfmask += blai;
					pdiag += blaj+1;
				};
				pfmask1 += blai+1;
			};
		};

	};

	_nlist1 = nlist2;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Compute diagonal pivot and scale current row
//========================================================================================
void CSMatrixRS::Ich2PivScaleRow (int _irow, double *_delem, // Compute diagonal pivot and scale current row
								int _nlist1, int *_lstupd, double *_fmask, double *_dpiv, 
								double *_aloc, double *_eig, int _lwork, double *_work,
								int &_iops) const {

//	const char *funcname = "Ich2PivScaleRow";

	int kii, kjj, kkk, ibs;
	int blai, blaj, ibsv, jbsv;
	double aux;

// Copy diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsv = sprndr[_irow];
	ibs = ibsv * blai;

	for (kii=0;kii<blai*blai;kii++) _aloc[kii] = _fmask[ibs+kii];

// Compute spectral decomposition of the diagonal supernode

	int blal = blai;
	int info = 0;

	dsyev_ ("V", "U", &blal, _aloc, &blal, _eig, _work, &_lwork, &info);

	if (info != 0) throw " Error in the Lapack routine DSYEV";

	_iops += 8 * blai * blai * blai;

// Sort the eigenvalues/eigenvectors

	for (kii=0;kii<blai;kii++) {
		for (kjj=0;kjj<blai;kjj++) _delem[(blai-kii-1)*blai+kjj] = _aloc[kii*blai+kjj];
	};
	for (kii=0;kii<blai*blai;kii++) _aloc[kii] = _delem[kii];

	for (kii=0;kii<blai;kii++) _work[blai-kii-1] = _eig[kii];
	for (kii=0;kii<blai;kii++) _eig[kii] = _work[kii];

// Check for negative pivots

	for (kii=0;kii<blai;kii++) {
		if (_eig[kii] < 0.0e0) throw " Negative pivot in ICH2 ";
		_eig[kii] = sqrt(_eig[kii]);
		_dpiv[ibsv+kii] = _eig[kii];
		_eig[kii] = 1.0e0/_eig[kii];
	};

// Compute diagonal supernode

	int kiib;
	double *pdeloc;
	double *pfmloc;
	double *palloc;

	pdeloc = _delem;
	
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = _delem+kii;
		palloc = _aloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_delem[kjj*blai+kii] = _aloc[kii*blai+kjj] * _eig[kii];
			*pdeloc = *palloc++ * _eig[kii];
			pdeloc += blai;
		};
		kiib += blai;
	};

	_iops += blai * blai;

// Perform scaling of the data in the list

	int kjjb, blaij;

	for (int i=0;i<_nlist1;i++) {

		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		blaij = blai*blaj;

		ibs = jbsv*blai;

		for (kii=0;kii<blai;kii++) {
			palloc = _aloc+kii;
			for (kjj=0;kjj<blaj;kjj++) {
				kjjb = kjj*blai;
				pfmloc = _fmask+ibs+kjjb;
				pdeloc = _delem+kii;
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _delem[kkk*blai+kii] * _fmask[ibs+kjj*blai+kkk];
					aux += *pdeloc * *pfmloc++;
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = aux;
				palloc += blai;
			};
		};

		pfmloc = _fmask+ibs;
		palloc = _aloc;
//		for (kii=0;kii<blaij;kii++) _fmask[ibs+kii] = _aloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmloc++ = *palloc++;

		_iops += blai * blaij;

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Free R part of the data
//========================================================================================
void CSMatrixRS::Ich2FreeR (ofstream &_fout, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, int &_nznew, int *_addrg, double *_g) const {

//	const char *funcname = "Ich2FreeR";

	int iendfr = _irow;

	cout  << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;
	_fout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;

// Take initial address of the first block to be compressed

	int i, j;

	if (_ibegfr != 0) {

		i = _ibegfr;
		j = _ig[i];

		_nznew = _addrg[j];

	} else {
		_nznew = 0;
	};

// Init data for the main cycle

	int ibgfrn = _ibegfr;

	_nzjg = _ig[_ibegfr];

	int ibeg = _ig[_ibegfr];

	int nzjro = _nzjrl;
	_nzjrl = 0;

// Main cycle over the supernode rows

	int blai, blaj, blaij;
	double *pgold, *pgnew;

	for (i=_ibegfr;i<iendfr;i++) {

		blai = sprndr[i+1]-sprndr[i];

// Scan current row

		for (j=ibeg;j<_ig[i+1];j++) {

			int icolmn = _jg[j];

			if (icolmn < 0) {
				icolmn = -icolmn;
				if (icolmn < iendfr) goto label2;
				_nzjrl++;
			};

			blaj = sprndr[icolmn+1]-sprndr[icolmn];
			blaij = blai*blaj;

			if (_nzjg != j) {
				_jg[_nzjg] = _jg[j];
				int ibs = _addrg[j];
				pgold = _g+ibs;
				pgnew = _g+_nznew;
//				for (int kii=0;kii<blaij;kii++) _g [_nznew+kii] = _g [ibs+kii];
				for (int kii=0;kii<blaij;kii++) *pgnew++ = *pgold++;
			};
			_addrg[_nzjg] = _nznew;
			_nznew += blaij;
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

	cout  << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;
	_fout << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Store current row
//========================================================================================
void CSMatrixRS::Ich2StoreRow (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
								double *_delem, int _nlist1, int *_lstupd, double *_fmask, 
								int *_ig, int *_jg, int &_nznew, int *_addrg, double *_g) const {

//	const char *funcname = "Ich2StoreRow";

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	double *pg, *pfmask, *pd;

// Store current diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	_jg[_nzjg] = _irow;
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = _delem+kiib;
		pg = _g+_nznew+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_g [_nznew+kjj*blai+kii] = _delem[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};
	_addrg[_nzjg] = _nznew;
	_nzjg++;
	_nzju++;
	_nznew += blai * blai;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		blaij = blai*blaj;
		double aux = 0.0e0;
		pfmask = _fmask+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmask[jbsv*blai+kii];
			double auxloc = *pfmask++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		_addrg[_nzjg] = _nznew;
		if (aux >= _tau1) {
			_jg[_nzjg] = j;
			pg = _g+_nznew;
			pfmask = _fmask+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _g [_nznew+kii] = _fmask[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
			pg = _g+_nznew;
			pfmask = _fmask+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _g [_nznew+kii] = _fmask[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			_nzjg++;
			_nzjr++;
			_nzjrl++;
			_nzrtot += blaij;
		};
		_nznew += blaij;
	};

	_ig[_irow+1] = _nzjg;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Store the resulting U and scale it
//========================================================================================
void CSMatrixRS::Ich2StoreScaleU (int *_bsdia, double *_scl, double *_sclinv, // Store the resulting U and scale it
								int *_ig, int *_jg, int *_addrg, double *_g,
								CSMatrixRS &_temp, 
								double *_aloc,
								double &_ops) const {

//	const char *funcname = "Ich2StoreScaleU";

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj, blaij;
	double aux;

	double *pscloc;
	double *pgloc;
	double *paloc;

	_temp.m = n;
	_temp.n = n;
	_temp.nsupr = nsupr;
	_temp.nsupc = nsupc;
	_temp.nlist = nsupr;

	int i;

	for (i=0;i<nsupr;i++) _temp.list[i] = i;
	for (i=0;i<=nsupr;i++) {
		_temp.sprndr[i] = sprndr[i];
		_temp.sprndc[i] = sprndc[i];
	};

	_temp.blamx = 0;
	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1] - sprndr[i];
		if (blai > _temp.blamx) _temp.blamx = blai;
	};

	_temp.ia[0] = 0;

	int nzj = 0;
	int nz = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		ibssc = _bsdia[i];
		for (int j=_ig[i];j<_ig[i+1];j++) {
			jbs = _addrg[j];
			int jj = _jg[j];
			int jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = sprndr[jjloc+1]-sprndr[jjloc];
			blaij = blai*blaj;
			jjbssc = _bsdia[jjloc];
			if (jj >= 0) {
				_temp.ja[nzj] = jj;
				_temp.bsa[nzj] = nz;
				if (jj == i) {
					kiib = 0;
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = _aloc+kii;
						for (kjj=0;kjj<blai;kjj++) {
							pgloc = _g+jbs+kjjb;
							pscloc = _scl+ibssc+kiib;
							aux = 0.0e0;
							for (kkk=0;kkk<blai;kkk++) {
//								aux += _scl[ibssc+kii*blai+kkk] * _g[jbs+kjj*blai+kkk];
								aux += *pscloc * *pgloc;
								pgloc++;
								pscloc++;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = aux;
							paloc += blai;
							kjjb += blai;
						};
						kiib += blai;
					};
					paloc = _aloc;
					pgloc = _temp.a+nz;
//					for (kii=0;kii<blaij;kii++) _temp.a[nz+kii] = _aloc[kii];
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;
				} else {
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = _aloc+kii;
						for (kjj=0;kjj<blaj;kjj++) {
							aux = 0.0e0;
							pgloc = _g+jbs+kii;
							pscloc = _sclinv+jjbssc+kjj;
							for (kkk=0;kkk<blaj;kkk++) {
//								aux += _g[jbs+kkk*blai+kii] * _sclinv[jjbs+kkk*blaj+kjj];
								aux += *pgloc * *pscloc;
								pgloc += blai;
								pscloc += blaj;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = aux;
							paloc += blai;
							kjjb += blai;
						};
					};
					paloc = _aloc;
					pgloc = _temp.a+nz;
//					for (kii=0;kii<blaij;kii++) _temp.a[nz+kii] = _aloc[kii];
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;
				};
				nzj++;
				nz += blaij;
 				_ops += blaij*(blai+blaj);
			};
		};
		_temp.ia[i+1] = nzj;
	};

	_temp.nzja = nzj;
	_temp.nza = nz;

};


// Author: Kharchenko S.A.
// CSMatrixRS: Store the resulting U and scale it
//========================================================================================
void CSMatrixRS::Ich2StoreScaleUDynamic (int *_bsdia, double *_scl, double *_sclinv, // Store the resulting U and scale it
								int *_ig, int **_jg, int **_addrg, double **_g,
								CSMatrixRS &_temp, 
								double *_aloc,
								double &_ops) const {

//	const char *funcname = "Ich2StoreScaleU";

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj, blaij;
	double aux;
	int *jgloc, *addrgloc;
	double *gloc;

	double *pscloc;
	double *pgloc;
	double *paloc;

	_temp.m = n;
	_temp.n = n;
	_temp.nsupr = nsupr;
	_temp.nsupc = nsupc;
	_temp.nlist = nsupr;

	int i;

	for (i=0;i<nsupr;i++) _temp.list[i] = i;
	for (i=0;i<=nsupr;i++) {
		_temp.sprndr[i] = sprndr[i];
		_temp.sprndc[i] = sprndc[i];
	};

	_temp.blamx = 0;
	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1] - sprndr[i];
		if (blai > _temp.blamx) _temp.blamx = blai;
	};

	_temp.ia[0] = 0;

	int nzj = 0;
	int nz = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		ibssc = _bsdia[i];
		jgloc = _jg[i];
		addrgloc = _addrg[i];
		gloc = _g[i];
		for (int j=_ig[i];j<_ig[i+1];j++) {
			int jloc = j-_ig[i];
			jbs = addrgloc[jloc];
			int jj = jgloc[jloc];
			int jjloc = jj;
			if (jjloc < 0) jjloc = -jjloc;
			blaj = sprndr[jjloc+1]-sprndr[jjloc];
			blaij = blai*blaj;
			jjbssc = _bsdia[jjloc];
			if (jj >= 0) {
				_temp.ja[nzj] = jj;
				_temp.bsa[nzj] = nz;
				if (jj == i) {
					kiib = 0;
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = _aloc+kii;
						for (kjj=0;kjj<blai;kjj++) {
							pgloc = gloc+jbs+kjjb;
							pscloc = _scl+ibssc+kiib;
							aux = 0.0e0;
							for (kkk=0;kkk<blai;kkk++) {
//								aux += _scl[ibssc+kii*blai+kkk] * _g[jbs+kjj*blai+kkk];
								aux += *pscloc * *pgloc;
								pgloc++;
								pscloc++;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = aux;
							paloc += blai;
							kjjb += blai;
						};
						kiib += blai;
					};
					paloc = _aloc;
					pgloc = _temp.a+nz;
//					for (kii=0;kii<blaij;kii++) _temp.a[nz+kii] = _aloc[kii];
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;
				} else {
					for (kii=0;kii<blai;kii++) {
						kjjb = 0;
						paloc = _aloc+kii;
						for (kjj=0;kjj<blaj;kjj++) {
							aux = 0.0e0;
							pgloc = gloc+jbs+kii;
							pscloc = _sclinv+jjbssc+kjj;
							for (kkk=0;kkk<blaj;kkk++) {
//								aux += _g[jbs+kkk*blai+kii] * _sclinv[jjbs+kkk*blaj+kjj];
								aux += *pgloc * *pscloc;
								pgloc += blai;
								pscloc += blaj;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = aux;
							paloc += blai;
							kjjb += blai;
						};
					};
					paloc = _aloc;
					pgloc = _temp.a+nz;
//					for (kii=0;kii<blaij;kii++) _temp.a[nz+kii] = _aloc[kii];
					for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;
				};
				nzj++;
				nz += blaij;
				_ops += blaij*(blai+blaj);
			};
		};
		_temp.ia[i+1] = nzj;
	};

	_temp.nzja = nzj;
	_temp.nza = nz;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Second order Incomplete Cholessky decomposition
//========================================================================================
CSMatrixRS CSMatrixRS::Ich2 (ofstream &_fout, CSlvParam &_param) const { // Second order Incomplete Cholessky decomposition

//	const char *funcname = "Ich2";

// Local variables

	int irow, i, j, jj, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzgmax, nzjgmx;
	int nz, nznew, nzju, nzjr, nzjrl, nzrtot, icheck, ibegfr;
	int kii, blai, blaj, ibs;
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
	int *bsdia;
	double *fmask, *dpiv, *diagg, *scl, *sclinv;
	double *lelem, *delem;
	double *aloc, *vloc, *eig, *work;

// Allocate and init working arrays including the scaling

// Integer arrays

	Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	Ich2AllocScale (_fout, fmask, dpiv, bsdia, diagg, 
					scl, sclinv, lelem, delem, 
					aloc, vloc, eig, lwork, work, ops);

// Allocate factorization arrays

	Ich2AllocG (_param.memory, nzjgmx, nzgmax, ig, jg, addrg, g);

// Main cycle over the rows

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;
	nznew = 0;
	nzrtot = 0;

	ig[0] = 0;

	for (irow=0;irow<nsupr;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			cout  << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			_fout << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;

// Output hystogram of the scaling values

//			int itype = 2;
			int nloc = sprndr[irow];

			csvhyst (cout, "Ich2Dia", nloc, dpiv);
			csvhyst (_fout,"Ich2Dia", nloc, dpiv);

		};

// Init current scaled row

		Ich2InitRow (irow, nlist1, lstloc, icycle, imask, 
						bsdia, diagg, fmask, scl, aloc, vloc, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];
			ibs = addrg[j];

			blai = sprndr[irwprv+1]-sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = sprndr[icolmnl+1]-sprndr[icolmnl];

			for (kii=0;kii<blai*blaj;kii++) lelem[kii] = g[ibs+kii];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMask (elmisr, irow, irwprv, 
							iv, ig, jg, 
							nupd, lstupd, indlst, nlist1, lstloc, 
							icycle, imask, fmask);

// Perform updates

			Ich2UpdateRow (irow, irwprv, nupd, lstupd, indlst, 
							lelem, fmask, addrg, g, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTransp (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow (irow, _param.theta, _param.tau2, 
						nlist1, lstloc, lstupd, 
						fmask, bsdia, diagg, iops);

// Compute and store the pivot value and scale current row

		Ich2PivScaleRow (irow, delem, nlist1, lstupd, 
							fmask, dpiv, 
							aloc, eig, lwork, work, iops);

// Count the local number of elements to be stored

		blai = sprndr[irow+1]-sprndr[irow];

		int nzloc=blai*blai;

		for (j=0;j<nlist1;j++) {
			jj = lstupd[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			nzloc += blai*blaj;
		};

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx || nznew+nzloc > nzgmax) {

// Free the memory containing R part that is not in use any more

			Ich2FreeR (_fout, irow, ibegfr, nzjg, nzju, nzjrl, 
						iv, ig, jg, 
						nznew, addrg, g);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx || nznew+nzloc > nzgmax) {
			cout  << " Insufficient workspace in Ich2, Irow = " << irow << endl;
			_fout << " Insufficient workspace in Ich2, Irow = " << irow << endl;
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

		Ich2StoreRow (irow, _param.tau1, nzjg, nzju, nzjr, nzjrl, nzrtot,
						delem, nlist1, lstupd, fmask, 
						ig, jg, nznew, addrg, g);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = n;

	csvhyst (cout, "Ich2Dia", nloc, dpiv);
	csvhyst (_fout,"Ich2Dia", nloc, dpiv);

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

// Compute the size of the resulting matrix

	int nztot = 0;
	nzju = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=ig[i];j<ig[i+1];j++) {
			jj = jg[j];
			if (jj >= 0) {
				blaj = sprndr[jj+1]-sprndr[jj];
				nzju++;
				nztot += blai*blaj;
			};
		};
	};

// Perform the final filtering and scaling and compute matrix on return

	CSMatrixRS temp (nsupr,nzju,nztot);

	Ich2StoreScaleU (bsdia, scl, sclinv, 
						ig, jg, addrg, g, 
						temp, 
						aloc, ops);

// Free factorization arrays 

	blai = sprndr[nsupr]-sprndr[nsupr-1];

	int nzmema = 2*nza-bsdia[nsupr-1]-blai*blai;

	delete [] bsdia;
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

	nz = nza;

	densu = 1.0e2 * (double) nztot / (double) nz;
	densr = 1.0e2 * (double) nzrtot / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	cout  << " Ich2 preconditioner generation statistics: " << endl;
	cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Ich2 preconditioner generation statistics: " << endl;
	_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Return the result

	return temp;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Allocate double data and init scaling for Ich2
//========================================================================================
void CSMatrixRS::Ich2AllocScale (ofstream &_fout, // Allocate double data and init scaling for Ich2
									int _nblks, int *_sp2blk, FILE **_mtrafiles,
									double *&_gread, int &_nzgread,
									double *&_fmask, double *&_dpiv, 
									int *&_bsdia, double *&_diagg, 
									int *&_bsblk, FILE **_mtrlufiles,
									int *&_addrscl, int *&_addrsclinv,
									double *&_lelem, double *&_delem,
									double *&_aloc, double *&_vloc, 
									double *&_eig, int &_lwork, double *&_work,
									double &_ops) const {

	const char *funcname = "Ich2AllocScale";

// Allocate working arrays

	_bsdia = new int [nsupr];
	if (!_bsdia) MemoryFail (funcname);
	_addrscl = new int [nsupr];
	if (!_addrscl) MemoryFail (funcname);
	_addrsclinv = new int [nsupr];
	if (!_addrsclinv) MemoryFail (funcname);
	_bsblk = new int [_nblks];
	if (!_bsblk) MemoryFail (funcname);

	_fmask = new double [n*blamx];
	if (!_fmask) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_lelem = new double [blamx*blamx];
	if (!_lelem) MemoryFail (funcname);
	_delem = new double [blamx*blamx];
	if (!_delem) MemoryFail (funcname);

	int info=0;

	_lwork = 10 * blamx;

	_aloc = new double [blamx*blamx];
	if (!_aloc) MemoryFail (funcname);
	_vloc = new double [blamx*blamx];
	if (!_vloc) MemoryFail (funcname);
	_eig = new double [blamx];
	if (!_eig) MemoryFail (funcname);
	_work = new double [_lwork];
	if (!_work) MemoryFail (funcname);

// Compute array bsdia

	int isup;
	int j, kii, kjj, kkk, ibssc;
	int blai, ibs, ibsv, ibsd;
	int iblk;

	int nz=0;

	for (isup=0;isup<nsupr;isup++) {
		_bsdia[isup] = nz;
		blai = sprndr[isup+1]-sprndr[isup];
		nz += blai * blai;
	};

	_diagg = new double [nz];
	if (!_diagg) MemoryFail (funcname);

// Determine the overestimate of the maximal work space required to read g data

	_nzgread = blamx * blamx;

	_gread = new double [_nzgread];
	if (!_gread) MemoryFail (funcname);

// Compute scaling and its inverse and store it to disk

	double aux;

	for (iblk=0;iblk<_nblks;iblk++) {
		_bsblk[iblk] = 0;
	};

	for (isup=0;isup<nsupr;isup++) {

		blai = sprndr[isup+1]-sprndr[isup];
		ibsv = sprndr[isup];
		ibsd = _bsdia[isup];

// Assign diagonal supernode

		j = ia[isup];
		ibs = bsa[j];
		iblk = _sp2blk[isup];

//		for (kii=0;kii<blai*blai;kii++) _aloc[kii] = a[ibs+kii];

		FGet (_mtrafiles[iblk],blai*blai,_aloc,ibs);

// Compute its symmetric spectral decomposition

		int blal = blai;

//		DSYEVDUM (&blal, _aloc, &blal, _eig, _work, &_lwork, &info);
		dsyev_ ("V", "U", &blal, _aloc, &blal, _eig, _work, &_lwork, &info);

		if (info != 0) throw " Error in the Lapack routine DSYEV";

		_ops += 8 * blai*blai*blai;

// Check for negative pivots

		for (kii=0;kii<blai;kii++) {
			if (_eig[kii] < 0.0e0) {
				cout << " Negative eigenvalue when scaling Isup = " << isup << endl;
				throw " Negative eigenvalue when scaling ";
			};
		};

// Store scaling factors

		for (kii=0;kii<blai;kii++) {
			aux = _eig[kii];
			aux = sqrt(aux);
			_eig[kii] = aux;
			_dpiv[ibsv+kii] = _eig[kii];
		};

		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinv[ibsd+kjj*blai+kii] = _aloc[kjj*blai+kii] * _eig[kjj];
				_delem[kjj*blai+kii] = _aloc[kjj*blai+kii] * _eig[kjj];
			};
		};

		ibs = _bsblk[iblk];

		FPut (_mtrlufiles[iblk],blai*blai,_delem,ibs);

		_addrsclinv[isup] = ibs;

		_bsblk[iblk] += blai*blai;

		for (kii=0;kii<blai;kii++) {
			aux = _eig[kii];
			aux = 1.0e0 / aux;
			_eig[kii] = aux;
		};

		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<blai;kjj++) {
				_delem[kjj*blai+kii] = _aloc[kii*blai+kjj] * _eig[kii];
			};
		};

		ibs = _bsblk[iblk];

		FPut (_mtrlufiles[iblk],blai*blai,_delem,ibs);

		_addrscl[isup] = ibs;

		_bsblk[iblk] += blai*blai;

		_ops += 2*blai*blai + 2*blai;

	};

// Init working arrays

	for (kii=0;kii<n*blamx;kii++) _fmask[kii] = 0.0e0;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		j = ia[isup];
		iblk = _sp2blk[isup];

		ibs = bsa[j];

		FGet (_mtrafiles[iblk],blai*blai,_gread,ibs);

		ibs = _addrscl[isup];

		FGet (_mtrlufiles[iblk],blai*blai,_delem,ibs);

		ibssc = _bsdia[isup];
		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<blai;kjj++) {
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scl[ibssc+kkk*blai+kii];
					aux += _gread[kjj*blai+kkk] * _delem[kkk*blai+kii];
				};
				_aloc[kjj*blai+kii] = aux;
			};
		};
		for (kii=0;kii<blai;kii++) {
			for (kjj=0;kjj<blai;kjj++) {
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
					aux += _aloc[kkk*blai+kii] * _delem[kkk*blai+kjj];
				};
				_vloc[kjj*blai+kii] = aux;
			};
		};

		for (kii=0;kii<blai*blai;kii++) _diagg[ibssc+kii] = _vloc[kii];

		_ops += 2 * blai * blai * blai;

	};

// Output hystogram of the scaling values

//	int itype = 1;
	int nloc = n;

	csvhyst (cout, "ScDia",nloc, _dpiv);
	csvhyst (_fout,"ScDia",nloc, _dpiv);

};

// Author: Kharchenko S.A.
// CSMatrixRS: Allocate G data for Ich2
//========================================================================================
void CSMatrixRS::Ich2AllocG (double _memory, int &_nzjgmx, // Allocate G data for Ich2
							int *&_ig, int *&_jg, int *&_addrg) const {

	const char *funcname = "Ich2AllocG";

	int nz = ia[nsupr];

	_nzjgmx = int (_memory * double(nz));

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int [_nzjgmx];
	if (!_jg) MemoryFail (funcname);
	_addrg = new int [_nzjgmx];
	if (!_addrg) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixRS: Init current scaled row
//========================================================================================
void CSMatrixRS::Ich2InitRow (int _irow, // Init current scaled row
								int *_sp2blk, FILE *_mtrafile, FILE **_mtrlufiles,
								int &_nlist1, 
								int *_lstloc, int _icycle, int *_imask, 
								int *_bsdia, double *_diagg, double *_fmask, int *_addrscl, double *_delem,
								double *_aloc, double *_vloc,
								int &_iops) const {

//	const char *funcname = "Ich2InitRow";

	int kii, kjj, kkk, ibs, ibssc, jbssc;
	int kjjb, blai, blaj, ibsd, ibsv, jbsv;
	int iblk, jblk;
	double aux;

	double *paloc;
	double *pscloc;

	iblk = _sp2blk[_irow];

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsd = _bsdia[_irow];
	ibsv = sprndr[_irow];

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	for (kii=0;kii<blai*blai;kii++) {
		_fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	};
	_nlist1++;

	for (int j=ia[_irow]+1;j<ia[_irow+1];j++) {
		int jj = ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;

// Local scaling

		jblk = _sp2blk[jj];

		ibssc = _bsdia[_irow];
		jbssc = _bsdia[jj];
		blaj = sprndr[jj+1]-sprndr[jj];
		jbsv = sprndr[jj];

		ibs = bsa[j];

		FGet (_mtrafile,blai*blaj,_vloc,ibs);

		ibs = _addrscl[_irow];

		FGet (_mtrlufiles[iblk],blai*blai,_delem,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				paloc = _vloc+kjjb;
				pscloc = _delem+kii;
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scl[ibssc+kkk*blai+kii];
					aux += *paloc++ * *pscloc;
					pscloc += blai;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};

		ibs = _addrscl[jj];

		FGet (_mtrlufiles[jblk],blaj*blaj,_delem,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				aux = 0.0e0;
				paloc = _aloc+kii;
				pscloc = _delem+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _scl[jbssc+kkk*blaj+kjj];
					aux += *paloc * *pscloc;
					pscloc += blaj;
					paloc += blai;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		for (kii=0;kii<blai*blaj;kii++) _fmask[jbsv*blai+kii] = _vloc[kii];

		_nlist1++;

		_iops += blai*(blai*blaj+blaj*blaj);

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Update current row
//========================================================================================
void CSMatrixRS::Ich2UpdateRow (int _irow, // Update current row
								FILE *_mtrlufile, double *_gread,
								int _irwprv, int _nupd, int *_lstupd, int *_indlst, 
								double *_lelem, double *_fmask, int *_addrg, int &_iops) const {

//	const char *funcname = "Ich2UpdateRow";

	int kii, kjj, ibs, jbs;
	int kikb, kjkb;
	int blai, blaj, blak, jbsv;

	double *pfmloc;
	double *plloc;
	double *pgloc;

	blai = sprndr[_irow+1]-sprndr[_irow];
	blak = sprndr[_irwprv+1]-sprndr[_irwprv];

	for (int k=0;k<_nupd;k++) {
		int icol = _lstupd[k];
		int indg = _indlst[k];
		blaj = sprndr[icol+1]-sprndr[icol];
		jbsv = sprndr[icol];
		ibs = jbsv*blai;
		jbs = _addrg[indg];

		FGet (_mtrlufile,blaj*blak,_gread,jbs);

		kikb = 0;
		for (kii=0;kii<blai;kii++) {
			pfmloc = _fmask+ibs+kii;
			kjkb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				plloc = _lelem+kikb;
				pgloc = _gread+kjkb;
				double *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
					*pfmloc -= *plloc * *pgloc;
				};
				pfmloc += blai;
				kjkb += blak;
			};
			kikb += blak;
		};

		_iops += blai * blaj * blak;

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Free R part of the data
//========================================================================================
void CSMatrixRS::Ich2FreeR (int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, int &_nznew, int *_addrg) const {

//	const char *funcname = "Ich2FreeR";

	int iendfr = _irow;

	cout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;

// Take initial address of the first block to be compressed

	int i, j;

	if (_ibegfr != 0) {

		i = _ibegfr;
		j = _ig[i];

		_nznew = 0;

		for (i=0;i<_ibegfr;i++) {
			int blai = sprndr[i+1]-sprndr[i];
			for (j=_ig[i];j<_ig[i+1];j++) {
				int icolmn = _jg[j];
				if (icolmn < 0) icolmn = -icolmn;
				int blaj = sprndr[icolmn+1]-sprndr[icolmn];
				_nznew += blai*blaj;
			};
		};

	} else {
		_nznew = 0;
	};

// Init data for the main cycle

	int ibgfrn = _ibegfr;

	_nzjg = _ig[_ibegfr];

	int ibeg = _ig[_ibegfr];

	int nzjro = _nzjrl;
	_nzjrl = 0;

	int blai, blaj;

// Main cycle over the supernode rows

	for (i=_ibegfr;i<iendfr;i++) {

		blai = sprndr[i+1]-sprndr[i];

// Scan current row

		for (j=ibeg;j<_ig[i+1];j++) {

			int icolmn = _jg[j];

			if (icolmn < 0) {
				icolmn = -icolmn;
				if (icolmn < iendfr) goto label2;
				_nzjrl++;
			};

			blaj = sprndr[icolmn+1]-sprndr[icolmn];

			if (_nzjg != j) {
				_jg[_nzjg] = _jg[j];
//				int ibs = _addrg[j];
//				for (int kii=0;kii<blai*blaj;kii++) _g [_nznew+kii] = _g [ibs+kii];
			};
//			_addrg[_nzjg] = _nznew;
			_addrg[_nzjg] = _addrg[j];
			_nznew += blai*blaj;
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
// CSMatrixRS: Store current row
//========================================================================================
void CSMatrixRS::Ich2StoreRow (int _irow, // Store current row
								FILE *_mtrlufile, double *_aloc,
								double _tau1,
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
								double *_delem, int _nlist1, int *_lstupd, double *_fmask, 
								int *_ig, int *_jg, int &_nznew, int &_nzgblk, int *_addrg) const {

//	const char *funcname = "Ich2StoreRow";

	int kii, kjj;
	int blai, blaj, jbsv;

// Store current diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	_jg[_nzjg] = _irow;
//	for (kii=0;kii<blai;kii++) {
//		for (kjj=0;kjj<blai;kjj++) {
//			_g [_nznew+kjj*blai+kii] = _delem[kii*blai+kjj];
//		};
//	};
	for (kii=0;kii<blai;kii++) {
		for (kjj=0;kjj<blai;kjj++) {
			_aloc[kjj*blai+kii] = _delem[kii*blai+kjj];
		};
	};

	FPut (_mtrlufile,blai*blai,_aloc,_nzgblk);

	_addrg[_nzjg] = _nzgblk;
	_nzjg++;
	_nzju++;
	_nznew += blai * blai;
	_nzgblk += blai * blai;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		double aux = 0.0e0;
		for (kii=0;kii<blai*blaj;kii++) {
			double auxloc = _fmask[jbsv*blai+kii];
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		_addrg[_nzjg] = _nzgblk;
		if (aux >= _tau1) {
			_jg[_nzjg] = j;
//			for (kii=0;kii<blai*blaj;kii++) _g [_nznew+kii] = _fmask[jbsv*blai+kii];
			FPut (_mtrlufile,blai*blaj,_fmask+jbsv*blai,_nzgblk);
			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
//			for (kii=0;kii<blai*blaj;kii++) _g [_nznew+kii] = _fmask[jbsv*blai+kii];
			FPut (_mtrlufile,blai*blaj,_fmask+jbsv*blai,_nzgblk);
			_nzjg++;
			_nzjr++;
			_nzjrl++;
			_nzrtot += blai*blaj;
		};
		_nznew += blai*blaj;
		_nzgblk += blai*blaj;
	};

	_ig[_irow+1] = _nzjg;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Store the resulting U and scale it
//========================================================================================
void CSMatrixRS::Ich2StoreScaleU (int _nblks, int *_blks, FILE **_mtrlufiles, // Store the resulting U and scale it
									int *_sp2blk, int *_bsdia, int *_addrscl, int *_addrsclinv, int *_bsblk,
									int *_ig, int *_jg, int *_addrg,
									CSMatrixRS &_temp, 
									double *_aloc, double *_vloc, double *_delem, 
									double &_ops) const {

//	const char *funcname = "Ich2StoreScaleU";

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj;
	double aux;

	double *pscloc;
	double *pgloc;
	double *paloc;

	_temp.m = n;
	_temp.n = n;
	_temp.nsupr = nsupr;
	_temp.nsupc = nsupc;
	_temp.nlist = nsupr;

	int i;

	for (i=0;i<nsupr;i++) _temp.list[i] = i;

	for (i=0;i<=nsupr;i++) {
		_temp.sprndr[i] = sprndr[i];
		_temp.sprndc[i] = sprndc[i];
	};

	_temp.blamx = 0;
	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1] - sprndr[i];
		if (blai > _temp.blamx) _temp.blamx = blai;
	};

	_temp.ia[0] = 0;

	int nzj = 0;
	int nz = 0;

	for (int iblk=0;iblk<_nblks;iblk++) {

		int nzgblk = _bsblk[iblk];

		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibssc = _bsdia[i];
			for (int j=_ig[i];j<_ig[i+1];j++) {
				jbs = _addrg[j];
				int jj = _jg[j];
				int jjloc = jj;
				if (jjloc < 0) jjloc = -jjloc;
				blaj = sprndr[jjloc+1]-sprndr[jjloc];

				FGet (_mtrlufiles[iblk],blai*blaj,_vloc,jbs);

				jjbssc = _bsdia[jjloc];

				int jblk = _sp2blk[jjloc];

				if (jj >= 0) {
					_temp.ja[nzj] = jj;
					_temp.bsa[nzj] = nzgblk;

					if (jj == i) {

						int ibs = _addrscl[i];
						FGet (_mtrlufiles[iblk],blai*blai,_delem,ibs);

						kiib = 0;
						for (kii=0;kii<blai;kii++) {
							kjjb = 0;
							paloc = _aloc+kii;
							for (kjj=0;kjj<blai;kjj++) {
								pgloc = _vloc+kjjb;
								pscloc = _delem+kiib;
								aux = 0.0e0;
								for (kkk=0;kkk<blai;kkk++) {
//									aux += _scl[ibssc+kii*blai+kkk] * _g[jbs+kjj*blai+kkk];
									aux += *pscloc * *pgloc;
									pgloc++;
									pscloc++;
								};
//								_aloc[kjjb+kii] = aux;
								*paloc = aux;
								paloc += blai;
								kjjb += blai;
							};
							kiib += blai;
						};

//						for (kii=0;kii<blai*blai;kii++) _temp.a[nz+kii] = _aloc[kii];
						FPut (_mtrlufiles[iblk],blai*blai,_aloc,nzgblk);

					} else {

						int ibs = _addrsclinv[jjloc];
						FGet (_mtrlufiles[jblk],blaj*blaj,_delem,ibs);

						for (kii=0;kii<blai;kii++) {
							kjjb = 0;
							paloc = _aloc+kii;
							for (kjj=0;kjj<blaj;kjj++) {
								aux = 0.0e0;
								pgloc = _vloc+kii;
								pscloc = _delem+kjj;
								for (kkk=0;kkk<blaj;kkk++) {
//									aux += _g[jbs+kkk*blai+kii] * _sclinv[jjbs+kkk*blaj+kjj];
									aux += *pgloc * *pscloc;
									pgloc += blai;
									pscloc += blaj;
								};
//								_aloc[kjjb+kii] = aux;
								*paloc = aux;
								paloc += blai;
								kjjb += blai;
							};
						};

//						for (kii=0;kii<blai*blaj;kii++) _temp.a[nz+kii] = _aloc[kii];
						FPut (_mtrlufiles[iblk],blai*blaj,_aloc,nzgblk);

					};
					nzj++;
					nz += blai*blaj;
					nzgblk += blai*blaj;
					_ops += blai*blaj*(blai+blaj);
				};
			};
			_temp.ia[i+1] = nzj;
		};
	};

	_temp.nzja = nzj;
	_temp.nza = 0;
	_temp.nzatot = nz;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Second order Incomplete Cholessky decomposition in the out-of-core
//========================================================================================
CSMatrixRS CSMatrixRS::Ich2 (ofstream &_fout, // Second order Incomplete Cholessky decomposition in the out-of-core
								int _nblks, int *_blks, 
								FILE **_mtrafiles, FILE **_mtrlufiles, 
								CSlvParam &_param) const {

	const char *funcname = "Ich2";

// Local variables

	int irow, i, j, jj, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzjgmx;
	int nz, nznew, nzju, nzjr, nzjrl, nzrtot, icheck, ibegfr;
	int nzgread, nzgblk;
	int blai, blaj, ibs, iblk, jblk;
	int *ig, *jg, *addrg;
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
	int *bsdia, *sp2blk, *bsblk, *addrscl, *addrsclinv;
	double *fmask, *dpiv, *diagg;
	double *lelem, *delem;
	double *gread, *aloc, *vloc, *eig, *work;

// Allocate and init array sp2blk

	sp2blk = new int [nsupc];
	if (!sp2blk) MemoryFail (funcname);

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			sp2blk[i] = iblk;
		};
	};

// Allocate and init working arrays including the scaling

// Integer arrays

	Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	Ich2AllocScale (_fout, 
					_nblks, sp2blk, _mtrafiles,
					gread, nzgread,
					fmask, dpiv, bsdia, diagg, 
					bsblk, _mtrlufiles, 
					addrscl, addrsclinv,
					lelem, delem, 
					aloc, vloc, eig, lwork, work, ops);

// Allocate factorization arrays

	Ich2AllocG (_param.memory, nzjgmx, ig, jg, addrg);

// Main cycle over the rows

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;
	nznew = 0;
	nzrtot = 0;
	nzgblk = 0;

	ig[0] = 0;

	for (irow=0;irow<nsupr;irow++) {

		iblk = sp2blk[irow];

		if (irow == _blks[iblk]) nzgblk = bsblk[iblk];

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			cout  << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			_fout << " Ich2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;

// Output hystogram of the scaling values

//			int itype = 2;
			int nloc = sprndr[irow];

//			SVHYSTDUM (&nloc, dpiv, &itype);
			csvhyst (cout, "Ich2Dia", nloc, dpiv);
			csvhyst (_fout,"Ich2Dia", nloc, dpiv);

		};

// Init current scaled row

		Ich2InitRow (irow, 
						sp2blk, _mtrafiles[iblk], _mtrlufiles,
						nlist1, lstloc, icycle, imask, 
						bsdia, diagg, fmask, addrscl, delem, aloc, vloc, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];
			ibs = addrg[j];

			blai = sprndr[irwprv+1]-sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = sprndr[icolmnl+1]-sprndr[icolmnl];

			jblk = sp2blk[irwprv];

//			for (kii=0;kii<blai*blaj;kii++) lelem[kii] = g[ibs+kii];

			FGet(_mtrlufiles[jblk],blai*blaj,lelem,ibs);

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			Ich2InitMask (elmisr, irow, irwprv, 
							iv, ig, jg, 
							nupd, lstupd, indlst, nlist1, lstloc, 
							icycle, imask, fmask);

// Perform updates

			Ich2UpdateRow (irow, 
							_mtrlufiles[jblk], gread,
							irwprv, nupd, lstupd, indlst, 
							lelem, fmask, addrg, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTransp (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		Ich2FiltrRow (irow, _param.theta, _param.tau2, 
						nlist1, lstloc, lstupd, 
						fmask, bsdia, diagg, iops);

// Compute and store the pivot value and scale current row

		Ich2PivScaleRow (irow, delem, nlist1, lstupd, 
							fmask, dpiv, 
							aloc, eig, lwork, work, iops);

// Count the local number of elements to be stored

		blai = sprndr[irow+1]-sprndr[irow];

		int nzloc=blai*blai;

		for (j=0;j<nlist1;j++) {
			jj = lstupd[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			nzloc += blai*blaj;
		};

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx) {

// Free the memory containing R part that is not in use any more

			Ich2FreeR (irow, ibegfr, nzjg, nzju, nzjrl, 
						iv, ig, jg, 
						nznew, addrg);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx) {
			cout  << " Insufficient workspace in Ich2, Irow = " << irow << endl;
			_fout << " Insufficient workspace in Ich2, Irow = " << irow << endl;
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

		Ich2StoreRow (irow,
						_mtrlufiles[iblk], aloc,
						_param.tau1, nzjg, nzju, nzjr, nzjrl, nzrtot,
						delem, nlist1, lstupd, fmask, 
						ig, jg, nznew, nzgblk, addrg);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = n;

	csvhyst (cout, "Ich2Dia", nloc, dpiv);
	csvhyst (_fout,"Ich2Dia", nloc, dpiv);

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

// Compute the size of the resulting matrix

	int nztot = 0;
	nzju = 0;

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		for (j=ig[i];j<ig[i+1];j++) {
			jj = jg[j];
			if (jj >= 0) {
				blaj = sprndr[jj+1]-sprndr[jj];
				nzju++;
				nztot += blai*blaj;
			};
		};
	};

// Perform the final filtering and scaling and compute matrix on return

	i = 0;

	CSMatrixRS temp (nsupr,nzju,i);

	Ich2StoreScaleU (_nblks, _blks, _mtrlufiles,
						sp2blk, bsdia, addrscl, addrsclinv, bsblk,
						ig, jg, addrg,
						temp, 
						aloc, vloc, delem, ops);

// Free factorization arrays 

	blai = sprndr[nsupr]-sprndr[nsupr-1];

	int nzmema = 2*nzatot-bsdia[nsupr-1]-blai*blai;

	delete [] addrscl;
	delete [] addrsclinv;
	delete [] bsdia;
	delete [] bsblk;
	delete [] sp2blk;
	delete [] lelem;
	delete [] delem;
	delete [] aloc;
	delete [] vloc;
	delete [] eig;
	delete [] work;

	delete [] ig;
	delete [] jg;
	delete [] addrg;
	delete [] gread;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Ich2 factorization statistics

	nz = nzatot;

	densu = 1.0e2 * (double) nztot / (double) nz;
	densr = 1.0e2 * (double) nzrtot / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	cout  << " Ich2 preconditioner generation statistics: " << endl;
	cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Ich2 preconditioner generation statistics: " << endl;
	_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

// Return the result

	return temp;

};
