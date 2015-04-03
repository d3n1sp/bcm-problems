//------------------------------------------------------------------------------------------------
// File: ilu2cs.cpp
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

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSMatrixCS: Allocate double complex data and init scaling for Ilu2
//========================================================================================
void CSMatrixCS::Ilu2AllocScale (ofstream &_fout, double _sclmin, // Allocate double complex data and init scaling for Ilu2
									dcmplx *&_fmaskl, dcmplx *&_fmasku, double *&_dpiv, 
									int *&_bsdia, dcmplx *&_diagg, 
									dcmplx *&_scll, dcmplx *&_sclu, dcmplx *&_sclinvl, dcmplx *&_sclinvu,
									dcmplx *&_lelem, dcmplx *&_uelem, dcmplx *&_deleml, dcmplx *&_delemu,
									dcmplx *&_aloc, dcmplx *&_uloc, dcmplx *&_vloc, 
									double *&_eig, int &_lwork, dcmplx *&_work, double *&_rwork, 
									double &_ops) const {

	const char *funcname = "Ilu2AllocScale";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone  (1.0e0,0.0e0);

	_fmaskl = new dcmplx [n*blamx];
	if (!_fmaskl) MemoryFail (funcname);
	_fmasku = new dcmplx [n*blamx];
	if (!_fmasku) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_bsdia = new int [nsupr];
	if (!_bsdia) MemoryFail (funcname);
	_diagg = new dcmplx [n*blamx];
	if (!_diagg) MemoryFail (funcname);
	_scll = new dcmplx [n*blamx];
	if (!_scll) MemoryFail (funcname);
	_sclu = new dcmplx [n*blamx];
	if (!_sclu) MemoryFail (funcname);
	_sclinvl = new dcmplx [n*blamx];
	if (!_sclinvl) MemoryFail (funcname);
	_sclinvu = new dcmplx [n*blamx];
	if (!_sclinvu) MemoryFail (funcname);
	_lelem = new dcmplx [blamx*blamx];
	if (!_lelem) MemoryFail (funcname);
	_uelem = new dcmplx [blamx*blamx];
	if (!_uelem) MemoryFail (funcname);
	_deleml = new dcmplx [blamx*blamx];
	if (!_deleml) MemoryFail (funcname);
	_delemu = new dcmplx [blamx*blamx];
	if (!_delemu) MemoryFail (funcname);

// Allocate working arrays

	int info=0;

	_lwork = 10 * blamx;

	_aloc = new dcmplx [blamx*blamx];
	if (!_aloc) MemoryFail (funcname);
	_uloc = new dcmplx [blamx*blamx];
	if (!_uloc) MemoryFail (funcname);
	_vloc = new dcmplx [blamx*blamx];
	if (!_vloc) MemoryFail (funcname);
	_eig = new double [blamx];
	if (!_eig) MemoryFail (funcname);
	_work = new dcmplx [_lwork];
	if (!_work) MemoryFail (funcname);
	_rwork = new double [_lwork];
	if (!_rwork) MemoryFail (funcname);

// Compute array bsdia

	int isup;
	int j, kii, kjj, kkk, kiib, ibssc;
	int blai, ibs, ibsv, ibsd, blai_2, kjjb;

	int nz=0;

	for (isup=0;isup<nsupr;isup++) {
		_bsdia[isup] = nz;
		blai = sprndr[isup+1]-sprndr[isup];
		nz += blai * blai;
	};

// Compute scaling and its inverse

	double *peig;

	dcmplx *paloc, *pa, *pscl;
	dcmplx *pfmask;

	dcmplx caux;
	double daux;

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

// Compute its singular value decomposition

		int blal = blai;

		if (blal != 1) {
			zgesvd_ ("A", "A", &blal, &blal,
						_aloc, &blal, _eig,
						_uloc, &blal, _vloc, &blal,
						_work, &_lwork, _rwork, &info);
			if (info != 0) throw " Error in the Lapack routine ZGESVD";
		} else {
			_eig[0] = abs(_aloc[0]);
			if (_eig[0] > 0.0e0) {
				_uloc[0] = _aloc[0] / _eig[0];
				_vloc[0] = cone;
			} else {
				_uloc[0] = cone;
				_vloc[0] = cone;
			};
		};

		_ops += 8 * blai_2*blai;

// Modify the pivots

		for (kii=0;kii<blai;kii++) {
			if (_eig[kii] < _sclmin) _eig[kii] = _sclmin;
		};

// Store scaling factors

		for (kii=0;kii<blai;kii++) {
			daux = _eig[kii];
			daux = sqrt(daux);
			_eig[kii] = daux;
			_dpiv[ibsv+kii] = _eig[kii];
		};

		for (kii=0;kii<blai;kii++) {
			peig = _eig;
			paloc = _uloc+kii;
			pscl = _sclinvl+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvl[ibsd+kjj*blai+kii] = _uloc[kjj*blai+kii] * _eig[kjj];
				*pscl = *paloc * *peig++;
				paloc += blai;
				pscl += blai;
			};
		};

		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			peig = _eig;
			paloc = _vloc+kiib;
			pscl = _sclinvu+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvu[ibsd+kjj*blai+kii] = _vloc[kii*blai+kjj] * _eig[kjj];
				*pscl = *paloc++ * *peig++;
				pscl += blai;
			};
			kiib += blai;
		};

		for (kii=0;kii<blai;kii++) {
			daux = _eig[kii];
			daux = 1.0e0 / daux;
			_eig[kii] = daux;
		};

		peig = _eig;
		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			paloc = _uloc+kiib;
			pscl = _scll+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_scll[ibsd+kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//				*pscl = conj(*paloc++) * *peig;
				dcmplx temp = conj(*paloc++);
				*pscl = temp * *peig;
				pscl += blai;
			};
			peig++;
			kiib += blai;
		};

		peig = _eig;
		for (kii=0;kii<blai;kii++) {
			paloc = _vloc+kii;
			pscl = _sclu+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclu[ibsd+kjj*blai+kii] = _vloc[kjj*blai+kii] * _eig[kii];
//				*pscl = conj(*paloc) * *peig;
				dcmplx temp = conj(*paloc);
				*pscl = temp * *peig;
				paloc += blai;
				pscl += blai;
			};
			peig++;
		};

		_ops += 4*blai_2 + 2*blai;

	};

// Reassign scaling factors if necessary

	if (false) {

		for (kii=0;kii<n*blamx;kii++) _scll[kii] = czero;
		for (kii=0;kii<n*blamx;kii++) _sclu[kii] = czero;
		for (kii=0;kii<n*blamx;kii++) _sclinvl[kii] = czero;
		for (kii=0;kii<n*blamx;kii++) _sclinvu[kii] = czero;

		for (isup=0;isup<nsupr;isup++) {
			blai = sprndr[isup+1]-sprndr[isup];
			ibsd = _bsdia[isup];
			for (kii=0;kii<blai;kii++) {
				_scll[ibsd+kii*blai+kii] = cone;
				_sclu[ibsd+kii*blai+kii] = cone;
				_sclinvl[ibsd+kii*blai+kii] = cone;
				_sclinvu[ibsd+kii*blai+kii] = cone;
			};
		};

	};

// Init working arrays

	pfmask = _fmaskl;
	for (kii=0;kii<n*blamx;kii++) *pfmask++ = czero;

	pfmask = _fmasku;
	for (kii=0;kii<n*blamx;kii++) *pfmask++ = czero;

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
				caux = czero;
				pa = a+ibs+kjjb;
				pscl = _scll+ibssc+kii;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//					caux += *pa++ * *pscl;
					caux.PlEq (*pa++ * *pscl);
					pscl += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*paloc = caux;
				paloc += blai;
				kjjb += blai;
			};
		};
		for (kii=0;kii<blai;kii++) {
			pa = _vloc+kii;
			for (kjj=0;kjj<blai;kjj++) {
				caux = czero;
				paloc = _aloc+kii;
				pscl = _sclu+ibssc+kjj;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _sclu[ibssc+kkk*blai+kjj];
//					caux += *paloc * *pscl;
					caux.PlEq (*paloc * *pscl);
					paloc += blai;
					pscl += blai;
				};
//				_vloc[kjj*blai+kii] = aux;
				*pa = caux;
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
// CSMatrixCS: Allocate G data for Ilu2
//========================================================================================
void CSMatrixCS::Ilu2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ilu2
							int *&_ig, int *&_jg, int *&_addrg, dcmplx *&_gl, dcmplx *&_gu) const {

	const char *funcname = "Ilu2AllocG";

	int nz = ia[nsupr];

	_nzjgmx = int (_memory * double(nz));
	_nzgmax = int (_memory * double(nza));

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_jg = new int [_nzjgmx];
	if (!_jg) MemoryFail (funcname);
	_addrg = new int [_nzjgmx];
	if (!_addrg) MemoryFail (funcname);
	_gl = new dcmplx [_nzgmax];
	if (!_gl) MemoryFail (funcname);
	_gu = new dcmplx [_nzgmax];
	if (!_gu) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixCS: Init current scaled row
//========================================================================================
void CSMatrixCS::Ilu2InitRow (const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau, // Init current scaled row
								int _irow, int &_nlist1, 
								int *_lstloc, int _icycle, int *_imask, 
								int *_bsdia, dcmplx *_diagg, dcmplx *_fmaskl, dcmplx *_fmasku, 
								dcmplx *_scll, dcmplx *_sclu,
								dcmplx *_aloc, dcmplx *_vloc,
								int &_iops) const {

//	const char *funcname = "Ilu2InitRow";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk, ibs, ibssc, jbssc;
	int kjjb, blai, blaj, ibsd, ibsv, jbsv, blaij, blai_2, blaj_2;
	dcmplx aux;

	dcmplx *paloc;
	dcmplx *pscloc;
	dcmplx *pfmask;
	dcmplx *pvloc;

	blai = _mtral.sprndr[_irow+1]-_mtral.sprndr[_irow];
	blai_2 = blai*blai;
	ibsd = _bsdia[_irow];
	ibsv = _mtral.sprndr[_irow];

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	pfmask = _fmaskl+ibsv*blai;
	pvloc = _diagg+ibsd;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	for (kii=0;kii<blai_2;kii++) *pfmask++ = *pvloc++;

	pfmask = _fmasku+ibsv*blai;
	pvloc = _diagg+ibsd;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	for (kii=0;kii<blai_2;kii++) *pfmask++ = *pvloc++;

	_nlist1++;

	for (int j=_mtral.ia[_irow]+1;j<_mtral.ia[_irow+1];j++) {
		int jj = _mtral.ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;

// Local scaling

		ibs = _mtral.bsa[j];
		ibssc = _bsdia[_irow];
		jbssc = _bsdia[jj];
		blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
		jbsv = _mtral.sprndr[jj];
		blaij = blai * blaj;
		blaj_2 = blaj*blaj;
		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				paloc = _mtrau.a+ibs+kjjb;
				pscloc = _scll+ibssc+kii;
				aux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//					aux += *paloc++ * *pscloc;
					aux.PlEq (*paloc++ * *pscloc);
					pscloc += blai;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				aux = czero;
				paloc = _aloc+kii;
				pscloc = _sclu+jbssc+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
//					aux += *paloc * *pscloc;
					aux.PlEq (*paloc * *pscloc);
					pscloc += blaj;
					paloc += blai;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		pfmask = _fmasku+jbsv*blai;
		pvloc = _vloc;
//		for (kii=0;kii<blai*blaj;kii++) _fmasku[jbsv*blai+kii] = _vloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmask++ = *pvloc++;

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				paloc = _mtral.a+ibs+kjjb;
				pscloc = _sclu+ibssc+kii;
				aux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _sclu[ibssc+kkk*blai+kii];
//					aux += *paloc++ * *pscloc;
					aux.PlEq (*paloc++ * *pscloc);
					pscloc += blai;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				aux = czero;
				paloc = _aloc+kii;
				pscloc = _scll+jbssc+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _scll[jbssc+kkk*blaj+kjj];
//					aux += *paloc * *pscloc;
					aux.PlEq (*paloc * *pscloc);
					pscloc += blaj;
					paloc += blai;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		pfmask = _fmaskl+jbsv*blai;
		pvloc = _vloc;
//		for (kii=0;kii<blai*blaj;kii++) _fmaskl[jbsv*blai+kii] = _vloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmask++ = *pvloc++;

		_nlist1++;

		_iops += 2*blai*(blaij+blaj_2);

	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Init mask and update arrays
//========================================================================================
void CSMatrixCS::Ilu2InitMask (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, dcmplx *_fmaskl, dcmplx *_fmasku) const {

//	const char *funcname = "Ilu2InitMask";

	dcmplx czero (0.0e0,0.0e0);

	int kii;
	int blai, blaj, jbsv, blaij;
	dcmplx *pfmask;

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
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blai*blaj;kii++) _fmaskl[jbsv*blai+kii] = 0.0e0;
			for (kii=0;kii<blaij;kii++) *pfmask++ = czero;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blai*blaj;kii++) _fmasku[jbsv*blai+kii] = 0.0e0;
			for (kii=0;kii<blaij;kii++) *pfmask++ = czero;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = k;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Update current row
//========================================================================================
void CSMatrixCS::Ich2UpdateRow (int _irow, int _irwprv, int _nupd, int *_lstupd, int *_indlst, // Update current row
							dcmplx *_lelem, dcmplx *_fmask, int *_addrg, dcmplx *_g, int &_iops) const {

//	const char *funcname = "Ich2UpdateRow";

	int kii, kjj, ibs, jbs;
	int kikb, kjkb;
	int blai, blaj, blak, jbsv, blaik;

	dcmplx *pfmloc;
	dcmplx *plloc;
	dcmplx *pgloc;

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
				dcmplx *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
//					*pfmloc -= *plloc * *pgloc;
					(*pfmloc).MnEq (*plloc * *pgloc);
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
// CSMatrixCS: Perform filtering of the current row
//========================================================================================
void CSMatrixCS::Ilu2FiltrRow (int _irow, int _fcttyp, double _theta, double _tau2, // Perform filtering of the current row
								int *_iaini, int *_jaini, int &_icycleini, int *_imaskini,
								int &_nlist1, int *_lstloc, int *_lstupd,
								dcmplx *_fmaskl, dcmplx* _fmasku, int *_bsdia, dcmplx *_diagg, int &_iops) const {

//	const char *funcname = "Ilu2FiltrRow";

	int nlist2 = 0;
	int kii, kjj, blai, blaj, ibsv, jbsv, jbsd, blaij;
	dcmplx *pfmaskl, *pfmasku, *pfmaskuu, *pdiag;

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsv = sprndr[_irow];

	bool large_elem;

	if (_fcttyp == 0) {
		int j, jj;
		_icycleini++;
		for (j=_iaini[_irow];j<_iaini[_irow+1];j++) {
			jj = _jaini[j];
			_imaskini[jj] = _icycleini;
		};
	};

	for (int i=1;i<_nlist1;i++) {

		int j = _lstloc[i];

		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		jbsd = _bsdia[j];
		blaij = blai*blaj;

		double aux = 0.0e0;

		pfmaskl = _fmaskl+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmaskl[jbsv*blai+kii];
			double auxloc = abs(*pfmaskl++);
			if (auxloc > aux) aux = auxloc;
		};
		pfmasku = _fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = abs(*pfmasku++);
			if (auxloc > aux) aux = auxloc;
		};

		large_elem = false;
		if (_fcttyp == 0) {
			if (_imaskini[j] == _icycleini) large_elem = true;
		} else {
			if (aux >= _tau2) {
				large_elem = true;
			};
		};

		if (large_elem) {
			_lstupd[nlist2] = j;
			nlist2++;
		} else {
			pfmaskuu = _fmasku+ibsv*blai;
			for (kii=0;kii<blai;kii++) {
				pfmaskl = _fmaskl+jbsv*blai+kii;
				pfmasku = _fmasku+jbsv*blai+kii;
				pdiag = _diagg+jbsd;
				for (kjj=0;kjj<blaj;kjj++) {
//					double auxloc = _fmaskl[jbsv*blai+kjj*blai+kii];
					double auxloc = abs(*pfmaskl);
//					double auxloc1 = _fmasku[jbsv*blai+kjj*blai+kii];
					double auxloc1 = abs(*pfmasku);
					if (auxloc1 > auxloc) auxloc = auxloc1;
					auxloc *= _theta;
//					_fmasku[ibsv*blai+kii*blai+kii] += auxloc;
					*pfmaskuu += auxloc;
//					_diagg[jbsd+kjj*blaj+kjj] += auxloc;
					*pdiag += auxloc;
					if (kjj<(blaj-1)) {
						pfmaskl += blai;
						pfmasku += blai;
						pdiag += blaj+1;
					};
				};
				if (kii<(blai-1)) {
					pfmaskuu += blai+1;
				};
			};
		};

	};

	_nlist1 = nlist2;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Compute diagonal pivot and scale current row
//========================================================================================
void CSMatrixCS::Ilu2PivScaleRow (double _pivmin, // Compute diagonal pivot and scale current row
									int _irow, dcmplx *_deleml, dcmplx *_delemu, 
									int _nlist1, int *_lstupd, dcmplx *_fmaskl, dcmplx *_fmasku, double *_dpiv, 
									dcmplx *_aloc, double *_eig, dcmplx *_uloc, dcmplx *_vloc, 
									int _lwork, dcmplx *_work, double *_rwork, int &_iops) const {

//	const char *funcname = "Ilu2PivScaleRow";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone  (1.0e0,0.0e0);

	int kii, kjj, kkk, ibs;
	int blai, blaj, ibsv, jbsv;
	dcmplx caux;

// Copy diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsv = sprndr[_irow];
	ibs = ibsv * blai;

	for (kii=0;kii<blai*blai;kii++) _aloc[kii] = _fmasku[ibs+kii];

// Compute singular value decomposition of the diagonal supernode

	int blal = blai;
	int info = 0;

	if (blal != 1) {
		zgesvd_ ("A", "A", &blal, &blal,
					_aloc, &blal, _eig, 
					_uloc, &blal, _vloc, &blal, 
					_work, &_lwork, _rwork, &info);
		if (info != 0) throw " Error in the Lapack routine ZGESVD";
	} else {
		_eig[0] = abs(_aloc[0]);
		if (_eig[0] > 0.0e0) {
			_uloc[0] = _aloc[0] / _eig[0];
			_vloc[0] = cone;
		} else {
			_uloc[0] = cone;
			_vloc[0] = cone;
		};
	};

	_iops += 8 * blai * blai * blai;

// Modify singular values if necessary

	for (kii=0;kii<blai;kii++) {
		if (_eig[kii] < _pivmin) _eig[kii] = _pivmin;
		_eig[kii] = sqrt(_eig[kii]);
		_dpiv[ibsv+kii] = _eig[kii];
		_eig[kii] = 1.0e0/_eig[kii];
	};

// Compute diagonal supernodes

	int kiib;
	dcmplx *pdeloc;
	dcmplx *pfmloc;
	dcmplx *palloc;

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = _deleml+kii;
		palloc = _uloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_deleml[kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//			*pdeloc = conj(*palloc++) * _eig[kii];
			dcmplx temp = conj(*palloc++);
			*pdeloc = temp * _eig[kii];
			pdeloc += blai;
		};
		kiib += blai;
	};

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = _delemu+kii;
		palloc = _vloc+kii;
//		palloc = _vloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_delemu[kjj*blai+kii] = _vloc[kjj*blai+kii] * _eig[kii];
//			*pdeloc = conj(*palloc) * _eig[kii];
			dcmplx temp = conj(*palloc);
			*pdeloc = temp * _eig[kii];
			pdeloc += blai;
			palloc += blai;
//			palloc++;
		};
		kiib += blai;
	};

	_iops += 2 * blai * blai;

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
				pfmloc = _fmaskl+ibs+kjjb;
				pdeloc = _delemu+kii;
				caux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _delemu[kkk*blai+kii] * _fmaskl[ibs+kjj*blai+kkk];
//					caux += *pdeloc * *pfmloc++;
					caux.PlEq (*pdeloc * *pfmloc++);
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = caux;
				palloc += blai;
			};
		};

		pfmloc = _fmaskl+ibs;
		palloc = _aloc;
//		for (kii=0;kii<blaij;kii++) _fmaskl[ibs+kii] = _aloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmloc++ = *palloc++;

		for (kii=0;kii<blai;kii++) {
			palloc = _aloc+kii;
			for (kjj=0;kjj<blaj;kjj++) {
				kjjb = kjj*blai;
				pfmloc = _fmasku+ibs+kjjb;
				pdeloc = _deleml+kii;
				caux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _deleml[kkk*blai+kii] * _fmasku[ibs+kjj*blai+kkk];
//					caux += *pdeloc * *pfmloc++;
					caux.PlEq (*pdeloc * *pfmloc++);
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = caux;
				palloc += blai;
			};
		};

		pfmloc = _fmasku+ibs;
		palloc = _aloc;
//		for (kii=0;kii<blaij;kii++) _fmasku[ibs+kii] = _aloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmloc++ = *palloc++;

		_iops += 2 * blai * blaij;

	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Free R part of the data
//========================================================================================
void CSMatrixCS::Ilu2FreeR (ofstream &_fout, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, int &_nznew, 
							int *_addrg, dcmplx *_gl, dcmplx *_gu) const {

//	const char *funcname = "Ilu2FreeR";

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
	dcmplx *pgold, *pgnew;

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
				pgold = _gl+ibs;
				pgnew = _gl+_nznew;
//				for (int kii=0;kii<blaij;kii++) _gl[_nznew+kii] = _gl[ibs+kii];
				int kii;
				for (kii=0;kii<blaij;kii++) *pgnew++ = *pgold++;
				pgold = _gu+ibs;
				pgnew = _gu+_nznew;
//				for (int kii=0;kii<blaij;kii++) _gu[_nznew+kii] = _gu[ibs+kii];
				for (kii=0;kii<blaij;kii++) *pgnew++ = *pgold++;
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
// CSMatrixCS: Store current row
//========================================================================================
void CSMatrixCS::Ilu2StoreRow (int _irow, int _fcttyp, double _tau1, // Store current row
											int *_iaini, int *_jaini, int &_icycleini, int *_imaskini,
											int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
											dcmplx *_deleml, dcmplx *_delemu, int _nlist1, int *_lstupd, dcmplx *_fmaskl, dcmplx *_fmasku, 
											int *_ig, int *_jg, int &_nznew, int *_addrg, dcmplx *_gl, dcmplx *_gu) const {

//	const char *funcname = "Ilu2StoreRow";

//	int isupchk=11;

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	dcmplx *pg, *pfmask, *pd;

// Store current diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	_jg[_nzjg] = _irow;
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = _deleml+kiib;
		pg = _gl+_nznew+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gl [_nznew+kjj*blai+kii] = _deleml[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = _delemu+kiib;
		pg = _gu+_nznew+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gu [_nznew+kjj*blai+kii] = _delemu[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};
	_addrg[_nzjg] = _nznew;
	_nzjg++;
	_nzju++;
	_nznew += blai * blai;

	int i;

	if (_fcttyp == 0) {
		_icycleini++;
		int j, jj;
		for (j=_iaini[_irow];j<_iaini[_irow+1];j++) {
			jj = _jaini[j];
			_imaskini[jj] = _icycleini;
		};
	};

	bool firstorder;

	for (i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		blaij = blai*blaj;
		double aux = 0.0e0;
		pfmask = _fmaskl+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmaskl[jbsv*blai+kii];
			double auxloc = abs(*pfmask++);
			if (auxloc > aux) aux = auxloc;
		};
		pfmask = _fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = abs(*pfmask++);
			if (auxloc > aux) aux = auxloc;
		};
		_addrg[_nzjg] = _nznew;
		firstorder = false;
		if (_fcttyp == 0) {
			if (_imaskini[j] == _icycleini) firstorder = true;
		} else {
			if (aux >= _tau1) firstorder = true;
		};
		if (firstorder) {
			_jg[_nzjg] = j;
			pg = _gl+_nznew;
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			pg = _gu+_nznew;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
			pg = _gl+_nznew;
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			pg = _gu+_nznew;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
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
// CSMatrixCS: Store the resulting U and scale it
//========================================================================================
void CSMatrixCS::Ich2StoreScaleU (int *_bsdia, dcmplx *_scl, dcmplx *_sclinv, // Store the resulting U and scale it
									int *_ig, int *_jg, int *_addrg, dcmplx *_g,
									CSMatrixCS &_temp, 
									dcmplx *_aloc,
									double &_ops) const {

//	const char *funcname = "Ich2StoreScaleU";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj, blaij;
	dcmplx caux;

	dcmplx *pscloc;
	dcmplx *pgloc;
	dcmplx *paloc;

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
							caux = czero;
							for (kkk=0;kkk<blai;kkk++) {
//								aux += _scl[ibssc+kii*blai+kkk] * _g[jbs+kjj*blai+kkk];
//								caux += *pscloc * *pgloc;
								caux.PlEq (*pscloc * *pgloc);
								pgloc++;
								pscloc++;
							};
//							_aloc[kjjb+kii] = caux;
							*paloc = caux;
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
							caux = czero;
							pgloc = _g+jbs+kii;
							pscloc = _sclinv+jjbssc+kjj;
							for (kkk=0;kkk<blaj;kkk++) {
//								aux += _g[jbs+kkk*blai+kii] * _sclinv[jjbs+kkk*blaj+kjj];
//								caux += *pgloc * *pscloc;
								caux.PlEq (*pgloc * *pscloc);
								pgloc += blai;
								pscloc += blaj;
							};
//							_aloc[kjjb+kii] = aux;
							*paloc = caux;
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
// CSMatrixCS: Second order Incomplete LU decomposition
//========================================================================================
void CSMatrixCS::Ilu2 (ofstream &_fout, const CSlvParam _param, // Second order Incomplete LU decomposition
								const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau,
								CSMatrixCS &_mtrl, CSMatrixCS &_mtru) { 

	const char *funcname = "Ilu2";

// Local variables

	int irow, i, j, jj, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzgmax, nzjgmx;
	int nz, nznew, nzju, nzjr, nzjrl, nzrtot, icheck, ibegfr;
	int kii, blai, blaj, ibs;
	int *ig, *jg, *addrg;
	dcmplx *gl, *gu;
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
	int icycleini, *imaskini;
	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	int *bsdia;
	dcmplx *fmaskl, *fmasku, *diagg, *scll, *sclu, *sclinvl, *sclinvu;
	dcmplx *lelem, *uelem, *deleml, *delemu;
	dcmplx *aloc, *uloc, *vloc, *work;
	double *dpiv, *eig, *rwork;

// Create dummy matrix

	CSMatrixCS mtrdummy;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

	imaskini = new int [_mtral.nsupr];
	if (!imaskini) MemoryFail (funcname);

	icycleini = -1;
	for (irow=0;irow<_mtral.nsupr;irow++) imaskini[irow] = -1;

// Double arrays and compute scaling

	_mtrau.Ilu2AllocScale (_fout, _param.sclmin, 
							fmaskl, fmasku, dpiv, bsdia, diagg, 
							scll, sclu, sclinvl, sclinvu, 
							lelem, uelem, deleml, delemu, 
							aloc, uloc, vloc, eig, 
							lwork, work, rwork, ops);

// Allocate factorization arrays

	_mtral.Ilu2AllocG (_param.memory, nzjgmx, nzgmax, ig, jg, addrg, gl, gu);

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

	int *pia = _mtral.GetIa ();
	int *pja = _mtral.GetJa ();

	for (irow=0;irow<_mtral.nsupr;irow++) {

		blai = _mtral.sprndr[irow+1]-_mtral.sprndr[irow];

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			cout  << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			_fout << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;

// Output hystogram of the scaling values

//			int itype = 2;
			int nloc = _mtral.sprndr[irow];

			csvhyst (cout, "Ilu2Piv", nloc, dpiv);
			csvhyst (_fout,"Ilu2Piv", nloc, dpiv);

		};

// Init current scaled row

		mtrdummy.Ilu2InitRow (_mtral, _mtrau, 
								irow, nlist1, lstloc, icycle, imask, 
								bsdia, diagg, fmaskl, fmasku, 
								scll, sclu, 
								aloc, vloc, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];
			ibs = addrg[j];

			int blailoc = _mtral.sprndr[irwprv+1]-_mtral.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _mtral.sprndr[icolmnl+1]-_mtral.sprndr[icolmnl];

			for (kii=0;kii<blailoc*blaj;kii++) lelem[kii] = gl[ibs+kii];
			for (kii=0;kii<blailoc*blaj;kii++) uelem[kii] = gu[ibs+kii];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			_mtral.Ilu2InitMask (elmisr, irow, irwprv, 
									iv, ig, jg, 
									nupd, lstupd, indlst, nlist1, lstloc, 
									icycle, imask, fmaskl, fmasku);

// Perform updates

			_mtral.Ich2UpdateRow (irow, irwprv, nupd, lstupd, indlst, 
									lelem, fmasku, addrg, gu, iops);
			_mtral.Ich2UpdateRow (irow, irwprv, nupd, lstupd, indlst, 
									uelem, fmaskl, addrg, gl, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTransp (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		_mtral.Ilu2FiltrRow (irow, _param.fcttyp, _param.theta, _param.tau2, 
									pia, pja, icycleini, imaskini, 
									nlist1, lstloc, lstupd, 
									fmaskl, fmasku, bsdia, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin,
								irow, deleml, delemu, nlist1, lstupd, 
								fmaskl, fmasku, dpiv, 
								aloc, eig, uloc, vloc, lwork, work, rwork, iops);

// Count the local number of elements to be stored

		blai = _mtral.sprndr[irow+1]-_mtral.sprndr[irow];

		int nzloc=blai*blai;

		for (j=0;j<nlist1;j++) {
			jj = lstupd[j];
			blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
			nzloc += blai*blaj;
		};

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx || nznew+nzloc > nzgmax) {

// Free the memory containing R part that is not in use any more

			_mtral.Ilu2FreeR (_fout, irow, ibegfr, nzjg, nzju, nzjrl, 
								iv, ig, jg, 
								nznew, addrg, gl, gu);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx || nznew+nzloc > nzgmax) {
			cout  << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
			_fout << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
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

		_mtral.Ilu2StoreRow (irow, _param.fcttyp, _param.tau1, 
								pia, pja, icycleini, imaskini, 
								nzjg, nzju, nzjr, nzjrl, nzrtot,
								deleml, delemu, nlist1, lstupd, fmaskl, fmasku, 
								ig, jg, nznew, addrg, gl, gu);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = _mtral.n;

	csvhyst (cout, "Ilu2Piv", nloc, dpiv);
	csvhyst (_fout,"Ilu2Piv", nloc, dpiv);

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] imaskini;

	delete [] fmaskl;
	delete [] fmasku;
	delete [] dpiv;
	delete [] diagg;

// Compute the size of the resulting matrix

	int nztot = 0;
	nzju = 0;

	for (i=0;i<_mtral.nsupr;i++) {
		blai = _mtral.sprndr[i+1]-_mtral.sprndr[i];
		for (j=ig[i];j<ig[i+1];j++) {
			jj = jg[j];
			if (jj >= 0) {
				blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
				nzju++;
				nztot += blai*blaj;
			};
		};
	};

// Perform the final filtering and scaling and compute matrix on return

	CSMatrixCS temp (_mtral.nsupr,nzju,nztot);

	_mtrl = temp;

	temp = mtrdummy;

	_mtru = _mtrl;

	_mtral.Ich2StoreScaleU (bsdia, scll, sclinvl, 
							ig, jg, addrg, gl, 
							_mtrl, 
							aloc, ops);
	_mtral.Ich2StoreScaleU (bsdia, sclu, sclinvu, 
							ig, jg, addrg, gu, 
							_mtru, 
							aloc, ops);

// Free factorization arrays 

	blai = _mtral.sprndr[_mtral.nsupr]-_mtral.sprndr[_mtral.nsupr-1];

	int nzmema = 2*_mtral.nza-bsdia[_mtral.nsupr-1]-blai*blai;

	delete [] bsdia;
	delete [] scll;
	delete [] sclu;
	delete [] sclinvl;
	delete [] sclinvu;
	delete [] lelem;
	delete [] uelem;
	delete [] deleml;
	delete [] delemu;
	delete [] aloc;
	delete [] vloc;
	delete [] uloc;
	delete [] eig;
	delete [] work;
	delete [] rwork;

	delete [] ig;
	delete [] jg;
	delete [] addrg;
	delete [] gl;
	delete [] gu;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Ich2 factorization statistics

	nz = _mtral.nza;

	densu = 1.0e2 * (double) nztot / (double) nz;
	densr = 1.0e2 * (double) nzrtot / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	cout  << " Ilu2 preconditioner generation statistics: " << endl;
	cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Ilu2 preconditioner generation statistics: " << endl;
	_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Allocate double complex data and init scaling for Ilu2
//========================================================================================
void CSMatrixCS::Ilu2AllocScale (ofstream &_fout, double _sclmin, // Allocate double complex data and init scaling for Ilu2
									int _nblks, int *_sp2blk, FILE **_mtraufiles,
									dcmplx *&_fmaskl, dcmplx *&_fmasku, double *&_dpiv, 
									int *&_bsdia, dcmplx *&_diagg, 
									int *&_bsblk, FILE **_mtrlfiles, FILE **_mtrufiles,
									int *&_addrscll, int *&_addrsclu, int *&_addrsclinvl, int *&_addrsclinvu,
									dcmplx *&_gread, dcmplx *&_lelem, dcmplx *&_uelem, dcmplx *&_deleml, dcmplx *&_delemu,
									dcmplx *&_aloc, dcmplx *&_uloc, dcmplx *&_vloc, 
									double *&_eig, int &_lwork, dcmplx *&_work, double *&_rwork, 
									double &_ops) const {

	const char *funcname = "Ilu2AllocScale";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone  (1.0e0,0.0e0);

	_fmaskl = new dcmplx [n*blamx];
	if (!_fmaskl) MemoryFail (funcname);
	_fmasku = new dcmplx [n*blamx];
	if (!_fmasku) MemoryFail (funcname);
	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_bsdia = new int [nsupr];
	if (!_bsdia) MemoryFail (funcname);
	_diagg = new dcmplx [n*blamx];
	if (!_diagg) MemoryFail (funcname);

	_bsblk = new int [_nblks];
	if (!_bsblk) MemoryFail (funcname);

	_addrscll = new int [nsupr];
	if (!_addrscll) MemoryFail (funcname);
	_addrsclu = new int [nsupr];
	if (!_addrsclu) MemoryFail (funcname);
	_addrsclinvl = new int [nsupr];
	if (!_addrsclinvl) MemoryFail (funcname);
	_addrsclinvu = new int [nsupr];
	if (!_addrsclinvu) MemoryFail (funcname);

	_gread = new dcmplx [blamx*blamx];
	if (!_gread) MemoryFail (funcname);
	_lelem = new dcmplx [blamx*blamx];
	if (!_lelem) MemoryFail (funcname);
	_uelem = new dcmplx [blamx*blamx];
	if (!_uelem) MemoryFail (funcname);
	_deleml = new dcmplx [blamx*blamx];
	if (!_deleml) MemoryFail (funcname);
	_delemu = new dcmplx [blamx*blamx];
	if (!_delemu) MemoryFail (funcname);

// Allocate working arrays

	int info=0;

	_lwork = 10 * blamx;

	_aloc = new dcmplx [blamx*blamx];
	if (!_aloc) MemoryFail (funcname);
	_uloc = new dcmplx [blamx*blamx];
	if (!_uloc) MemoryFail (funcname);
	_vloc = new dcmplx [blamx*blamx];
	if (!_vloc) MemoryFail (funcname);
	_eig = new double [blamx];
	if (!_eig) MemoryFail (funcname);
	_work = new dcmplx [_lwork];
	if (!_work) MemoryFail (funcname);
	_rwork = new double [_lwork];
	if (!_rwork) MemoryFail (funcname);

// Compute array bsdia

	int isup, iblk;
	int j, kii, kjj, kkk, kiib, ibssc;
	int blai, ibs, ibsv, ibsd, blai_2, kjjb;

	int nz=0;

	for (isup=0;isup<nsupr;isup++) {
		_bsdia[isup] = nz;
		blai = sprndr[isup+1]-sprndr[isup];
		nz += blai * blai;
	};

	for (iblk=0;iblk<_nblks;iblk++) _bsblk[iblk] = 0;

// Compute scaling and its inverse

	double *peig;

	dcmplx *paloc, *pa, *pscl;
	dcmplx *pfmask;

	dcmplx caux;
	double daux;

	for (isup=0;isup<nsupr;isup++) {

		blai = sprndr[isup+1]-sprndr[isup];
		blai_2 = blai*blai;
		ibsv = sprndr[isup];
		ibsd = _bsdia[isup];

// Assign diagonal supernode

		iblk = _sp2blk[isup];
		j = ia[isup];
		ibs = bsa[j];

		paloc = _aloc;
//		pa = a+ibs;
//		for (kii=0;kii<blai_2;kii++) _aloc[kii] = a[ibs+kii];
//		for (kii=0;kii<blai_2;kii++) *paloc++ = *pa++;

		FGet (_mtraufiles[iblk],blai*blai,_aloc,ibs);

// Compute its singular value decomposition

		int blal = blai;

		zgesvd_ ("A", "A", &blal, &blal,
					_aloc, &blal, _eig, 
					_uloc, &blal, _vloc, &blal, 
					_work, &_lwork, _rwork, &info);

		if (info != 0) throw " Error in the Lapack routine ZGESVD";

		_ops += 8 * blai_2*blai;

// Modify the pivots

		for (kii=0;kii<blai;kii++) {
			if (_eig[kii] < _sclmin) _eig[kii] = _sclmin;
		};

// Store scaling factors

		for (kii=0;kii<blai;kii++) {
			daux = _eig[kii];
			daux = sqrt(daux);
			_eig[kii] = daux;
			_dpiv[ibsv+kii] = _eig[kii];
		};

		for (kii=0;kii<blai;kii++) {
			peig = _eig;
			paloc = _uloc+kii;
			pscl = _deleml+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvl[ibsd+kjj*blai+kii] = _uloc[kjj*blai+kii] * _eig[kjj];
				*pscl = *paloc * *peig++;
				paloc += blai;
				pscl += blai;
			};
		};

		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			peig = _eig;
			paloc = _vloc+kiib;
			pscl = _delemu+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclinvu[ibsd+kjj*blai+kii] = _vloc[kii*blai+kjj] * _eig[kjj];
				*pscl = *paloc++ * *peig++;
				pscl += blai;
			};
			kiib += blai;
		};

		ibs = _bsblk[iblk];

		FPut (_mtrlfiles[iblk],blai*blai,_deleml,ibs);
		FPut (_mtrufiles[iblk],blai*blai,_delemu,ibs);

		_addrsclinvl[isup] = ibs;
		_addrsclinvu[isup] = ibs;

		_bsblk[iblk] += blai*blai;

		for (kii=0;kii<blai;kii++) {
			daux = _eig[kii];
			daux = 1.0e0 / daux;
			_eig[kii] = daux;
		};

		peig = _eig;
		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			paloc = _uloc+kiib;
			pscl = _deleml+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_scll[ibsd+kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
//				*pscl = conj(*paloc++) * *peig;
				dcmplx temp = conj (*paloc++);
				*pscl = temp * *peig;
				pscl += blai;
			};
			peig++;
			kiib += blai;
		};

		peig = _eig;
		for (kii=0;kii<blai;kii++) {
			paloc = _vloc+kii;
			pscl = _delemu+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_sclu[ibsd+kjj*blai+kii] = _vloc[kjj*blai+kii] * _eig[kii];
//				*pscl = conj(*paloc) * *peig;
				dcmplx temp = conj(*paloc);
				*pscl = temp * *peig;
				paloc += blai;
				pscl += blai;
			};
			peig++;
		};

		ibs = _bsblk[iblk];

		FPut (_mtrlfiles[iblk],blai*blai,_deleml,ibs);
		FPut (_mtrufiles[iblk],blai*blai,_delemu,ibs);

		_addrscll[isup] = ibs;
		_addrsclu[isup] = ibs;

		_bsblk[iblk] += blai*blai;

		_ops += 4*blai_2 + 2*blai;

	};

// Init working arrays

	pfmask = _fmaskl;
	for (kii=0;kii<n*blamx;kii++) *pfmask++ = czero;

	pfmask = _fmasku;
	for (kii=0;kii<n*blamx;kii++) *pfmask++ = czero;

	for (isup=0;isup<nsupr;isup++) {
		blai = sprndr[isup+1]-sprndr[isup];
		blai_2 = blai*blai;
		j = ia[isup];
		ibs = bsa[j];
		ibssc = _bsdia[isup];

		iblk = _sp2blk[isup];

		ibs = bsa[j];

		FGet (_mtraufiles[iblk],blai*blai,_gread,ibs);

		ibs = _addrscll[isup];

		FGet (_mtrlfiles[iblk],blai*blai,_deleml,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			paloc = _aloc+kii;
			for (kjj=0;kjj<blai;kjj++) {
				caux = czero;
				pa = _gread+kjjb;
				pscl = _deleml+kii;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//					caux += *pa++ * *pscl;
					caux.PlEq (*pa++ * *pscl);
					pscl += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*paloc = caux;
				paloc += blai;
				kjjb += blai;
			};
		};

		ibs = _addrsclu[isup];

		FGet (_mtrufiles[iblk],blai*blai,_delemu,ibs);

		for (kii=0;kii<blai;kii++) {
			pa = _vloc+kii;
			for (kjj=0;kjj<blai;kjj++) {
				caux = czero;
				paloc = _aloc+kii;
				pscl = _delemu+kjj;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _sclu[ibssc+kkk*blai+kjj];
//					caux += *paloc * *pscl;
					caux.PlEq (*paloc * *pscl);
					paloc += blai;
					pscl += blai;
				};
//				_vloc[kjj*blai+kii] = aux;
				*pa = caux;
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
// CSMatrixCS: Allocate G data for Ilu2
//========================================================================================
void CSMatrixCS::Ilu2AllocG (double _memory, int &_nzjgmx, // Allocate G data for Ilu2
							int *&_ig, int *&_jg, int *&_addrg) const {

	const char *funcname = "Ilu2AllocG";

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
// CSMatrixCS: Init current scaled row
//========================================================================================
void CSMatrixCS::Ilu2InitRow (int *_sp2blk, FILE *_mtralfile, FILE *_mtraufile, // Init current scaled row
								const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau, 
								FILE **_mtrlfiles, FILE **_mtrufiles,
								int _irow, int &_nlist1, 
								int *_lstloc, int _icycle, int *_imask, 
								int *_bsdia, dcmplx *_diagg, dcmplx *_fmaskl, dcmplx *_fmasku, dcmplx *_gread,
								int *_addrscll, int *_addrsclu,
								dcmplx *_aloc, dcmplx *_vloc,
								int &_iops) const {

//	const char *funcname = "Ilu2InitRow";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk, ibs, ibssc, jbssc, iblk, jblk;
	int kjjb, blai, blaj, ibsd, ibsv, jbsv, blaij, blai_2, blaj_2;
	dcmplx aux;

	dcmplx *paloc;
	dcmplx *pscloc;
	dcmplx *pfmask;
	dcmplx *pvloc;

	iblk = _sp2blk[_irow];

	blai = _mtral.sprndr[_irow+1]-_mtral.sprndr[_irow];
	blai_2 = blai*blai;
	ibsd = _bsdia[_irow];
	ibsv = _mtral.sprndr[_irow];

	_nlist1 = 0;

	_lstloc[_nlist1] = _irow;
	_imask[_irow] = _icycle;
	pfmask = _fmaskl+ibsv*blai;
	pvloc = _diagg+ibsd;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	for (kii=0;kii<blai_2;kii++) *pfmask++ = *pvloc++;

	pfmask = _fmasku+ibsv*blai;
	pvloc = _diagg+ibsd;
//	for (kii=0;kii<blai_2;kii++) _fmask[ibsv*blai+kii] = _diagg[ibsd+kii];
	for (kii=0;kii<blai_2;kii++) *pfmask++ = *pvloc++;

	_nlist1++;

	for (int j=_mtral.ia[_irow]+1;j<_mtral.ia[_irow+1];j++) {
		int jj = _mtral.ja[j];
		_lstloc[_nlist1] = jj;
		_imask[jj] = _icycle;

// Local scaling

		jblk = _sp2blk[jj];

		ibssc = _bsdia[_irow];
		jbssc = _bsdia[jj];
		blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
		jbsv = _mtral.sprndr[jj];
		blaij = blai * blaj;
		blaj_2 = blaj*blaj;

		ibs = _mtrau.bsa[j];

		FGet (_mtraufile,blaij,_vloc,ibs);

		ibs = _addrscll[_irow];

		FGet (_mtrlfiles[iblk],blai*blai,_gread,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				paloc = _vloc+kjjb;
				pscloc = _gread+kii;
				aux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
//					aux += *paloc++ * *pscloc;
					aux.PlEq (*paloc++ * *pscloc);
					pscloc += blai;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};

		ibs = _addrsclu[jj];

		FGet (_mtrufiles[jblk],blaj*blaj,_gread,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				aux = czero;
				paloc = _aloc+kii;
				pscloc = _gread+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
//					aux += *paloc * *pscloc;
					aux.PlEq (*paloc * *pscloc);
					pscloc += blaj;
					paloc += blai;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		pfmask = _fmasku+jbsv*blai;
		pvloc = _vloc;
//		for (kii=0;kii<blai*blaj;kii++) _fmasku[jbsv*blai+kii] = _vloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmask++ = *pvloc++;

		ibs = _mtral.bsa[j];

		FGet (_mtralfile,blaij,_vloc,ibs);

		ibs = _addrsclu[_irow];

		FGet (_mtrufiles[iblk],blai*blai,_gread,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				paloc = _vloc+kjjb;
				pscloc = _gread+kii;
				aux = czero;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _sclu[ibssc+kkk*blai+kii];
//					aux += *paloc++ * *pscloc;
					aux.PlEq (*paloc++ * *pscloc);
					pscloc += blai;
				};
				_aloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};

		ibs = _addrscll[jj];

		FGet (_mtrlfiles[jblk],blaj*blaj,_gread,ibs);

		for (kii=0;kii<blai;kii++) {
			kjjb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				aux = czero;
				paloc = _aloc+kii;
				pscloc = _gread+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _scll[jbssc+kkk*blaj+kjj];
//					aux += *paloc * *pscloc;
					aux.PlEq (*paloc * *pscloc);
					pscloc += blaj;
					paloc += blai;
				};
				_vloc[kjjb+kii] = aux;
				kjjb += blai;
			};
		};
		pfmask = _fmaskl+jbsv*blai;
		pvloc = _vloc;
//		for (kii=0;kii<blai*blaj;kii++) _fmaskl[jbsv*blai+kii] = _vloc[kii];
		for (kii=0;kii<blaij;kii++) *pfmask++ = *pvloc++;

		_nlist1++;

		_iops += 2*blai*(blaij+blaj_2);

	};

};

// Author: Kharchenko S.A.
// CSMatrixCS: Update current row
//========================================================================================
void CSMatrixCS::Ich2UpdateRow (int _irow, // Update current row
								FILE *_mtrlufile, dcmplx *_gread,
								int _irwprv, int _nupd, int *_lstupd, int *_indlst,
								dcmplx *_lelem, dcmplx *_fmask, int *_addrg, int &_iops) const {

//	const char *funcname = "Ich2UpdateRow";

	int kii, kjj, ibs, jbs;
	int kikb, kjkb;
	int blai, blaj, blak, jbsv, blaik;

	dcmplx *pfmloc;
	dcmplx *plloc;
	dcmplx *pgloc;

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

		FGet (_mtrlufile,blaj*blak,_gread,jbs);

		kikb = 0;
		for (kii=0;kii<blai;kii++) {
			pfmloc = _fmask+ibs+kii;
			kjkb = 0;
			for (kjj=0;kjj<blaj;kjj++) {
				plloc = _lelem+kikb;
				pgloc = _gread+kjkb;
				dcmplx *end = plloc+blak;
				for (;plloc<end;plloc++,pgloc++) {
//					_fmask[ibs+kjj*blai+kii] -= _lelem[kii*blak+kkk] * _g[jbs+kjj*blak+kkk];
//					*pfmloc -= *plloc * *pgloc;
					(*pfmloc).MnEq (*plloc * *pgloc);
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
// CSMatrixCS: Free R part of the data
//========================================================================================
void CSMatrixCS::Ilu2FreeR (ofstream &_fout, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, int &_nznew, 
							int *_addrg) const {

//	const char *funcname = "Ilu2FreeR";

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
//				int ibs = _addrg[j];
//				pgold = _gl+ibs;
//				pgnew = _gl+_nznew;
//				for (int kii=0;kii<blaij;kii++) _gl[_nznew+kii] = _gl[ibs+kii];
//				for (int kii=0;kii<blaij;kii++) *pgnew++ = *pgold++;
//				pgold = _gu+ibs;
//				pgnew = _gu+_nznew;
//				for (int kii=0;kii<blaij;kii++) _gu[_nznew+kii] = _gu[ibs+kii];
//				for (kii=0;kii<blaij;kii++) *pgnew++ = *pgold++;
			};
//			_addrg[_nzjg] = _nznew;
			_addrg[_nzjg] = _addrg[j];
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
// CSMatrixCS: Store current row
//========================================================================================
void CSMatrixCS::Ilu2StoreRow (int _irow, double _tau1, // Store current row
								FILE *_mtrlfile, FILE *_mtrufile, dcmplx *_gread,
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
								dcmplx *_deleml, dcmplx *_delemu, int _nlist1, int *_lstupd, dcmplx *_fmaskl, dcmplx *_fmasku, 
								int *_ig, int *_jg, int &_nznew, int &_nzgblk, int *_addrg) const {

//	const char *funcname = "Ilu2StoreRow";

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	dcmplx *pg, *pfmask, *pd;

// Store current diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	_jg[_nzjg] = _irow;
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = _deleml+kiib;
		pg = _gread+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gl [_nznew+kjj*blai+kii] = _deleml[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};

	FPut (_mtrlfile,blai*blai,_gread,_nzgblk);

	_addrg[_nzjg] = _nzgblk;

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = _delemu+kiib;
		pg = _gread+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gu [_nznew+kjj*blai+kii] = _delemu[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};

	FPut (_mtrufile,blai*blai,_gread,_nzgblk);

	_addrg[_nzjg] = _nzgblk;

	_nzjg++;
	_nzju++;
	_nznew += blai * blai;
	_nzgblk += blai * blai;

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		blaij = blai*blaj;
		double aux = 0.0e0;
		pfmask = _fmaskl+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmaskl[jbsv*blai+kii];
			double auxloc = abs(*pfmask++);
			if (auxloc > aux) aux = auxloc;
		};
		pfmask = _fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = abs(*pfmask++);
			if (auxloc > aux) aux = auxloc;
		};
		_addrg[_nzjg] = _nzgblk;
		if (aux >= _tau1) {
			_jg[_nzjg] = j;
//			pg = _gl+_nznew;
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
//			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;

			FPut (_mtrlfile,blai*blaj,_fmaskl+jbsv*blai,_nzgblk);

//			pg = _gu+_nznew;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
//			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;

			FPut (_mtrufile,blai*blaj,_fmasku+jbsv*blai,_nzgblk);

			_nzjg++;
			_nzju++;
		} else {
			_jg[_nzjg] = -j;
//			pg = _gl+_nznew;
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
//			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;

			FPut (_mtrlfile,blai*blaj,_fmaskl+jbsv*blai,_nzgblk);

//			pg = _gu+_nznew;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
//			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;

			FPut (_mtrufile,blai*blaj,_fmasku+jbsv*blai,_nzgblk);

			_nzjg++;
			_nzjr++;
			_nzjrl++;
			_nzrtot += blaij;
		};
		_nznew += blaij;
		_nzgblk += blaij;
	};

	_ig[_irow+1] = _nzjg;

};

// Author: Kharchenko S.A.
// CSMatrixCS: Store the resulting U and scale it
//========================================================================================
void CSMatrixCS::Ich2StoreScaleU (int _nblks, int *_blks, FILE **_mtrlufiles, // Store the resulting U and scale it
									int *_sp2blk, int *_bsdia, int *_addrscl, int *_addrsclinv, int *_bsblk,
									int *_ig, int *_jg, int *_addrg,
									CSMatrixCS &_temp, 
									dcmplx *_aloc, dcmplx *_gread, dcmplx *_delem,
									double &_ops) const {

//	const char *funcname = "Ich2StoreScaleU";

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj, kkk, kiib, kjjb;
	int ibssc, jbs, jjbssc;
	int blai, blaj, blaij;
	dcmplx caux;

	dcmplx *pscloc;
	dcmplx *pgloc;
	dcmplx *paloc;

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
				blaij = blai*blaj;
				jjbssc = _bsdia[jjloc];

				FGet (_mtrlufiles[iblk],blai*blaj,_gread,jbs);

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
								pgloc = _gread+kjjb;
								pscloc = _delem+kiib;
								caux = czero;
								for (kkk=0;kkk<blai;kkk++) {
//									aux += _scl[ibssc+kii*blai+kkk] * _g[jbs+kjj*blai+kkk];
//									caux += *pscloc * *pgloc;
									caux.PlEq (*pscloc * *pgloc);
									pgloc++;
									pscloc++;
								};
//								_aloc[kjjb+kii] = caux;
								*paloc = caux;
								paloc += blai;
								kjjb += blai;
							};
							kiib += blai;
						};
						paloc = _aloc;
						pgloc = _temp.a+nz;
//						for (kii=0;kii<blaij;kii++) _temp.a[nz+kii] = _aloc[kii];
//						for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;

						FPut (_mtrlufiles[iblk],blai*blai,_aloc,nzgblk);

					} else {

						int ibs = _addrsclinv[jjloc];
						FGet (_mtrlufiles[jblk],blaj*blaj,_delem,ibs);

						for (kii=0;kii<blai;kii++) {
							kjjb = 0;
							paloc = _aloc+kii;
							for (kjj=0;kjj<blaj;kjj++) {
								caux = czero;
								pgloc = _gread+kii;
								pscloc = _delem+kjj;
								for (kkk=0;kkk<blaj;kkk++) {
//									aux += _g[jbs+kkk*blai+kii] * _sclinv[jjbs+kkk*blaj+kjj];
//									caux += *pgloc * *pscloc;
									caux.PlEq (*pgloc * *pscloc);
									pgloc += blai;
									pscloc += blaj;
								};
//								_aloc[kjjb+kii] = aux;
								*paloc = caux;
								paloc += blai;
								kjjb += blai;
							};
						};
						paloc = _aloc;
						pgloc = _temp.a+nz;
//						for (kii=0;kii<blaij;kii++) _temp.a[nz+kii] = _aloc[kii];
//						for (kii=0;kii<blaij;kii++) *pgloc++ = *paloc++;
						FPut (_mtrlufiles[iblk],blai*blaj,_aloc,nzgblk);
					};
					nzj++;
					nz += blaij;
					nzgblk += blaij;
					_ops += blaij*(blai+blaj);
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
// CSMatrixCS: Second order Incomplete LU decomposition
//========================================================================================
void CSMatrixCS::Ilu2 (ofstream &_fout, const CSlvParam _param, // Second order Incomplete LU decomposition
						int _nblks, int *_blks, 
						FILE **_mtralfiles, FILE **_mtraufiles, 
						const CSMatrixCS &_mtral, const CSMatrixCS &_mtrau,
						FILE **_mtrlfiles, FILE **_mtrufiles, 
						CSMatrixCS &_mtrl, CSMatrixCS &_mtru) { 

	const char *funcname = "Ilu2";

// Local variables

	int irow, i, j, jj, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzjgmx;
	int nz, nznew, nzju, nzjr, nzjrl, nzrtot, icheck, ibegfr;
	int blai, blaj, ibs;
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
	int icycleini, *imaskini;
	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	int *bsdia, *sp2blk, *bsblk;
	dcmplx *fmaskl, *fmasku, *diagg;
	int *addrscll, *addrsclu, *addrsclinvl, *addrsclinvu;
	dcmplx *gread, *lelem, *uelem, *deleml, *delemu;
	dcmplx *aloc, *uloc, *vloc, *work;
	double *dpiv, *eig, *rwork;

// Allocate and init array sp2blk

	sp2blk = new int [_mtral.nsupc];
	if (!sp2blk) MemoryFail (funcname);

	int iblk;
	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			sp2blk[i] = iblk;
		};
	};

// Create dummy matrix

	CSMatrixCS mtrdummy;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double complex arrays and compute scaling

	_mtrau.Ilu2AllocScale (_fout, _param.sclmin, 
							_nblks, sp2blk, _mtraufiles,
							fmaskl, fmasku, dpiv, 
							bsdia, diagg, 
							bsblk, _mtrlfiles, _mtrufiles,
							addrscll, addrsclu, addrsclinvl, addrsclinvu, 
							gread, lelem, uelem, deleml, delemu, 
							aloc, uloc, vloc, eig, 
							lwork, work, rwork, ops);

	imaskini = new int [_mtral.nsupr];
	if (!imaskini) MemoryFail (funcname);

	icycleini = -1;
	for (irow=0;irow<_mtral.nsupr;irow++) imaskini[irow] = -1;

// Allocate factorization arrays

	_mtral.Ilu2AllocG (_param.memory, nzjgmx, ig, jg, addrg);

// Main cycle over the rows

	int nzgblk;

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;
	nznew = 0;
	nzrtot = 0;

	ig[0] = 0;

	int *pia = _mtral.GetIa ();
	int *pja = _mtral.GetJa ();

	for (irow=0;irow<_mtral.nsupr;irow++) {

		iblk = sp2blk[irow];

		if (irow == _blks[iblk]) nzgblk = bsblk[iblk];

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			cout  << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			_fout << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;

// Output hystogram of the scaling values

//			int itype = 2;
			int nloc = _mtral.sprndr[irow];

			csvhyst (cout, "Ilu2Piv", nloc, dpiv);
			csvhyst (_fout,"Ilu2Piv", nloc, dpiv);

		};

// Init current scaled row

		mtrdummy.Ilu2InitRow (sp2blk, 
								_mtralfiles[iblk], _mtraufiles[iblk],
								_mtral, _mtrau, 
								_mtrlfiles, _mtrufiles,
								irow, nlist1, lstloc, icycle, imask, 
								bsdia, diagg, fmaskl, fmasku, gread,
								addrscll, addrsclu, 
								aloc, vloc, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			j = iv[irwprv];

			icolmn = jg[j];
			ibs = addrg[j];

			blai = _mtral.sprndr[irwprv+1]-_mtral.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _mtral.sprndr[icolmnl+1]-_mtral.sprndr[icolmnl];

			int jblk = sp2blk[irwprv];

//			for (kii=0;kii<blai*blaj;kii++) lelem[kii] = gl[ibs+kii];
//			for (kii=0;kii<blai*blaj;kii++) uelem[kii] = gu[ibs+kii];

			FGet(_mtrlfiles[jblk],blai*blaj,lelem,ibs);
			FGet(_mtrufiles[jblk],blai*blaj,uelem,ibs);

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			_mtral.Ilu2InitMask (elmisr, irow, irwprv, 
									iv, ig, jg, 
									nupd, lstupd, indlst, nlist1, lstloc, 
									icycle, imask, fmaskl, fmasku);

// Perform updates

			_mtral.Ich2UpdateRow (irow, 
									_mtrufiles[jblk], gread,
									irwprv, nupd, lstupd, indlst, 
									lelem, fmasku, addrg, iops);
			_mtral.Ich2UpdateRow (irow, 
									_mtrlfiles[jblk], gread,
									irwprv, nupd, lstupd, indlst, 
									uelem, fmaskl, addrg, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		Ich2UpdateTransp (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		_mtral.Ilu2FiltrRow (irow, _param.fcttyp, _param.theta, _param.tau2, 
								pia, pja, icycleini, imaskini,
								nlist1, lstloc, lstupd, 
								fmaskl, fmasku, bsdia, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin,
								irow, deleml, delemu, nlist1, lstupd, 
								fmaskl, fmasku, dpiv, 
								aloc, eig, uloc, vloc, lwork, work, rwork, iops);

// Count the local number of elements to be stored

		blai = _mtral.sprndr[irow+1]-_mtral.sprndr[irow];

		int nzloc=blai*blai;

		for (j=0;j<nlist1;j++) {
			jj = lstupd[j];
			blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
			nzloc += blai*blaj;
		};

// Check available memory

		if (nzjg+nlist1+1 > nzjgmx) {

// Free the memory containing R part that is not in use any more

			_mtral.Ilu2FreeR (_fout, irow, ibegfr, nzjg, nzju, nzjrl, 
								iv, ig, jg, 
								nznew, addrg);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx) {
			cout  << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
			_fout << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
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
								_mtrlfiles[iblk], _mtrufiles[iblk], gread,
								nzjg, nzju, nzjr, nzjrl, nzrtot,
								deleml, delemu, nlist1, lstupd, fmaskl, fmasku, 
								ig, jg, nznew, nzgblk, addrg);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = _mtral.n;

	csvhyst (cout, "Ilu2Piv", nloc, dpiv);
	csvhyst (_fout,"Ilu2Piv", nloc, dpiv);

// Free working arrays

	delete [] imask;
	delete [] lstloc;
	delete [] iv;
	delete [] madj;
	delete [] ibegm;
	delete [] lstupd;
	delete [] indlst;
	delete [] imaskini;

	delete [] fmaskl;
	delete [] fmasku;
	delete [] dpiv;
	delete [] diagg;

// Compute the size of the resulting matrix

	int nztot = 0;
	nzju = 0;

	for (i=0;i<_mtral.nsupr;i++) {
		blai = _mtral.sprndr[i+1]-_mtral.sprndr[i];
		for (j=ig[i];j<ig[i+1];j++) {
			jj = jg[j];
			if (jj >= 0) {
				blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
				nzju++;
				nztot += blai*blaj;
			};
		};
	};

// Perform the final filtering and scaling and compute matrix on return

	int nztotloc = 0;

	CSMatrixCS temp (_mtral.nsupr,nzju,nztotloc);

	_mtrl = temp;

	temp = mtrdummy;

	_mtru = _mtrl;

	_mtral.Ich2StoreScaleU (_nblks, _blks, _mtrlfiles,
							sp2blk, bsdia, addrscll, addrsclinvl, bsblk,
							ig, jg, addrg,
							_mtrl, 
							aloc, gread, deleml, 
							ops);
	_mtral.Ich2StoreScaleU (_nblks, _blks, _mtrufiles,
							sp2blk, bsdia, addrsclu, addrsclinvu, bsblk,
							ig, jg, addrg,
							_mtru, 
							aloc, gread, delemu, 
							ops);

// Free factorization arrays 

	blai = _mtral.sprndr[_mtral.nsupr]-_mtral.sprndr[_mtral.nsupr-1];

	int nzmema = 2*_mtral.nzatot-bsdia[_mtral.nsupr-1]-blai*blai;

	delete [] bsdia;
	delete [] bsblk;
	delete [] sp2blk;
	delete [] addrscll;
	delete [] addrsclu;
	delete [] addrsclinvl;
	delete [] addrsclinvu;
	delete [] lelem;
	delete [] deleml;
	delete [] delemu;
	delete [] aloc;
	delete [] vloc;
	delete [] uloc;
	delete [] eig;
	delete [] work;
	delete [] rwork;

	delete [] ig;
	delete [] jg;
	delete [] addrg;
	delete [] gread;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

// Output Ich2 factorization statistics

	nz = _mtral.nzatot;

	densu = 1.0e2 * (double) nztot / (double) nz;
	densr = 1.0e2 * (double) nzrtot / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	cout  << " Ilu2 preconditioner generation statistics: " << endl;
	cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	cout  << "     Costs = " << ops << " MvmA flops. " << endl;
	cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

	_fout << " Ilu2 preconditioner generation statistics: " << endl;
	_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
	_fout << "     Costs = " << ops << " MvmA flops. " << endl;
	_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif// 
