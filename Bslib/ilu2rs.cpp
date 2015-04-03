//------------------------------------------------------------------------------------------------
// File: ilu2rs.cpp
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
// CSMatrixRS: Allocate double data and init scaling for Ilu2
//========================================================================================
void CSMatrixRS::Ilu2AllocScale (ofstream &_fout, int _msglev, double _sclmin, // Allocate double data and init scaling for Ilu2
									double *&_fmaskl, double *&_fmasku, double *&_dpiv, 
									int *&_bsdia, double *&_diagg, 
									double *&_scll, double *&_sclu, double *&_sclinvl, double *&_sclinvu,
									double *&_lelem, double *&_uelem, double *&_deleml, double *&_delemu,
									double *&_aloc, double *&_uloc, double *&_vloc, 
									double *&_eig, int &_lwork, double *&_work,
									double &_ops) const {

	const char *funcname = "Ilu2AllocScale";

	_bsdia = new int [nsupr];
	if (!_bsdia) MemoryFail (funcname);

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

	_fmaskl = new double [n*blamx];
	if (!_fmaskl) MemoryFail (funcname);
	_fmasku = new double [n*blamx];
	if (!_fmasku) MemoryFail (funcname);

	_dpiv = new double [n];
	if (!_dpiv) MemoryFail (funcname);
	_diagg = new double [nz];
	if (!_diagg) MemoryFail (funcname);
	_scll = new double [nz];
	if (!_scll) MemoryFail (funcname);
	_sclu = new double [nz];
	if (!_sclu) MemoryFail (funcname);
	_sclinvl = new double [nz];
	if (!_sclinvl) MemoryFail (funcname);
	_sclinvu = new double [nz];
	if (!_sclinvu) MemoryFail (funcname);
	_lelem = new double [blamx*blamx];
	if (!_lelem) MemoryFail (funcname);
	_uelem = new double [blamx*blamx];
	if (!_uelem) MemoryFail (funcname);
	_deleml = new double [blamx*blamx];
	if (!_deleml) MemoryFail (funcname);
	_delemu = new double [blamx*blamx];
	if (!_delemu) MemoryFail (funcname);

// Allocate working arrays

	int info=0;

	_lwork = 10 * blamx;

	_aloc = new double [blamx*blamx];
	if (!_aloc) MemoryFail (funcname);
	_uloc = new double [blamx*blamx];
	if (!_uloc) MemoryFail (funcname);
	_vloc = new double [blamx*blamx];
	if (!_vloc) MemoryFail (funcname);
	_eig = new double [blamx];
	if (!_eig) MemoryFail (funcname);
	_work = new double [_lwork];
	if (!_work) MemoryFail (funcname);

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

// Compute its singular value decomposition

		int blal = blai;

		if (blal > 1) {
			dgesvd_ ("A", "A", &blal, &blal,
						_aloc, &blal, _eig, 
						_uloc, &blal, _vloc, &blal, 
						_work, &_lwork, &info);

			if (info != 0) throw " Error in the Lapack routine DGESVD";

			_ops += 8 * blai_2*blai;
		} else {
			_uloc[0] = 1.0e0;
			if (_aloc[0] >= 0.0e0) {
				_vloc[0] = 1.0e0;
				_eig[0] = _aloc[0];
			} else {
				_vloc[0] = -1.0e0;
				_eig[0] = -_aloc[0];
			};
		};

// Modify the pivots

		for (kii=0;kii<blai;kii++) {
			if (_eig[kii] < _sclmin) _eig[kii] = _sclmin;
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
			aux = _eig[kii];
			aux = 1.0e0 / aux;
			_eig[kii] = aux;
		};

		peig = _eig;
		kiib = 0;
		for (kii=0;kii<blai;kii++) {
			paloc = _uloc+kiib;
			pscl = _scll+ibsd+kii;
			for (kjj=0;kjj<blai;kjj++) {
//				_scll[ibsd+kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
				*pscl = *paloc++ * *peig;
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
				*pscl = *paloc * *peig;
				paloc += blai;
				pscl += blai;
			};
			peig++;
		};

		_ops += 4*blai_2 + 2*blai;

	};

// Reassign scaling factors if necessary

	if (false) {

		for (kii=0;kii<n*blamx;kii++) _scll[kii] = 0.0e0;
		for (kii=0;kii<n*blamx;kii++) _sclu[kii] = 0.0e0;
		for (kii=0;kii<n*blamx;kii++) _sclinvl[kii] = 0.0e0;
		for (kii=0;kii<n*blamx;kii++) _sclinvu[kii] = 0.0e0;

		for (isup=0;isup<nsupr;isup++) {
			blai = sprndr[isup+1]-sprndr[isup];
			ibsd = _bsdia[isup];
			for (kii=0;kii<blai;kii++) {
				_scll[ibsd+kii*blai+kii] = 1.0e0;
				_sclu[ibsd+kii*blai+kii] = 1.0e0;
				_sclinvl[ibsd+kii*blai+kii] = 1.0e0;
				_sclinvu[ibsd+kii*blai+kii] = 1.0e0;
			};
		};

	};

// Init working arrays

	pfmask = _fmaskl;
	for (kii=0;kii<n*blamx;kii++) *pfmask++ = 0.0e0;

	pfmask = _fmasku;
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
				pscl = _scll+ibssc+kii;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
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
				pscl = _sclu+ibssc+kjj;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _sclu[ibssc+kkk*blai+kjj];
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

	if (_msglev > 1) {
		csvhyst (cout, "ScDia",nloc, _dpiv);
	};
	if (_msglev >= 1) {
		csvhyst (_fout,"ScDia",nloc, _dpiv);
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Allocate G data for Ilu2
//========================================================================================
void CSMatrixRS::Ilu2AllocG (double _memory, int &_nzjgmx, int &_nzgmax, // Allocate G data for Ilu2
							int *&_ig, int *&_jg, int *&_addrg, double *&_gl, double *&_gu) const {

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
	_gl = new double [_nzgmax];
	if (!_gl) MemoryFail (funcname);
	_gu = new double [_nzgmax];
	if (!_gu) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixRS: Allocate G data for Ilu2
//========================================================================================
void CSMatrixRS::Ilu2AllocGDynamicByBlocks (int *&_ig, int *&_blks, int *&_sp2blk, // Allocate G data for Ilu2
											int **&_jg, int **&_addrg, double **&_gl, double **&_gu,
											int **&_jgblk, int **&_addrgblk, double **&_glblk, double **&_gublk) const {

	const char *funcname = "Ilu2AllocGDynamicByBlocks";

	_ig = new int [nsupr+1];
	if (!_ig) MemoryFail (funcname);
	_blks = new int [nsupr+1];
	if (!_blks) MemoryFail (funcname);
	_sp2blk = new int [nsupr+1];
	if (!_sp2blk) MemoryFail (funcname);
	_jg = new int * [nsupr+1];
	if (!_jg) MemoryFail (funcname);
	_addrg = new int * [nsupr+1];
	if (!_addrg) MemoryFail (funcname);
	_gl = new double * [nsupr+1];
	if (!_gl) MemoryFail (funcname);
	_gu = new double * [nsupr+1];
	if (!_gu) MemoryFail (funcname);
	_jgblk = new int * [nsupr+1];
	if (!_jgblk) MemoryFail (funcname);
	_addrgblk = new int * [nsupr+1];
	if (!_addrgblk) MemoryFail (funcname);
	_glblk = new double * [nsupr+1];
	if (!_glblk) MemoryFail (funcname);
	_gublk = new double * [nsupr+1];
	if (!_gublk) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CSMatrixRS: Init current scaled row
//========================================================================================
void CSMatrixRS::Ilu2InitRow (const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau, // Init current scaled row
								int _irow, int &_nlist1, 
								int *_lstloc, int _icycle, int *_imask, 
								int *_bsdia, double *_diagg, double *_fmaskl, double *_fmasku, 
								double *_scll, double *_sclu,
								double *_aloc, double *_vloc,
								int &_iops) const {

//	const char *funcname = "Ilu2InitRow";

	int kii, kjj, kkk, ibs, ibssc, jbssc;
	int kjjb, blai, blaj, ibsd, ibsv, jbsv, blaij, blai_2, blaj_2;
	double aux;

	double *paloc;
	double *pscloc;
	double *pfmask;
	double *pvloc;

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
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _scll[ibssc+kkk*blai+kii];
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
				pscloc = _sclu+jbssc+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _sclu[jbssc+kkk*blaj+kjj];
					aux += *paloc * *pscloc;
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
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += a[ibs+kjj*blai+kkk] * _sclu[ibssc+kkk*blai+kii];
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
				pscloc = _scll+jbssc+kjj;
				for (kkk=0;kkk<blaj;kkk++) {
//					aux += _aloc[kkk*blai+kii] * _scll[jbssc+kkk*blaj+kjj];
					aux += *paloc * *pscloc;
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
// CSMatrixRS: Init mask and update arrays
//========================================================================================
void CSMatrixRS::Ilu2InitMask (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const {

//	const char *funcname = "Ilu2InitMask";

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
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blai*blaj;kii++) _fmaskl[jbsv*blai+kii] = 0.0e0;
			for (kii=0;kii<blaij;kii++) *pfmask++ = 0.0e0;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blai*blaj;kii++) _fmasku[jbsv*blai+kii] = 0.0e0;
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
// CSMatrixRS: Init mask and update arrays
//========================================================================================
void CSMatrixRS::Ilu2InitMaskDynamicByBlocks (bool _elmisr, int _irow, int _irwprv, // Init mask and update arrays
							int *_iv, int *_ig, int *_jg, 
							int &_nupd, int *_lstupd, int *_indlst,
							int &_nlist1, int *_lstloc, 
							int _icycle, int *_imask, double *_fmaskl, double *_fmasku) const {

//	const char *funcname = "Ilu2InitMaskDynamicByBlocks";

	int kii;
	int blai, blaj, jbsv, blaij;
	double *pfmask;

	_nupd = 0;

	blai = sprndr[_irow+1]-sprndr[_irow];

	for (int k=_iv[_irwprv];k<_ig[_irwprv+1];k++) {

		int kloc = k-_ig[_irwprv];

		int jcolmn = _jg[kloc];

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
			for (kii=0;kii<blaij;kii++) *pfmask++ = 0.0e0;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blai*blaj;kii++) _fmasku[jbsv*blai+kii] = 0.0e0;
			for (kii=0;kii<blaij;kii++) *pfmask++ = 0.0e0;
			_nlist1++;
		};
		_lstupd[_nupd] = jcolmn;
		_indlst[_nupd] = kloc;
		_nupd++;

label1:;

	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Perform filtering of the current row
//========================================================================================
void CSMatrixRS::Ilu2FiltrRow (int _irow, double _theta, double _tau2, // Perform filtering of the current row
								int &_nlist1, int *_lstloc, int *_lstupd,
								double *_fmaskl, double* _fmasku, int *_bsdia, double *_diagg, int &_iops) const {

//	const char *funcname = "Ilu2FiltrRow";

	int nlist2 = 0;
	int kii, kjj, blai, blaj, ibsv, jbsv, jbsd, blaij;
	double *pfmaskl, *pfmasku, *pfmaskuu, *pdiag;

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsv = sprndr[_irow];

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
			double auxloc = *pfmaskl++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		pfmasku = _fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = *pfmasku++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};

		if (aux >= _tau2) {
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
					double auxloc = *pfmaskl;
					if (auxloc < 0.0e0) auxloc = -auxloc;
//					double auxloc1 = _fmasku[jbsv*blai+kjj*blai+kii];
					double auxloc1 = *pfmasku;
					if (auxloc1 < 0.0e0) auxloc1 = -auxloc1;
					if (auxloc1 > auxloc) auxloc = auxloc1;
					auxloc *= _theta;
//					_fmasku[ibsv*blai+kii*blai+kii] += auxloc;
					*pfmaskuu += auxloc;
//					_diagg[jbsd+kjj*blaj+kjj] += auxloc;
					*pdiag += auxloc;
					pfmaskl += blai;
					pfmasku += blai;
					pdiag += blaj+1;
				};
				pfmaskuu += blai+1;
			};
		};

	};

	_nlist1 = nlist2;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Compute diagonal pivot and scale current row
//========================================================================================
void CSMatrixRS::Ilu2PivScaleRow (double _pivmin, // Compute diagonal pivot and scale current row
									int _irow, double *_deleml, double *_delemu, 
									int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, double *_dpiv, 
									double *_aloc, double *_eig, double *_uloc, double *_vloc, int _lwork, double *_work,
									int &_iops) const {

//	const char *funcname = "Ilu2PivScaleRow";

	int kii, kjj, kkk, ibs;
	int blai, blaj, ibsv, jbsv;
	double aux;

// Copy diagonal supernode

	blai = sprndr[_irow+1]-sprndr[_irow];
	ibsv = sprndr[_irow];
	ibs = ibsv * blai;

	for (kii=0;kii<blai*blai;kii++) _aloc[kii] = _fmasku[ibs+kii];

// Compute singular value decomposition of the diagonal supernode

	int blal = blai;
	int info = 0;

	if (blal > 1) {

		dgesvd_ ("A", "A", &blal, &blal,
					_aloc, &blal, _eig, 
					_uloc, &blal, _vloc, &blal, 
					_work, &_lwork, &info);

		if (info != 0) throw " Error in the Lapack routine DGESVD";

		_iops += 8 * blai * blai * blai;

	} else {
		_uloc[0] = 1.0e0;
		if (_aloc[0] >= 0.0e0) {
			_vloc[0] = 1.0e0;
			_eig[0] = _aloc[0];
		} else {
			_vloc[0] = -1.0e0;
			_eig[0] = -_aloc[0];
		};
	};

// Modify singular values if necessary

	for (kii=0;kii<blai;kii++) {
		if (_eig[kii] < _pivmin) _eig[kii] = _pivmin;
		_eig[kii] = sqrt(_eig[kii]);
		_dpiv[ibsv+kii] = _eig[kii];
		_eig[kii] = 1.0e0/_eig[kii];
	};

// Compute diagonal supernodes

	int kiib;
	double *pdeloc;
	double *pfmloc;
	double *palloc;

	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pdeloc = _deleml+kii;
		palloc = _uloc+kiib;
		for (kjj=0;kjj<blai;kjj++) {
//			_deleml[kjj*blai+kii] = _uloc[kii*blai+kjj] * _eig[kii];
			*pdeloc = *palloc++ * _eig[kii];
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
			*pdeloc = *palloc * _eig[kii];
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
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _delemu[kkk*blai+kii] * _fmaskl[ibs+kjj*blai+kkk];
					aux += *pdeloc * *pfmloc++;
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = aux;
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
				aux = 0.0e0;
				for (kkk=0;kkk<blai;kkk++) {
//					aux += _deleml[kkk*blai+kii] * _fmasku[ibs+kjj*blai+kkk];
					aux += *pdeloc * *pfmloc++;
					pdeloc += blai;
				};
//				_aloc[kjj*blai+kii] = aux;
				*palloc = aux;
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
// CSMatrixRS: Free R part of the data
//========================================================================================
void CSMatrixRS::Ilu2FreeR (ofstream &_fout, int _msglev, int _irow, // Free R part of the data
							int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
							int *_iv, int *_ig, int *_jg, int &_nznew, 
							int *_addrg, double *_gl, double *_gu) const {

//	const char *funcname = "Ilu2FreeR";

	int iendfr = _irow;

	if (_msglev > 1) {
		cout  << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;
	};
	if (_msglev >= 1) {
		_fout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << endl;
	};

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

	if (_msglev > 1) {
		cout  << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;
	};
	if (_msglev >= 1) {
		_fout << " Nzju = " << _nzju << " Nzjr old = " << nzjro << " Nzjr new = " << _nzjrl << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Free R part of the data
//========================================================================================
void CSMatrixRS::Ilu2FreeRDynamicByBlocks (int _msglev, int _irow, // Free R part of the data
										int &_ibegfr, int &_nzjg, int _nzju, int &_nzjrl,
										int *_iv, int *_ig, int **_jg, int **_addrg, double **_gl, double **_gu,
										int _iblk, 
										int *_blks, int *_sp2blk, int **_jgblk, int **_addrgblk, double **_glblk, double **_gublk) const {

	const char *funcname = "Ilu2FreeRDynamicByBlocks";

	int iendfr = _irow;
	int ibegblk = _sp2blk[_ibegfr];

	if (_msglev > 1) {
		cout << " Free R data: Ibegfr = " << _ibegfr << " Iendfr = " << iendfr << " Ibegblk = " << ibegblk << " Iendblk = " << _iblk << endl;
	};

	_ibegfr = _blks[ibegblk];

	int ibgfrn = _ibegfr;

	_nzjg = _ig[_ibegfr];

	int ibeg = _ig[_ibegfr];

	int nzjro = _nzjrl;
	_nzjrl = 0;

	int i, j, jloc, iblkloc, ibegnd, iendnd, nzjblk, nzblk, ibegloc;
	int blai, blaj, blaij, nzjloc, nzloc, ibs, kii;
	int *jgloc;
	int *addrgloc;
	double *glloc;
	double *guloc;
	int *jgloc1;
	int *addrgloc1;
	double *glloc1;
	double *guloc1;
	int *jgarr;
	int *addrgarr;
	double *glarr;
	double *guarr;

	for (iblkloc=ibegblk;iblkloc<_iblk;iblkloc++) {

		ibegnd = _blks[iblkloc];
		iendnd = _blks[iblkloc+1];

// Initial scan of the current block

		nzjblk = 0;
		nzblk = 0;

		for (i=ibegnd;i<iendnd;i++) {
			blai = sprndr[i+1]-sprndr[i];
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
				blaj = sprndr[icolmn+1]-sprndr[icolmn];
				nzjblk++;
				nzblk += blai*blaj;
label:;
			};

		};

// Allocate new memory and store data of the block row by row

		jgarr = new int [nzjblk];
		if (!jgarr) MemoryFail (funcname);
		addrgarr = new int [nzjblk];
		if (!addrgarr) MemoryFail (funcname);
		glarr = new double [nzblk];
		if (!glarr) MemoryFail (funcname);
		guarr = new double [nzblk];
		if (!guarr) MemoryFail (funcname);

		nzjblk = 0;
		nzblk = 0;

		for (i=ibegnd;i<iendnd;i++) {

			blai = sprndr[i+1]-sprndr[i];

// Initial scan of the current row

			jgloc1 = _jg[i];
			addrgloc1 = _addrg[i];
			glloc1 = _gl[i];
			guloc1 = _gu[i];

			jgloc = jgarr+nzjblk;
			addrgloc = addrgarr+nzjblk;
			glloc = glarr+nzblk;
			guloc = guarr+nzblk;

			nzjloc = 0;
			nzloc = 0;

			for (j=ibeg;j<_ig[i+1];j++) {

				jloc = j-ibeg;

				int icolmn = jgloc1[jloc];

				if (icolmn < 0) {
					icolmn = -icolmn;
					if (icolmn < iendfr) goto label1;
					_nzjrl++;
				};

				blaj = sprndr[icolmn+1]-sprndr[icolmn];

				jgloc[nzjloc] = jgloc1[jloc];
				addrgloc[nzjloc] = nzloc;
				ibs = addrgloc1[jloc];

				blaij = blai*blaj;

				for (kii=0;kii<blaij;kii++) glloc[nzloc+kii] = glloc1[ibs+kii];
				for (kii=0;kii<blaij;kii++) guloc[nzloc+kii] = guloc1[ibs+kii];

				nzjloc++;
				nzloc += blaij;

				_nzjg++;

label1:;

			};

			_jg[i] = jgloc;
			_addrg[i] = addrgloc;
			_gl[i] = glloc;
			_gu[i] = guloc;

// Final scan of the current row

			ibeg = _ig[i+1];

			_ig[i+1] = _nzjg;

			if (_nzjrl == 0 && ibgfrn == i) ibgfrn = i+1;

// Update iv

			_iv[i] += _nzjg-ibeg;

			nzjblk += nzjloc;
			nzblk += nzloc;

		};

// Free previous memory and store the new one

		delete [] _jgblk[iblkloc];
		delete [] _addrgblk[iblkloc];
		delete [] _glblk[iblkloc];
		delete [] _gublk[iblkloc];

		_jgblk[iblkloc] = jgarr;
		_addrgblk[iblkloc] = addrgarr;
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
// CSMatrixRS: Store current row
//========================================================================================
void CSMatrixRS::Ilu2StoreRow (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
								double *_deleml, double *_delemu, int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, 
								int *_ig, int *_jg, int &_nznew, int *_addrg, double *_gl, double *_gu) const {

//	const char *funcname = "Ilu2StoreRow";

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	double *pg, *pfmask, *pd;

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

	for (int i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		blaij = blai*blaj;
		double aux = 0.0e0;
		pfmask = _fmaskl+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmaskl[jbsv*blai+kii];
			double auxloc = *pfmask++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		pfmask = _fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = *pfmask++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		_addrg[_nzjg] = _nznew;
		if (aux >= _tau1) {
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
// CSMatrixRS: Store current row
//========================================================================================
void CSMatrixRS::Ilu2StoreRowDynamicByBlocks (int _irow, double _tau1, // Store current row
								int &_nzjg, int &_nzju, int &_nzjr, int &_nzjrl, int &_nzrtot,
								double *_deleml, double *_delemu, int _nlist1, int *_lstupd, double *_fmaskl, double *_fmasku, 
								int *_ig, int **_jg, int **_addrg, double **_gl, double **_gu,
								int &_iblk, int _nzjallocmax, int &_nzjalloc, int &_nzjused, int _nzallocmax, int &_nzalloc, int &_nzused,
								int *_blks, int *_sp2blk, int **_addrgblk, int **_jgblk, 
								double **_glblk, double **_gublk) const {

	const char *funcname = "Ilu2StoreRowDynamicByBlocks";

	int kii, kjj, kiib;
	int blai, blaj, jbsv, blaij;
	double *pg, *pfmask, *pd;

// Count the size of memory to be used

	int nzjloc, nzloc;

	nzjloc = _nlist1+1;
	blai = sprndr[_irow+1]-sprndr[_irow];

	nzloc = blai*blai;

	int i, j;
	for (i=0;i<_nlist1;i++) {
		j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		nzloc += blai*blaj;
	};

// Allocate the memory 

	int *jgloc, *addrgloc;
	double *glloc, *guloc;

	if (nzjloc > _nzjalloc-_nzjused || nzloc > _nzalloc-_nzused) {

		_nzjalloc = _nzjallocmax;
		if (_nzjalloc < nzjloc) _nzjalloc = nzjloc;
		_nzalloc = _nzallocmax;
		if (_nzalloc < nzloc) _nzalloc = nzloc;

		jgloc = new int [_nzjalloc];
		if (!jgloc) MemoryFail (funcname);
		addrgloc = new int [_nzjalloc];
		if (!addrgloc) MemoryFail (funcname);
		glloc = new double [_nzalloc];
		if (!glloc) MemoryFail (funcname);
		guloc = new double [_nzalloc];
		if (!guloc) MemoryFail (funcname);

		_iblk++;
		_jgblk[_iblk] = jgloc;
		_addrgblk[_iblk] = addrgloc;
		_glblk[_iblk] = glloc;
		_gublk[_iblk] = guloc;
		_nzjused = 0;
		_nzused = 0;

	} else {
		jgloc = _jgblk[_iblk]+_nzjused;
		addrgloc = _addrgblk[_iblk]+_nzjused;
		glloc = _glblk[_iblk]+_nzused;
		guloc = _gublk[_iblk]+_nzused;
	};

// Store current diagonal supernode

	nzjloc = 0;
	nzloc = 0;

	jgloc[nzjloc] = _irow;
	kiib = 0;
	for (kii=0;kii<blai;kii++) {
		pd = _deleml+kiib;
		pg = glloc+nzloc+kii;
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
		pg = guloc+nzloc+kii;
		for (kjj=0;kjj<blai;kjj++) {
//			_gu [_nznew+kjj*blai+kii] = _delemu[kii*blai+kjj];
			*pg = *pd++;
			pg += blai;
		};
		kiib += blai;
	};
	addrgloc[nzjloc] = nzloc;
	_nzjg++;
	_nzju++;
	nzjloc++;
	nzloc += blai*blai;

	for (i=0;i<_nlist1;i++) {
		int j = _lstupd[i];
		blaj = sprndr[j+1]-sprndr[j];
		jbsv = sprndr[j];
		blaij = blai*blaj;
		double aux = 0.0e0;
		pfmask = _fmaskl+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmaskl[jbsv*blai+kii];
			double auxloc = *pfmask++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		pfmask = _fmasku+jbsv*blai;
		for (kii=0;kii<blaij;kii++) {
//			double auxloc = _fmasku[jbsv*blai+kii];
			double auxloc = *pfmask++;
			if (auxloc < 0.0e0) auxloc = -auxloc;
			if (auxloc > aux) aux = auxloc;
		};
		addrgloc[nzjloc] = nzloc;
		if (aux >= _tau1) {
			jgloc[nzjloc] = j;
			pg = glloc+nzloc;
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			pg = guloc+nzloc;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			_nzjg++;
			_nzju++;
		} else {
			jgloc[nzjloc] = -j;
			pg = glloc+nzloc;
			pfmask = _fmaskl+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gl [_nznew+kii] = _fmaskl[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			pg = guloc+nzloc;
			pfmask = _fmasku+jbsv*blai;
//			for (kii=0;kii<blaij;kii++) _gu [_nznew+kii] = _fmasku[jbsv*blai+kii];
			for (kii=0;kii<blaij;kii++) *pg++ = *pfmask++;
			_nzjg++;
			_nzjr++;
			_nzjrl++;
			_nzrtot += blaij;
		};
		nzjloc++;
		nzloc += blaij;
	};

	_ig[_irow+1] = _nzjg;

	_jg[_irow] = jgloc;
	_addrg[_irow] = addrgloc;
	_gl[_irow] = glloc;
	_gu[_irow] = guloc;

	_blks[_iblk+1] = _irow+1;
	_sp2blk[_irow] = _iblk;
	_nzjused += nzjloc;
	_nzused += nzloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Second order Incomplete LU decomposition
//========================================================================================
void CSMatrixRS::Ilu2 (ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition
						const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
						CSMatrixRS &_mtrl, CSMatrixRS &_mtru) { 

//	const char *funcname = "Ilu2";

// Local variables

	int irow, i, j, jj, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg, nzgmax, nzjgmx;
	int nz, nznew, nzju, nzjr, nzjrl, nzrtot, icheck, msglev, ibegfr;
	int kii, blai, blaj, ibs;
	int *ig, *jg, *addrg;
	double *gl, *gu;
	bool elmisr;

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;
	int iops;

	icheck = _param.ichkfct;
	msglev = (int) _param.msglev;

	ops = 0.0e0;

// Init time measurement

	time0 = clock ();

// Define working arrays

	int lwork;
	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	int *bsdia;
	double *fmaskl, *fmasku, *dpiv, *diagg, *scll, *sclu, *sclinvl, *sclinvu;
	double *lelem, *uelem, *deleml, *delemu;
	double *aloc, *uloc, *vloc, *eig, *work;

// Create dummy matrix

	CSMatrixRS mtrdummy;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	_mtrau.Ilu2AllocScale (_fout, msglev, _param.sclmin, 
							fmaskl, fmasku, dpiv, bsdia, diagg, 
							scll, sclu, sclinvl, sclinvu, 
							lelem, uelem, deleml, delemu, 
							aloc, uloc, vloc, eig, 
							lwork, work, ops);

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

	for (irow=0;irow<_mtral.nsupr;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			if (msglev > 1) {
				cout  << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			};
			if (msglev >= 1) {
				_fout << " Ilu2: Irow = " << irow << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << " Nzjgmx = " << nzjgmx << endl;
			};

// Output hystogram of the scaling values

//			int itype = 2;
			int nloc = _mtral.sprndr[irow];

			if (msglev > 1) {
				csvhyst (cout, "Ilu2Piv", nloc, dpiv);
			};
			if (msglev >= 1) {
				csvhyst (_fout,"Ilu2Piv", nloc, dpiv);
			};

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

			blai = _mtral.sprndr[irwprv+1]-_mtral.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _mtral.sprndr[icolmnl+1]-_mtral.sprndr[icolmnl];

			for (kii=0;kii<blai*blaj;kii++) lelem[kii] = gl[ibs+kii];
			for (kii=0;kii<blai*blaj;kii++) uelem[kii] = gu[ibs+kii];

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

		_mtral.Ilu2FiltrRow (irow, _param.theta, _param.tau2, 
								nlist1, lstloc, lstupd, 
								fmaskl, fmasku, bsdia, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin,
								irow, deleml, delemu, nlist1, lstupd, 
								fmaskl, fmasku, dpiv, 
								aloc, eig, uloc, vloc, lwork, work, iops);

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

			_mtral.Ilu2FreeR (_fout, msglev, irow, ibegfr, nzjg, nzju, nzjrl, 
								iv, ig, jg, 
								nznew, addrg, gl, gu);

		};

// Check available memory one more time

		if (nzjg+nlist1+1 > nzjgmx || nznew+nzloc > nzgmax) {
			if (msglev > 1) {
				cout  << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
			};
			if (msglev >= 1) {
				_fout << " Insufficient workspace in Ilu2, Irow = " << irow << endl;
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
								nzjg, nzju, nzjr, nzjrl, nzrtot,
								deleml, delemu, nlist1, lstupd, fmaskl, fmasku, 
								ig, jg, nznew, addrg, gl, gu);

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = _mtral.n;

	if (msglev > 1) {
		csvhyst (cout, "Ilu2Piv", nloc, dpiv);
	};
	if (msglev >= 1) {
		csvhyst (_fout,"Ilu2Piv", nloc, dpiv);
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

	CSMatrixRS temp (_mtral.nsupr,nzju,nztot);

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
	delete [] deleml;
	delete [] delemu;
	delete [] aloc;
	delete [] vloc;
	delete [] uloc;
	delete [] eig;
	delete [] work;

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

	nz = _mtral.nza;

	densu = 1.0e2 * (double) nztot / (double) nz;
	densr = 1.0e2 * (double) nzrtot / (double) nz;

	_param.density  = (double) nztot / (double) nz;
	_param.density2 = (double) nzrtot / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	if (msglev > 1) {
		cout  << " Ilu2 preconditioner generation statistics: " << endl;
		cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (msglev >= 1) {
		_fout << " Ilu2 preconditioner generation statistics: " << endl;
		_fout << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		_fout << "     Costs = " << ops << " MvmA flops. " << endl;
		_fout << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Second order Incomplete LU decomposition with dynamic memory allocation by blocks
//========================================================================================
void CSMatrixRS::Ilu2DynamicByBlocks (ofstream &_fout, CSlvParam &_param, // Second order Incomplete LU decomposition with dynamic memory allocation by blocks
										const CSMatrixRS &_mtral, const CSMatrixRS &_mtrau,
										CSMatrixRS &_mtrl, CSMatrixRS &_mtru) { 

//	const char *funcname = "Ilu2DynamicByBlocks";

// Local variables

	int irow, i, j, jj, nlist1, irwprv, icolmn, nupd;
	int icycle, nzjg, nzg;
	int nz, nzju, nzjr, nzjrl, nzrtot, icheck, msglev, ibegfr;
	int kii, blai, blaj, ibs;
	int nblks, iblk, nzjallocmax, nzjalloc, nzjused, nzallocmax, nzalloc, nzused;
	int *ig, *blks, *sp2blk;
	int **jg, **addrg;
	double **gl, **gu;
	int **jgblk, **addrgblk;
	double **glblk, **gublk;
	bool elmisr;

// Statistics data

	double ops, densu, densr, perf, tottim;
	clock_t time0, time1;
	clock_t time2, time3;
	int iops;

	icheck = _param.ichkfct;
	msglev = (int) _param.msglev;

	ops = 0.0e0;

// Init time measurement

	time0 = clock ();

// Define working arrays

	int lwork;
	int *imask, *lstloc, *iv, *madj, *ibegm, *lstupd, *indlst;
	int *bsdia;
	double *fmaskl, *fmasku, *dpiv, *diagg, *scll, *sclu, *sclinvl, *sclinvu;
	double *lelem, *uelem, *deleml, *delemu;
	double *aloc, *uloc, *vloc, *eig, *work;

// Create dummy matrix

	CSMatrixRS mtrdummy;

// Allocate and init working arrays including the scaling

// Integer arrays

	_mtral.Ich2AllocInt (icycle, imask, lstloc, iv, madj, ibegm, lstupd, indlst);

// Double arrays and compute scaling

	_mtrau.Ilu2AllocScale (_fout, msglev, _param.sclmin, 
							fmaskl, fmasku, dpiv, bsdia, diagg, 
							scll, sclu, sclinvl, sclinvu, 
							lelem, uelem, deleml, delemu, 
							aloc, uloc, vloc, eig, 
							lwork, work, ops);

// Allocate factorization arrays

	_mtral.Ilu2AllocGDynamicByBlocks (ig, blks, sp2blk,
										jg, addrg, gl, gu,
										jgblk, addrgblk, glblk, gublk);

// Set free R parameters

	int nfreer = (int) _param.nfreer;

	int nrowsfree;

	if (nfreer <= 0) {
		nrowsfree = 0;
	} else {
		nrowsfree = _mtral.nsupr / nfreer;
	};

// Main cycle over the rows

//	int nrowstime = 10000;

	int jloc;
	int *jgloc, *addrgloc;
	double *glloc, *guloc;

	ibegfr = 0;

	nzjg = 0;
	nzg = 0;
	nzju = 0;
	nzjr = 0;
	nzjrl = 0;
	nzrtot = 0;

	iblk = -1;
	nzjalloc = 0;
	nzjallocmax = (int) ( (double)_mtral.nzja * (-_param.memory) );
	if (nzjallocmax < 0) nzjallocmax = 0;
	nzalloc = 0;
	nzallocmax = (int) ( (double)_mtral.nza * (-_param.memory) );
	if (nzallocmax < 0) nzallocmax = 0;
	nzjused = 0;
	nzused = 0;

	ig[0] = 0;
	blks[0] = 0;

	time2 = clock ();
	time3 = clock ();

	for (irow=0;irow<_mtral.nsupr;irow++) {

		icycle++;
		iops = 0;

// Output statistics if necessary

		if (irow > 0 && irow%icheck==0) {

			if (msglev > 1) {
				cout  << " Ilu2: Irow = " << irow << " Iblk = " << iblk << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << endl;
			};
			if (msglev >= 1) {
				_fout << " Ilu2: Irow = " << irow << " Iblk = " << iblk << " Nzjg = " << nzjg << " Nzjr = " << nzjrl << endl;
			};

// Output hystogram of the scaling values

//			int itype = 2;
			int nloc = _mtral.sprndr[irow];

			if (msglev > 1) {
				csvhyst (cout, "Ilu2Piv", nloc, dpiv);
			};
			if (msglev >= 1) {
				csvhyst (_fout,"Ilu2Piv", nloc, dpiv);
			};

		};

#ifdef __FlowVisionOut
		if (irow%nrowstime == 0) {

			time3 = clock ();

			tottim = (double) (time3-time2) / (double) CLOCKS_PER_SEC;

			if (tottim > _param.timefctout) {

				nz = _mtral.nza;

				densu = (double) nzju / (double) nz;
				densr = (double) nzjr / (double) nz;

				int prcent = 100*irow / _mtral.nsupr;

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
								bsdia, diagg, fmaskl, fmasku, 
								scll, sclu, 
								aloc, vloc, iops);

// Update current row by the previous ones

		irwprv = ibegm[irow];

		while (irwprv != -1) {

			jgloc = jg[irwprv];
			addrgloc = addrg[irwprv];
			glloc = gl[irwprv];
			guloc = gu[irwprv];

			j = iv[irwprv];

			jloc = j-ig[irwprv];

			icolmn = jgloc[jloc];
			ibs = addrgloc[jloc];

			blai = _mtral.sprndr[irwprv+1]-_mtral.sprndr[irwprv];
			int icolmnl = icolmn;
			if (icolmnl < 0) icolmnl = -icolmnl;
			blaj = _mtral.sprndr[icolmnl+1]-_mtral.sprndr[icolmnl];

			for (kii=0;kii<blai*blaj;kii++) lelem[kii] = glloc[ibs+kii];
			for (kii=0;kii<blai*blaj;kii++) uelem[kii] = guloc[ibs+kii];

			elmisr = icolmn < 0;

// Init and register if necessary the mask

			_mtral.Ilu2InitMaskDynamicByBlocks (elmisr, irow, irwprv, 
									iv, ig, jgloc, 
									nupd, lstupd, indlst, nlist1, lstloc, 
									icycle, imask, fmaskl, fmasku);

// Perform updates

			_mtral.Ich2UpdateRow (irow, irwprv, nupd, lstupd, indlst, 
									lelem, fmasku, addrgloc, guloc, iops);
			_mtral.Ich2UpdateRow (irow, irwprv, nupd, lstupd, indlst, 
									uelem, fmaskl, addrgloc, glloc, iops);

			irwprv = madj[irwprv];

		};

// Update arrays iv, madj and ibegm for the next column

		_mtral.Ich2UpdateTranspDynamic (irow, ibegm, madj, iv, ig, jg);

// Perform filtering of the row with diagonal modification if necessary 

		_mtral.Ilu2FiltrRow (irow, _param.theta, _param.tau2, 
								nlist1, lstloc, lstupd, 
								fmaskl, fmasku, bsdia, diagg, iops);

// Compute and store the pivot value and scale current row

		_mtral.Ilu2PivScaleRow (_param.pivmin,
								irow, deleml, delemu, nlist1, lstupd, 
								fmaskl, fmasku, dpiv, 
								aloc, eig, uloc, vloc, lwork, work, iops);

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
								nzjg, nzju, nzjr, nzjrl, nzrtot,
								deleml, delemu, nlist1, lstupd, fmaskl, fmasku, 
								ig, jg, addrg, gl, gu,
								iblk, nzjallocmax, nzjalloc, nzjused, nzallocmax, nzalloc, nzused,
								blks, sp2blk, addrgblk, jgblk, glblk, gublk);

// Free R part if necessary

		if (nrowsfree > 0 && irow > 0 && irow % nrowsfree == 0 && irow != _mtral.n-1) {

			int msglev = (int) _param.msglev;

			_mtral.Ilu2FreeRDynamicByBlocks (msglev, irow+1,
										ibegfr, nzjg, nzju, nzjrl,
										iv, ig, jg, addrg, gl, gu,
										iblk, blks, sp2blk, jgblk, addrgblk, glblk, gublk);

		};

// End of main cycle over the rows

		ops += iops;

	};

// Output hystogram of the scaling values

//	int itype = 2;
	int nloc = _mtral.n;

	if (msglev > 1) {
		csvhyst (cout, "Ilu2Piv", nloc, dpiv);
	};
	if (msglev >= 1) {
		csvhyst (_fout,"Ilu2Piv", nloc, dpiv);
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

// Compute the size of the resulting matrix

	int nztot = 0;
	nzju = 0;

	for (i=0;i<_mtral.nsupr;i++) {
		blai = _mtral.sprndr[i+1]-_mtral.sprndr[i];
		jgloc = jg[i];
		for (j=ig[i];j<ig[i+1];j++) {
			jloc = j-ig[i];
			jj = jgloc[jloc];
			if (jj >= 0) {
				blaj = _mtral.sprndr[jj+1]-_mtral.sprndr[jj];
				nzju++;
				nztot += blai*blaj;
			};
		};
	};

// Perform the final filtering and scaling and compute matrix on return

	nblks = iblk+1;

	CSMatrixRS temp (_mtral.nsupr,nzju,nztot);

	_mtrl = temp;

	temp = mtrdummy;

	_mtru = _mtrl;

	_mtral.Ich2StoreScaleUDynamic (bsdia, scll, sclinvl, 
							ig, jg, addrg, gl, 
							_mtrl, 
							aloc, ops);

	for (iblk=0;iblk<nblks;iblk++) {
		delete [] glblk[iblk];
	};

	_mtral.Ich2StoreScaleUDynamic (bsdia, sclu, sclinvu, 
							ig, jg, addrg, gu, 
							_mtru, 
							aloc, ops);

	for (iblk=0;iblk<nblks;iblk++) {
		delete [] gublk[iblk];
	};

// Free factorization arrays 

	blai = _mtral.sprndr[_mtral.nsupr]-_mtral.sprndr[_mtral.nsupr-1];

	int nzmema = 2*_mtral.nza-bsdia[_mtral.nsupr-1]-blai*blai;

	delete [] bsdia;
	delete [] scll;
	delete [] sclu;
	delete [] sclinvl;
	delete [] sclinvu;
	delete [] lelem;
	delete [] deleml;
	delete [] delemu;
	delete [] aloc;
	delete [] vloc;
	delete [] uloc;
	delete [] eig;
	delete [] work;

	for (iblk=0;iblk<nblks;iblk++) {
		delete [] jgblk[iblk];
		delete [] addrgblk[iblk];
	};

	delete [] ig;
	delete [] blks;
	delete [] sp2blk;
	delete [] jg;
	delete [] addrg;
	delete [] gl;
	delete [] gu;
	delete [] jgblk;
	delete [] addrgblk;
	delete [] glblk;
	delete [] gublk;

// Finalize time measurement

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	if (tottim <= 0.0e0) tottim = 1.0e0;

	_param.timefct = tottim;

// Output Ich2 factorization statistics

	nz = _mtral.nza;

	densu = 1.0e2 * (double) nztot / (double) nz;
	densr = 1.0e2 * (double) nzrtot / (double) nz;

	_param.density  = (double) nztot / (double) nz;
	_param.density2 = (double) nzrtot / (double) nz;

	perf = 1.0e-6 * ops / tottim;

	ops = ops / (double) nzmema;

	if (msglev > 1) {
		cout  << " Ilu2 preconditioner generation statistics: " << endl;
		cout  << "     DensU = " << densu << " %   DensR = " << densr << " % " << endl;
		cout  << "     Costs = " << ops << " MvmA flops. " << endl;
		cout  << "     Time  = " << tottim << " sec.   Perf = " << perf << " Mflops." << endl;
	};

	if (msglev >= 1) {
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
