//------------------------------------------------------------------------------------------------
// File: globals.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include "globals.h"
#include "hyst.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iosfwd>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <time.h>

using namespace std;

// Author: Kharchenko S.A.
// Description: Dummy integer function
//========================================================================================
void IntCompDummy (int &_n) {

};

static clock_t gltimes[100];
static double  glwtimes[100];
static int glcurTimerIndex = -1;

// Author: Alexander Dyadkin, Kharchenko S.A.
// StartTimer()
//========================================================================================
void StartTimer()
{
   glcurTimerIndex++;
   assert(glcurTimerIndex>=0 && glcurTimerIndex<100);
   //
   gltimes[glcurTimerIndex]  = clock();
   glwtimes[glcurTimerIndex] = CMPIExchange::WallTimeMPI();
}

// Author: Alexander Dyadkin, Kharchenko S.A.
// StopTimer()
//========================================================================================
void StopTimer(double& _time0, double& _time1)
{
   assert(glcurTimerIndex>=0 && glcurTimerIndex<100);
   //
   clock_t time1;
   double  wtime1;
   time1 = clock();
   wtime1 = CMPIExchange::WallTimeMPI();
   _time0 = (double) (time1-gltimes[glcurTimerIndex]) / (double) CLOCKS_PER_SEC;
   _time1 = wtime1-glwtimes[glcurTimerIndex];
   //
   glcurTimerIndex--;
}

// Author: Alexander Dyadkin, Kharchenko S.A.
// StopTimer()
//========================================================================================
void StopTimer(const char* _str)
{
   assert(glcurTimerIndex>=0 && glcurTimerIndex<100);
   //
   clock_t time1;
   double  wtime1;
   time1 = clock();
   wtime1 = CMPIExchange::WallTimeMPI();
   double tottim;
   tottim = (double) (time1-gltimes[glcurTimerIndex]) / (double) CLOCKS_PER_SEC;
   std::cout << _str << ' ' << tottim << ' ' << wtime1-glwtimes[glcurTimerIndex] << std::endl;
   glcurTimerIndex--;
}

// Author: Kharchenko S.A.
// Description: Solve system with dense matrix by SVD
// SolveDenseBySVD()
//========================================================================================
// Solve system with dense matrix by SVD
void SolveDenseBySVD (int _n, int _nrhs, // Solve system with dense matrix by SVD
								dcmplx *_rmatr, int _ldrmatr,
								dcmplx *_rhs, int _ldrhs,
								dcmplx *_sol, int _ldsol) {

	const char *funcname = "SolveDenseBySVD";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ntot = _n;
	int lworkloc = 10*ntot;
	int info = 0;

	dcmplx *aloc;
	dcmplx *uloc;
	dcmplx *vloc;
	double *svloc;
	dcmplx *workloc;
	double *rworkloc;

	aloc = new dcmplx [ntot*ntot];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [ntot*ntot];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [ntot*ntot];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [ntot];
	if (!svloc) MemoryFail (funcname);
	workloc = new dcmplx [lworkloc];
	if (!workloc) MemoryFail (funcname);
	rworkloc = new double [lworkloc];
	if (!rworkloc) MemoryFail (funcname);

// Compute SVD

	int kii, kjj, irhs;

	for (kii=0;kii<ntot*ntot;kii++) aloc[kii] = czero;

	for (kjj=0;kjj<ntot;kjj++) {
		for (kii=0;kii<ntot;kii++) {
			aloc[kjj*ntot+kii] = _rmatr[kjj*_ldrmatr+kii];
		};
	};

	zgesvd_ ("A", "A", &ntot, &ntot, 
				aloc, &ntot, svloc, uloc, &ntot, vloc, &ntot, 
				workloc, &lworkloc, rworkloc, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
		throw " Error inside Lapack routine ZGESVD";
	};

	CHyst hyst;

	hyst.UpdateHystogram (ntot, svloc);

	cout << " Singular values: " << hyst << endl;

// Solve

	dcmplx caux;
	double daux;

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (kjj=0;kjj<ntot;kjj++) {
			caux = czero;
			for (kii=0;kii<ntot;kii++) {
				caux += conj(uloc[kjj*ntot+kii]) * _rhs[irhs*_ldrhs+kii];
			};
			aloc[irhs*ntot+kjj] = caux;
		};
	};

	for (kjj=0;kjj<ntot;kjj++) {
		daux = svloc[kjj];
		daux = 1.0e0 / daux;
		for (irhs=0;irhs<_nrhs;irhs++) {
			caux = aloc[irhs*ntot+kjj] * daux;
			aloc[irhs*ntot+kjj] = caux;
		};
	};

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (kjj=0;kjj<ntot;kjj++) {
			caux = czero;
			for (kii=0;kii<ntot;kii++) {
				caux += conj(vloc[kjj*ntot+kii]) * aloc[irhs*ntot+kii];
			};
			_sol[irhs*_ldsol+kjj] = caux;
		};
	};

// Free work arrays

	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] workloc;
	delete [] rworkloc;

};

// Author: Kharchenko S.A.
// Description: Solve system with dense matrix by SVD
// SolveDenseBySVDTikhonov()
//========================================================================================
// Solve system with dense matrix by SVD
void SolveDenseBySVDTikhonov (double _alpha, int _n, int _nrhs, // Solve system with dense matrix by SVD
								dcmplx *_rmatr, int _ldrmatr,
								dcmplx *_rhs, int _ldrhs,
								dcmplx *_sol, int _ldsol) {

	const char *funcname = "SolveDenseBySVDTikhonov";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ntot = _n;
	int lworkloc = 10*ntot;
	int info = 0;

	dcmplx *aloc;
	dcmplx *uloc;
	dcmplx *vloc;
	double *svloc;
	dcmplx *workloc;
	double *rworkloc;

	aloc = new dcmplx [ntot*ntot];
	if (!aloc) MemoryFail (funcname);
	uloc = new dcmplx [ntot*ntot];
	if (!uloc) MemoryFail (funcname);
	vloc = new dcmplx [ntot*ntot];
	if (!vloc) MemoryFail (funcname);
	svloc = new double [ntot];
	if (!svloc) MemoryFail (funcname);
	workloc = new dcmplx [lworkloc];
	if (!workloc) MemoryFail (funcname);
	rworkloc = new double [lworkloc];
	if (!rworkloc) MemoryFail (funcname);

// Compute SVD

	int kii, kjj, irhs;

	for (kii=0;kii<ntot*ntot;kii++) aloc[kii] = czero;

	for (kjj=0;kjj<ntot;kjj++) {
		for (kii=0;kii<ntot;kii++) {
			aloc[kjj*ntot+kii] = _rmatr[kjj*_ldrmatr+kii];
		};
	};

	zgesvd_ ("A", "A", &ntot, &ntot, 
				aloc, &ntot, svloc, uloc, &ntot, vloc, &ntot, 
				workloc, &lworkloc, rworkloc, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
		throw " Error inside Lapack routine ZGESVD";
	};

	CHyst hyst;

	hyst.UpdateHystogram (ntot, svloc);

	cout << " Singular values: " << hyst << endl;

// Solve

	dcmplx caux;
	double daux;

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (kjj=0;kjj<ntot;kjj++) {
			caux = czero;
			for (kii=0;kii<ntot;kii++) {
				caux += conj(uloc[kjj*ntot+kii]) * _rhs[irhs*_ldrhs+kii];
			};
			aloc[irhs*ntot+kjj] = caux;
		};
	};

	double alpha_2 = _alpha * _alpha;

	double daux1;

	for (kjj=0;kjj<ntot;kjj++) {
		daux1 = svloc[kjj];
		daux = daux1 / (daux1*daux1 + alpha_2);
		for (irhs=0;irhs<_nrhs;irhs++) {
			caux = aloc[irhs*ntot+kjj] * daux;
			aloc[irhs*ntot+kjj] = caux;
		};
	};

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (kjj=0;kjj<ntot;kjj++) {
			caux = czero;
			for (kii=0;kii<ntot;kii++) {
				caux += conj(vloc[kjj*ntot+kii]) * aloc[irhs*ntot+kii];
			};
			_sol[irhs*_ldsol+kjj] = caux;
		};
	};

// Free work arrays

	delete [] aloc;
	delete [] uloc;
	delete [] vloc;
	delete [] svloc;
	delete [] workloc;
	delete [] rworkloc;

};

// Author: Kharchenko S.A.
// Description: Compute SVD of the current block
// SvdBlk()
//========================================================================================
void SvdBlk (int _m, int _n, // Compute SVD of the current block
				double *_a, double *_sv, double *_u, double *_v, 
				double *_vwork, double *_tau, 
				int _lwork, double *_work) {

// Compute QR decomposition in-place

	QrdBlk (_m, _n, _a, _tau);

// Move R part into vwork

	int i;

	for (i=0;i<_n*_n;i++) _vwork[i] = 0.0e0;

	for (i=0;i<_n;i++) {
		for (int j=0;j<=i;j++) {
			_vwork[i*_n+j] = _a[i*_m+j];
		};
	};

// Compute SVD

	int info = 0;

	dgesvd_ ("O", "A", &_n, &_n, 
				_vwork, &_n, 
				_sv, _vwork, &_n, _v, &_n, 
				_work, &_lwork, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine DGESVD, info = " << info << endl;
		throw " Error inside Lapack routine DGESVD";
	};

// Multiply the result by U

	MvmQBlk (_m, _n, _n, 
				_a, _tau, 
				_vwork, _u);

// Copy data to a

	for (i=0;i<_m*_n;i++) _a[i] = _u[i];

};

// Author: Kharchenko S.A.
// Description: Compute SVD of the current block
// SvdBlk()
//========================================================================================
void SvdBlk (int _m, int _n, // Compute SVD of the current block
				dcmplx *_a, double *_sv, dcmplx *_u, dcmplx *_v, 
				dcmplx *_vwork, dcmplx *_tau, 
				int _lwork, dcmplx *_work, double *_rwork) {

	dcmplx czero (0.0e0,0.0e0);

// Compute QR decomposition in-place

	QrdBlk (_m, _n, _a, _m, _tau);

// Move R part into vwork

	int i;

	for (i=0;i<_n*_n;i++) _vwork[i] = czero;

	for (i=0;i<_n;i++) {
		for (int j=0;j<=i;j++) {
			_vwork[i*_n+j] = _a[i*_m+j];
		};
	};

// Compute SVD

	int info = 0;

	zgesvd_ ("O", "A", &_n, &_n, 
				_vwork, &_n, 
				_sv, _vwork, &_n, _v, &_n, 
				_work, &_lwork, _rwork, &info);
	if (info != 0) {
		cout << " Error inside Lapack routine ZGESVD, info = " << info << endl;
		throw " Error inside Lapack routine ZGESVD";
	};

// Multiply the result by U

	MvmQBlk (_m, _n, _n, 
				_a, _m, _tau, 
				_vwork, _n, _u, _m);

// Copy data to a

	for (i=0;i<_m*_n;i++) _a[i] = _u[i];

};

// Author: Kharchenko S.A.
// Description: Compute QR decomposition of the block
// QrdBlk()
//========================================================================================
void QrdBlk (int _m, int _n, double *_a, double *_tau) { // Compute QR decomposition of the block

	int i, j, k, kk;
	double scprod;

	double *pa, *ph;

// Main cycle over columns

	for (i=0;i<_n;i++) {

// Apply previous columns to the current column

		for (j=0;j<i;j++) {

			scprod = _a[i*_m+j];

			pa = _a+i*_m+j+1;
			ph = _a+j*_m+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
				scprod += *pa++ * *ph++;
			};

			scprod *= _tau[j];

			_a[i*_m+j] -= scprod;

			pa = _a+i*_m+j+1;
			ph = _a+j*_m+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
				*pa++ -= scprod * *ph++;
			};

		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		dlarfg_ (&k, _a+i*_m+i, _a+i*_m+j, &kk, _tau+i);
//		DLARFG (&k, _a+i*_m+i, _a+i*_m+j, &kk, _tau+i);

	};

};

// Author: Kharchenko S.A.
// Description: Compute QR decomposition of the complex block
// QrdBlk()
//========================================================================================
void QrdBlk (int _m, int _n, double *_a, int _lda, double *_tau) { // Compute QR decomposition of the block

	int i, j, k, kk;
	double scprod;

	double *pa, *ph;

// Main cycle over columns

	for (i=0;i<_n;i++) {

// Apply previous columns to the current column

		for (j=0;j<i;j++) {

			scprod = _a[i*_lda+j];

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
				scprod += (*pa++ * *ph++);
//				scprod.PlEq (*pa++ * conj(*ph++));
			};

			scprod *= _tau[j];
//			scprod.MlEq (conj(_tau[j]));

			_a[i*_lda+j] -= scprod;

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
				*pa++ -= scprod * *ph++;
//				(*pa++).MnEq (scprod * *ph++);
			};

		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		dlarfg_ (&k, _a+i*_lda+i, _a+i*_lda+j, &kk, _tau+i);

	};

};

// Author: Kharchenko S.A.
// Description: Update QR decomposition of the complex block
// UpdateQrd()
//========================================================================================
void UpdateQrd (int _m, int _ibeg, int _iend, double *_a, int _lda, double *_tau) { // Update QR decomposition of the block

	int i, j, k, kk;
	double scprod;

	double *pa, *ph;

// Main cycle over columns

	for (i=_ibeg;i<_iend;i++) {

// Apply previous columns to the current column

		for (j=0;j<i;j++) {

			scprod = _a[i*_lda+j];

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
				scprod += *pa++ * *ph++;
//				scprod.PlEq (*pa++ * conj(*ph++));
			};

			scprod *= _tau[j];
//			scprod.MlEq (conj(_tau[j]));

			_a[i*_lda+j] -= scprod;

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
				*pa++ -= scprod * *ph++;
//				(*pa++).MnEq (scprod * *ph++);
			};

		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		dlarfg_ (&k, _a+i*_lda+i, _a+i*_lda+j, &kk, _tau+i);

	};

};

// Author: Kharchenko S.A.
// Description: Update QR decomposition of the complex block in the out-of-core mode
// UpdateQrd()
//========================================================================================
void UpdateQrd (int _m, int _iblk, int _blksiz, // Update QR decomposition of the block in the out-of-core mode
				FILE **_afiles, int *_bl2file, int *_bsblk, int _lda, double *_tau,
				double *_abuff) {

	int i, j, k, kk;
	double scprod;

	double *pa, *ph;

	double *pacurr;

	int ncolcomp = _blksiz;
	int ncolread = _blksiz;

	pacurr = _abuff+_lda*ncolread;

	int iblkrd;

// Main cycle over columns

	int ifile = _bl2file[_iblk];
	int ibs = _bsblk[_iblk];

	FGet (_afiles[ifile],_lda*ncolcomp,pacurr,ibs);

// Apply previous columns to the current column

	for (iblkrd=0;iblkrd<_iblk;iblkrd++) {

		int jfile = _bl2file[iblkrd];
		int jbs = _bsblk[iblkrd];

		FGet (_afiles[jfile],_lda*ncolread,_abuff,jbs);

		for (i=_iblk*_blksiz;i<(_iblk+1)*_blksiz;i++) {

			for (j=iblkrd*_blksiz;j<(iblkrd+1)*_blksiz;j++) {

				scprod = pacurr[(i-_iblk*_blksiz)*_lda+j];

				pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
				ph = _abuff+(j-iblkrd*_blksiz)*_lda+j+1;

				for (k=j+1;k<_m;k++) {
//					scprod += _a[i*_m+k] * _a[j*_m+k];
					scprod += *pa++ * *ph++;
//					scprod.PlEq (*pa++ * conj(*ph++));
				};

				scprod *= _tau[j];
//				scprod.MlEq (conj(_tau[j]));

				pacurr[(i-_iblk*_blksiz)*_lda+j] -= scprod;

				pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
				ph = _abuff+(j-iblkrd*_blksiz)*_lda+j+1;

				for (k=j+1;k<_m;k++) {
//					_a[i*_m+k] -= scprod*_a[j*_m+k];
//					*pa++ -= scprod * conj(*ph++);
					*pa++ -= scprod * *ph++;
//					(*pa++).MnEq (scprod * *ph++);
				};
			};

		};

	};

	for (i=_iblk*_blksiz;i<(_iblk+1)*_blksiz;i++) {

// Apply previous columns stored in-core to the current column

		for (j=_iblk*_blksiz;j<i;j++) {

			scprod = pacurr[(i-_iblk*_blksiz)*_lda+j];

			pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
			ph = pacurr+(j-_iblk*_blksiz)*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
				scprod += *pa++ * *ph++;
//				scprod.PlEq (*pa++ * conj(*ph++));
			};

			scprod *= _tau[j];
//			scprod.MlEq (conj(_tau[j]));

			pacurr[(i-_iblk*_blksiz)*_lda+j] -= scprod;

			pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
			ph = pacurr+(j-_iblk*_blksiz)*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
				*pa++ -= scprod * *ph++;
//				(*pa++).MnEq (scprod * *ph++);
			};
		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		dlarfg_ (&k, pacurr+(i-_iblk*_blksiz)*_lda+i, pacurr+(i-_iblk*_blksiz)*_lda+j, &kk, _tau+i);

	};

	ifile = _bl2file[_iblk];
	ibs = _bsblk[_iblk];

	FPut (_afiles[ifile],_lda*ncolcomp,pacurr,ibs);

};

// Author: Kharchenko S.A.
// Description: Multiply Q factor by the current block
// MvmQBlk()
//========================================================================================
void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block
				double *_q, double *_tau, 
				double *_x, double *_qx) { 

	int irhs, j, k;
	double scprod;

	double *pq, *pqx;

// Init qx

	pqx = _qx;

	for (j=0;j<_m*_nrhs;j++) *pqx++ = 0.0e0;

// Main cycle over columns

	for (irhs=0;irhs<_n;irhs++) {

// Init new column

		for (j=0;j<_n;j++) {
			_qx[irhs*_m+j] = _x[irhs*_n+j];
		};

// Apply Q to the current column

		for (j=_n-1;j>=0;j--) {

			scprod = _qx[irhs*_m+j];

			pq = _q+j*_m+j+1;
			pqx = _qx+irhs*_m+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _q[j*_m+k] * _qx[irhs*_m+k];
				scprod += *pq++ * *pqx++;
			};

			scprod = -scprod * _tau[j];

			_qx[irhs*_m+j] += scprod;

			pq = _q+j*_m+j+1;
			pqx = _qx+irhs*_m+j+1;

			for (k=j+1;k<_m;k++) {
//				_qx[irhs*_m+k] += scprod*_q[j*_m+k];
				*pqx++ += scprod * *pq++;
			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Multiply Q factor by the current block
// MvmQBlk()
//========================================================================================
void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block
				double *_q, int _ldq, double *_tau, 
				double *_x, int _ldx, double *_qx, int _ldqx) { 

	double dzero = 0.0e0;

	int irhs, j, k;
	double scprod;

	double *pq, *pqx;

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<_m;j++) *pqx++ = dzero;
	};

// Main cycle over columns

	for (irhs=0;irhs<_nrhs;irhs++) {

// Init new column

		for (j=0;j<_n;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};

// Apply Q to the current column

		for (j=_n-1;j>=0;j--) {

			scprod = _qx[irhs*_ldqx+j];

			pq = _q+j*_ldq+j+1;
			pqx = _qx+irhs*_ldqx+j+1;

			for (k=j+1;k<_m;k++) {
				scprod += *pq++ * *pqx++;
//				scprod.PlEq (*pqx++ * conj(*pq++));
			};

			scprod = -scprod * _tau[j];

			_qx[irhs*_ldqx+j] += scprod;

			pq = _q+j*_ldq+j+1;
			pqx = _qx+irhs*_ldqx+j+1;

			for (k=j+1;k<_m;k++) {
				*pqx++ += scprod * *pq++;
//				(*pqx++).PlEq (scprod * *pq++);
			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Multiply Q factor by the current block in the out-of-core mode
// MvmQBlk()
//========================================================================================
void MvmQBlk (int _m, int _nblk, int _blksiz, int _nrhs, // Multiply Q factor by the current block in the out-of-core mode
				FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, double *_tau, 
				double *_x, int _ldx, double *_qx, int _ldqx,
				double *_qbuff) { 

	double dzero = 0.0e0;

	int irhs, j, k;
	double scprod;

	double *pq, *pqx;

	if (_blksiz != _nrhs) {
		cout << " Block size and the number of right hand sides do not coincide " << endl;
		throw " Block size and the number of right hand sides do not coincide ";
	};

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<_m;j++) *pqx++ = dzero;
		for (j=0;j<_nblk*_blksiz;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};
	};

// Main cycle over columns

	int iblkrd;

	for (iblkrd=_nblk-1;iblkrd>=0;iblkrd--) {

		int ifile = _bl2file[iblkrd];
		int ibs = _bsblk[iblkrd];

		FGet (_qfiles[ifile],_ldq*_blksiz,_qbuff,ibs);

// Apply current part of Q to the set of columns

		for (j=(iblkrd+1)*_blksiz-1;j>=iblkrd*_blksiz;j--) {

			for (irhs=0;irhs<_nrhs;irhs++) {

				scprod = _qx[irhs*_ldqx+j];

				pq = _qbuff+(j-iblkrd*_blksiz)*_ldq+j+1;
				pqx = _qx+irhs*_ldqx+j+1;

				for (k=j+1;k<_m;k++) {
					scprod += *pq++ * *pqx++;
//					scprod.PlEq (*pqx++ * conj(*pq++));
				};

				scprod = -scprod * _tau[j];

				_qx[irhs*_ldqx+j] += scprod;

				pq = _qbuff+(j-iblkrd*_blksiz)*_ldq+j+1;
				pqx = _qx+irhs*_ldqx+j+1;

				for (k=j+1;k<_m;k++) {
					*pqx++ += scprod * *pq++;
//					(*pqx++).PlEq (scprod * *pq++);
				};

			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Get R part of the QR decomposition
// GetRPartQrd()
//========================================================================================
void GetRPartQrd (int _m, int _n, // Get R part of the QR decomposition
					double *_q, double *_r) {

	int kii, kjj;

	for (kii=0;kii<_n*_n;kii++) _r[kii] = 0.0e0;

	for (kii=0;kii<_n;kii++) {
		for (kjj=0;kjj<=kii;kjj++) {
			_r[kii*_n+kjj] = _q[kii*_m+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Get R part of the QR decomposition
// GetRPartQrd()
//========================================================================================
void GetRPartQrd (int _m, int _ibeg, int _iend, // Get R part of the QR decomposition
					double *_q, int _ldq, double *_r, int _ldr) {

	double dzero = 0.0e0;

	int kii, kjj;

	kii = _m;

	for (kii=0;kii<_iend-_ibeg;kii++) {
		for (kjj=0;kjj<_iend;kjj++) {
			_r[kii*_ldr+kjj] = dzero;
		};
	};

	for (kii=_ibeg;kii<_iend;kii++) {
		for (kjj=0;kjj<=kii;kjj++) {
			_r[(kii-_ibeg)*_ldr+kjj] = _q[kii*_ldq+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Get R part of the QR decomposition in the out-of-core mode
// GetRPartQrd()
//========================================================================================
void GetRPartQrd (int _m, int _iblkbeg, int _iblkend, int _blksiz, // Get R part of the QR decomposition in the out-of-core mode
					FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, double *_r, int _ldr) {

	double dzero = 0.0e0;

	int iblk, kii, kjj;

	kii = _m;

	for (kii=0;kii<(_iblkend-_iblkbeg)*_blksiz;kii++) {
		for (kjj=0;kjj<_iblkend*_blksiz;kjj++) {
			_r[kii*_ldr+kjj] = dzero;
		};
	};

	for (iblk=_iblkbeg;iblk<_iblkend;iblk++) {

		int ifile = _bl2file[iblk];
		int ibs = _bsblk[iblk];

		for (kii=iblk*_blksiz;kii<(iblk+1)*_blksiz;kii++) {
			FGet (_qfiles[ifile],kii+1,_r+(kii-_iblkbeg*_blksiz)*_ldr,ibs+(kii-iblk*_blksiz)*_ldq);
		};

	};

};

// Author: Kharchenko S.A.
// Description: Compute QR decomposition of the complex block
// QrdBlk()
//========================================================================================
void QrdBlk (int _m, int _n, dcmplx *_a, int _lda, dcmplx *_tau) { // Compute QR decomposition of the block

	int i, j, k, kk;
	dcmplx scprod;

	dcmplx *pa, *ph;

// Main cycle over columns

	for (i=0;i<_n;i++) {

// Apply previous columns to the current column

		for (j=0;j<i;j++) {

			scprod = _a[i*_lda+j];

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
//				scprod += (*pa++ * conj(*ph++));
				scprod.PlEq (*pa++ * conj(*ph++));
			};

//			scprod *= conj(_tau[j]);
			scprod.MlEq (conj(_tau[j]));

			_a[i*_lda+j] -= scprod;

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * *ph++;
				(*pa++).MnEq (scprod * *ph++);
			};

		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		zlarfg_ (&k, _a+i*_lda+i, _a+i*_lda+j, &kk, _tau+i);

	};

};

// Author: Kharchenko S.A.
// Description: Update QR decomposition of the complex block
// UpdateQrd()
//========================================================================================
void UpdateQrd (int _m, int _ibeg, int _iend, dcmplx *_a, int _lda, dcmplx *_tau) { // Update QR decomposition of the block

	int i, j, k, kk;
	dcmplx scprod;

	dcmplx *pa, *ph;

// Main cycle over columns

	for (i=_ibeg;i<_iend;i++) {

// Apply previous columns to the current column

		for (j=0;j<i;j++) {

			scprod = _a[i*_lda+j];

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
//				scprod += *pa++ * conj(*ph++);
				scprod.PlEq (*pa++ * conj(*ph++));
			};

//			scprod *= conj(_tau[j]);
			scprod.MlEq (conj(_tau[j]));

			_a[i*_lda+j] -= scprod;

			pa = _a+i*_lda+j+1;
			ph = _a+j*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
//				*pa++ -= scprod * *ph++;
				(*pa++).MnEq (scprod * *ph++);
			};

		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		zlarfg_ (&k, _a+i*_lda+i, _a+i*_lda+j, &kk, _tau+i);

	};

};

// Author: Kharchenko S.A.
// Description: Update QR decomposition of the block in the out-of-core mode
// UpdateQrd()
//========================================================================================
void UpdateQrd (int _m, int _ibeg, int _iend, // Update QR decomposition of the block in the out-of-core mode
				FILE *_afile, int _lda, dcmplx *_tau,
				dcmplx *_abuff, int _ncolbuff) {

	int i, j, k, kk;
	dcmplx scprod;

	dcmplx *pa, *ph;

	dcmplx *pacurr;

	if (_ncolbuff <= _iend-_ibeg) {
		cout << " Insufficient number of work vectors in the exchange buffer " << endl;
		throw " Insufficient number of work vectors in the exchange buffer ";
	};

	int ncolcomp = _iend-_ibeg;
	int ncolread = _ncolbuff-ncolcomp;

	pacurr = _abuff+_lda*ncolread;

	int ibegrd, iendrd, icolrd;

// Main cycle over columns

	FGet (_afile,_lda*ncolcomp,pacurr,_ibeg*_lda);

// Apply previous columns to the current column

	iendrd = -1;

	while (iendrd < _ibeg-1) {

		ibegrd = iendrd+1;
		iendrd = ibegrd+ncolread-1;

		if (iendrd > _ibeg-1) iendrd = _ibeg-1;

		icolrd = iendrd-ibegrd+1;

		FGet (_afile,_lda*icolrd,_abuff,ibegrd*_lda);

		for (i=_ibeg;i<_iend;i++) {

			for (j=ibegrd;j<=iendrd;j++) {

				scprod = pacurr[(i-_ibeg)*_lda+j];

				pa = pacurr+(i-_ibeg)*_lda+j+1;
				ph = _abuff+(j-ibegrd)*_lda+j+1;

				for (k=j+1;k<_m;k++) {
//					scprod += _a[i*_m+k] * _a[j*_m+k];
//					scprod += (*pa++ * conj(*ph++));
					scprod.PlEq (*pa++ * conj(*ph++));
				};

//				scprod *= conj(_tau[j]);
				scprod.MlEq (conj(_tau[j]));

				pacurr[(i-_ibeg)*_lda+j] -= scprod;

				pa = pacurr+(i-_ibeg)*_lda+j+1;
				ph = _abuff+(j-ibegrd)*_lda+j+1;

				for (k=j+1;k<_m;k++) {
//					_a[i*_m+k] -= scprod*_a[j*_m+k];
//					*pa++ -= scprod * conj(*ph++);
//					*pa++ -= scprod * *ph++;
					(*pa++).MnEq (scprod * *ph++);
				};
			};

		};

	};

	for (i=_ibeg;i<_iend;i++) {

// Apply previous columns stored in-core to the current column

		for (j=_ibeg;j<i;j++) {

			scprod = pacurr[(i-_ibeg)*_lda+j];

			pa = pacurr+(i-_ibeg)*_lda+j+1;
			ph = pacurr+(j-_ibeg)*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
//				scprod += *pa++ * conj(*ph++);
				scprod.PlEq (*pa++ * conj(*ph++));
			};

//			scprod *= conj(_tau[j]);
			scprod.MlEq (conj(_tau[j]));

			pacurr[(i-_ibeg)*_lda+j] -= scprod;

			pa = pacurr+(i-_ibeg)*_lda+j+1;
			ph = pacurr+(j-_ibeg)*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
//				*pa++ -= scprod * *ph++;
				(*pa++).MnEq (scprod * *ph++);
			};
		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		zlarfg_ (&k, pacurr+(i-_ibeg)*_lda+i, pacurr+(i-_ibeg)*_lda+j, &kk, _tau+i);

	};

	FPut (_afile,_lda*ncolcomp,pacurr,_ibeg*_lda);

};

// Author: Kharchenko S.A.
// Description: Update QR decomposition of the block in the out-of-core mode
// UpdateQrd()
//========================================================================================
void UpdateQrd (int _m, int _iblk, int _blksiz, // Update QR decomposition of the block in the out-of-core mode
				FILE **_afiles, int *_bl2file, int *_bsblk, int _lda, dcmplx *_tau,
				dcmplx *_abuff) {

	int i, j, k, kk;
	dcmplx scprod;

	dcmplx *pa, *ph;

	dcmplx *pacurr;

	int ncolcomp = _blksiz;
	int ncolread = _blksiz;

	pacurr = _abuff+_lda*ncolread;

	int iblkrd;

// Main cycle over columns

	int ifile = _bl2file[_iblk];
	int ibs = _bsblk[_iblk];

	FGet (_afiles[ifile],_lda*ncolcomp,pacurr,ibs);

// Apply previous columns to the current column

	for (iblkrd=0;iblkrd<_iblk;iblkrd++) {

		int jfile = _bl2file[iblkrd];
		int jbs = _bsblk[iblkrd];

		FGet (_afiles[jfile],_lda*ncolread,_abuff,jbs);

		for (i=_iblk*_blksiz;i<(_iblk+1)*_blksiz;i++) {

			for (j=iblkrd*_blksiz;j<(iblkrd+1)*_blksiz;j++) {

				scprod = pacurr[(i-_iblk*_blksiz)*_lda+j];

				pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
				ph = _abuff+(j-iblkrd*_blksiz)*_lda+j+1;

				for (k=j+1;k<_m;k++) {
//					scprod += _a[i*_m+k] * _a[j*_m+k];
//					scprod += *pa++ * conj(*ph++);
					scprod.PlEq (*pa++ * conj(*ph++));
				};

//				scprod *= conj(_tau[j]);
				scprod.MlEq (conj(_tau[j]));

				pacurr[(i-_iblk*_blksiz)*_lda+j] -= scprod;

				pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
				ph = _abuff+(j-iblkrd*_blksiz)*_lda+j+1;

				for (k=j+1;k<_m;k++) {
//					_a[i*_m+k] -= scprod*_a[j*_m+k];
//					*pa++ -= scprod * conj(*ph++);
//					*pa++ -= scprod * *ph++;
					(*pa++).MnEq (scprod * *ph++);
				};
			};

		};

	};

	for (i=_iblk*_blksiz;i<(_iblk+1)*_blksiz;i++) {

// Apply previous columns stored in-core to the current column

		for (j=_iblk*_blksiz;j<i;j++) {

			scprod = pacurr[(i-_iblk*_blksiz)*_lda+j];

			pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
			ph = pacurr+(j-_iblk*_blksiz)*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += _a[i*_m+k] * _a[j*_m+k];
//				scprod += *pa++ * conj(*ph++);
				scprod.PlEq (*pa++ * conj(*ph++));
			};

//			scprod *= conj(_tau[j]);
			scprod.MlEq (conj(_tau[j]));

			pacurr[(i-_iblk*_blksiz)*_lda+j] -= scprod;

			pa = pacurr+(i-_iblk*_blksiz)*_lda+j+1;
			ph = pacurr+(j-_iblk*_blksiz)*_lda+j+1;

			for (k=j+1;k<_m;k++) {
//				_a[i*_m+k] -= scprod*_a[j*_m+k];
//				*pa++ -= scprod * conj(*ph++);
//				*pa++ -= scprod * *ph++;
				(*pa++).MnEq (scprod * *ph++);
			};
		};

// Compute new transformation

		j = _m-1;
		if (i+1 < _m-1) j = i+1;

		k = _m-i;
		kk = 1;

		zlarfg_ (&k, pacurr+(i-_iblk*_blksiz)*_lda+i, pacurr+(i-_iblk*_blksiz)*_lda+j, &kk, _tau+i);

	};

	ifile = _bl2file[_iblk];
	ibs = _bsblk[_iblk];

	FPut (_afiles[ifile],_lda*ncolcomp,pacurr,ibs);

};

// Author: Kharchenko S.A.
// Description: Multiply Q factor by the current block
// MvmQBlk()
//========================================================================================
void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block
				dcmplx *_q, int _ldq, dcmplx *_tau, 
				dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx) { 

	dcmplx czero (0.0e0,0.0e0);

	int irhs, j, k;
	dcmplx scprod;

	dcmplx *pq, *pqx;

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<_m;j++) *pqx++ = czero;
	};

// Main cycle over columns

	for (irhs=0;irhs<_nrhs;irhs++) {

// Init new column

		for (j=0;j<_n;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};

// Apply Q to the current column

		for (j=_n-1;j>=0;j--) {

			scprod = _qx[irhs*_ldqx+j];

			pq = _q+j*_ldq+j+1;
			pqx = _qx+irhs*_ldqx+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += conj(*pq++) * *pqx++;
				scprod.PlEq (*pqx++ * conj(*pq++));
			};

			scprod = -scprod * _tau[j];

			_qx[irhs*_ldqx+j] += scprod;

			pq = _q+j*_ldq+j+1;
			pqx = _qx+irhs*_ldqx+j+1;

			for (k=j+1;k<_m;k++) {
//				*pqx++ += scprod * *pq++;
				(*pqx++).PlEq (scprod * *pq++);
			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Multiply QH factor by the current block
// MvmQHBlk()
//========================================================================================
void MvmQHBlk (int _m, int _n, int _nrhs, // Multiply QH factor by the current block
					dcmplx *_q, int _ldq, dcmplx *_tau, 
					dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx) { 

	dcmplx czero (0.0e0,0.0e0);

	int irhs, j, k;
	dcmplx scprod, caux;

	dcmplx *pq, *pqx;

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<_m;j++) *pqx++ = czero;
	};

// Main cycle over columns

	for (irhs=0;irhs<_nrhs;irhs++) {

// Init new column

		for (j=0;j<_m;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};

// Apply Q to the current column

		for (j=0;j<_n;j++) {

			scprod = _qx[irhs*_ldqx+j];

			pq = _q+j*_ldq+j+1;
			pqx = _qx+irhs*_ldqx+j+1;

			for (k=j+1;k<_m;k++) {
//				scprod += conj(*pq++) * *pqx++;
				scprod.PlEq (*pqx++ * conj(*pq++));
			};

			caux = scprod * conj(_tau[j]);
			scprod = -caux;

			_qx[irhs*_ldqx+j] += scprod;

			pq = _q+j*_ldq+j+1;
			pqx = _qx+irhs*_ldqx+j+1;

			for (k=j+1;k<_m;k++) {
//				*pqx++ += scprod * *pq++;
				(*pqx++).PlEq (scprod * *pq++);
			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Multiply QH factor by the current block in the out-of-core mode
// MvmQBlk()
//========================================================================================
void MvmQBlk (int _m, int _n, int _nrhs, // Multiply Q factor by the current block in the out-of-core mode
				FILE *_qfile, int _ldq, dcmplx *_tau, 
				dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx,
				dcmplx *_qbuff, int _ncolbuff) { 

	dcmplx czero (0.0e0,0.0e0);

	int irhs, j, k;
	dcmplx scprod;

	dcmplx *pq, *pqx;

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<_m;j++) *pqx++ = czero;
		for (j=0;j<_n;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};
	};

// Main cycle over columns

	int ibegrd, iendrd, icolrd;

	ibegrd = _n;

	while (ibegrd > 0) {

		iendrd = ibegrd - 1;
		ibegrd = iendrd-_ncolbuff+1;
		if (ibegrd <0) ibegrd = 0;

		icolrd = iendrd-ibegrd+1;

		FGet (_qfile,icolrd*_ldq,_qbuff,ibegrd*_ldq);

// Apply current part of Q to the set of columns

		for (j=iendrd;j>=ibegrd;j--) {

			for (irhs=0;irhs<_nrhs;irhs++) {

				scprod = _qx[irhs*_ldqx+j];

				pq = _qbuff+(j-ibegrd)*_ldq+j+1;
				pqx = _qx+irhs*_ldqx+j+1;

				for (k=j+1;k<_m;k++) {
//					scprod += conj(*pq++) * *pqx++;
					scprod.PlEq (*pqx++ * conj(*pq++));
				};

				scprod = -scprod * _tau[j];

				_qx[irhs*_ldqx+j] += scprod;

				pq = _qbuff+(j-ibegrd)*_ldq+j+1;
				pqx = _qx+irhs*_ldqx+j+1;

				for (k=j+1;k<_m;k++) {
//					*pqx++ += scprod * *pq++;
					(*pqx++).PlEq (scprod * *pq++);
				};

			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Multiply QH factor by the current block in the out-of-core mode
// MvmQBlk()
//========================================================================================
void MvmQBlk (int _m, int _nblk, int _blksiz, int _nrhs, // Multiply Q factor by the current block in the out-of-core mode
				FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, dcmplx *_tau, 
				dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx,
				dcmplx *_qbuff) { 

	dcmplx czero (0.0e0,0.0e0);

	int irhs, j, k;
	dcmplx scprod;

	dcmplx *pq, *pqx;

	if (_blksiz != _nrhs) {
		cout << " Block size and the number of right hand sides do not coincide " << endl;
		throw " Block size and the number of right hand sides do not coincide ";
	};

// Init qx

	for (irhs=0;irhs<_nrhs;irhs++) {
		pqx = _qx+irhs*_ldqx;
		for (j=0;j<_m;j++) *pqx++ = czero;
		for (j=0;j<_nblk*_blksiz;j++) {
			_qx[irhs*_ldqx+j] = _x[irhs*_ldx+j];
		};
	};

// Main cycle over columns

	int iblkrd;

	for (iblkrd=_nblk-1;iblkrd>=0;iblkrd--) {

		int ifile = _bl2file[iblkrd];
		int ibs = _bsblk[iblkrd];

		FGet (_qfiles[ifile],_ldq*_blksiz,_qbuff,ibs);

// Apply current part of Q to the set of columns

		for (j=(iblkrd+1)*_blksiz-1;j>=iblkrd*_blksiz;j--) {

			for (irhs=0;irhs<_nrhs;irhs++) {

				scprod = _qx[irhs*_ldqx+j];

				pq = _qbuff+(j-iblkrd*_blksiz)*_ldq+j+1;
				pqx = _qx+irhs*_ldqx+j+1;

				for (k=j+1;k<_m;k++) {
//					scprod += conj(*pq++) * *pqx++;
					scprod.PlEq (*pqx++ * conj(*pq++));
				};

				scprod = -scprod * _tau[j];

				_qx[irhs*_ldqx+j] += scprod;

				pq = _qbuff+(j-iblkrd*_blksiz)*_ldq+j+1;
				pqx = _qx+irhs*_ldqx+j+1;

				for (k=j+1;k<_m;k++) {
//					*pqx++ += scprod * *pq++;
					(*pqx++).PlEq (scprod * *pq++);
				};

			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Get R part of the QR decomposition
// GetRPartQrd()
//========================================================================================
void GetRPartQrd (int _m, int _ibeg, int _iend, // Get R part of the QR decomposition
					dcmplx *_q, int _ldq, dcmplx *_r, int _ldr) {

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj;

	kii = _m;

	for (kii=0;kii<_iend-_ibeg;kii++) {
		for (kjj=0;kjj<_iend;kjj++) {
			_r[kii*_ldr+kjj] = czero;
		};
	};

	for (kii=_ibeg;kii<_iend;kii++) {
		for (kjj=0;kjj<=kii;kjj++) {
			_r[(kii-_ibeg)*_ldr+kjj] = _q[kii*_ldq+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Get R part of the QR decomposition in the out-of-core mode
// GetRPartQrd()
//========================================================================================
void GetRPartQrd (int _m, int _ibeg, int _iend, // Get R part of the QR decomposition in the out-of-core mode
					FILE *_qfile, int _ldq, dcmplx *_r, int _ldr) {

	dcmplx czero (0.0e0,0.0e0);

	int kii, kjj;

	kii = _m;

	for (kii=0;kii<_iend-_ibeg;kii++) {
		for (kjj=0;kjj<_iend;kjj++) {
			_r[kii*_ldr+kjj] = czero;
		};
	};

	for (kii=_ibeg;kii<_iend;kii++) {

		FGet (_qfile,kii+1,_r+(kii-_ibeg)*_ldr,kii*_ldq);

	};

};

// Author: Kharchenko S.A.
// Description: Get R part of the QR decomposition in the out-of-core mode
// GetRPartQrd()
//========================================================================================
void GetRPartQrd (int _m, int _iblkbeg, int _iblkend, int _blksiz, // Get R part of the QR decomposition in the out-of-core mode
					FILE **_qfiles, int *_bl2file, int *_bsblk, int _ldq, dcmplx *_r, int _ldr) {

	dcmplx czero (0.0e0,0.0e0);

	int iblk, kii, kjj;

	kii = _m;

	for (kii=0;kii<(_iblkend-_iblkbeg)*_blksiz;kii++) {
		for (kjj=0;kjj<_iblkend*_blksiz;kjj++) {
			_r[kii*_ldr+kjj] = czero;
		};
	};

	for (iblk=_iblkbeg;iblk<_iblkend;iblk++) {

		int ifile = _bl2file[iblk];
		int ibs = _bsblk[iblk];

		for (kii=iblk*_blksiz;kii<(iblk+1)*_blksiz;kii++) {
			FGet (_qfiles[ifile],kii+1,_r+(kii-_iblkbeg*_blksiz)*_ldr,ibs+(kii-iblk*_blksiz)*_ldq);
		};

	};

};

// Author: Kharchenko S.A.
// Description: transposed sparsity 
// TransposeMatrix()
//========================================================================================
void TransposeMatrix (int _n, int *_ia, int *_ja, int *_iat, int *_jat) { // Compute transposed sparsity

//	const char *funcname = "TransposeMatrix";

	int i, j, jj;

	for (i=0;i<=_n;i++) _iat[i] = 0;

	for (i=0;i<_n;i++) {
		for (j=_ia[i];j<_ia[i+1];j++) {
			jj = _ja[j];
			_iat[jj+1]++;
		};
	};

	for (i=0;i<_n;i++) _iat[i+1] += _iat[i];

	int *pat;

	pat = new int [_n];
	if (!pat) throw " Memory allocation error";

	for (i=0;i<_n;i++) pat[i] = _iat[i];
	
	int k;
	for (i=0;i<_n;i++) {
		for (j=_ia[i];j<_ia[i+1];j++) {
			jj = _ja[j];
			pat[jj]++;
			k = pat[jj];
			_jat[k-1] = i;
		};
	};

	delete [] pat;

};

// Author: Kharchenko S.A.
// Description: Create the histogram of values
// csvhyst()
//========================================================================================
void csvhyst (ostream &_fout, const char *_name, int _n, double *_sv) { // Create the histogram of values

	const int nhyst = 16;

	int hyst[nhyst];
	double thre[nhyst];

	if (_n < 1) return;

	double smin, smax;

	smin = _sv[0];
	smax = _sv[0];

	int i;

	for (i=0;i<_n;i++) {
		if (_sv[i] > smax) smax = _sv[i];
		if (_sv[i] < smin) smin = _sv[i];
	};

	thre[nhyst-1] = smax;

	for (i=0;i<nhyst;i++) hyst[i] = 0;
	for (i=nhyst-2;i>=0;i--) thre[i] = thre[i+1] / 1.0e1;

	for (i=0;i<_n;i++) {
		if (_sv[i] <= thre[0]) hyst[0]++;
		for (int j=1;j<nhyst;j++) {
			if (_sv[i] > thre[j-1] && _sv[i] <= thre[j]) hyst[j]++;
		};
	};

	int iii;

	iii = 0;

	for (i=0;i<nhyst;i++) {
		if (hyst[i] != 0) {
			iii = i;
			goto exit;
		};
	};

exit:;

	_fout << "    " << _name << ": N = " << _n << "; SMin = " << smin << "; SMax = " << smax << ";" << endl;
	_fout << " H: ";
	for (i=iii;i<nhyst;i++) _fout << hyst[i] << " ";
	_fout << endl;

};

// Author: Kharchenko S.A.
// Description: Create the histogram of values
// csvhyst()
//========================================================================================
void csvhyst (const char *_name, int _n, double *_sv) { // Create the histogram of values

	const int nhyst = 16;

	int hyst[nhyst];
	double thre[nhyst];

	if (_n < 1) return;

	double smin, smax;

	smin = _sv[0];
	smax = _sv[0];

	int i;

	for (i=0;i<_n;i++) {
		if (_sv[i] > smax) smax = _sv[i];
		if (_sv[i] < smin) smin = _sv[i];
	};

	thre[nhyst-1] = smax;

	for (i=0;i<nhyst;i++) hyst[i] = 0;
	for (i=nhyst-2;i>=0;i--) thre[i] = thre[i+1] / 1.0e1;

	for (i=0;i<_n;i++) {
		if (_sv[i] <= thre[0]) hyst[0]++;
		for (int j=1;j<nhyst;j++) {
			if (_sv[i] > thre[j-1] && _sv[i] <= thre[j]) hyst[j]++;
		};
	};

	int iii;

	iii = 0;

	for (i=0;i<nhyst;i++) {
		if (hyst[i] != 0) {
			iii = i;
			goto exit;
		};
	};

exit:;

	cout << "    " << _name << ": N = " << _n << "; SMin = " << smin << "; SMax = " << smax << ";" << endl;
	cout << " H: ";
	for (i=iii;i<nhyst;i++) cout << hyst[i] << " ";
	cout << endl;

};

// Author: Kharchenko S.A.
// Description: Memory allocation failed
// MemoryFail()
//========================================================================================
void MemoryFail (const char *_funcname) { // Memory allocation failed

	cout << " Memory allocation failed in function " << _funcname << endl;

#ifdef __FV_2004__
   std::string str1 (" Memory allocation failed in function ");
   std::string str2 (_funcname);
   std::string str = str1 + str2;
   ErrSolverMessage (false, FVE_MATRIXSOLVERERROR,  str);
#endif

	throw " Memory allocation failed";

};

// Author: Kharchenko S.A.
// Description: Reorder boolean array
// OrdBool()
//========================================================================================
void OrdBool (char _ordtyp, int _n, int *_order, bool *_barr, bool *_barro) { // Reorder boolean array
	int i;
	if (_ordtyp == 'D') {
		for (i=0;i<_n;i++) _barro[_order[i]] = _barr[i];
	} else {
		for (i=0;i<_n;i++) _barro[i] = _barr[_order[i]];
	};
};

// Author: Kharchenko S.A.
// Description: Read boolean array from disk
// ReadBool()
//========================================================================================
void ReadBool (istream &_fin, int _n, bool *_barr) { // Read boolean array from disk

	int i;
	char buffer[20]="";
 
	_fin.getline(buffer,10);
	cout << buffer << endl;
	_fin.getline(buffer,10);
	cout << buffer << endl;

	for (i=0; i<_n; i++)
	{
		_fin >> _barr[i];
	}

};

// Author: Kharchenko S.A.
// Description: Read boolean array from disk in binary form
// ReadBoolBin()
//========================================================================================
void ReadBoolBin (istream &_fin, int _n, bool *_barr) { // Read boolean array from disk in binary form

	ReadArrAccV (_fin, _n, _barr);

};

// Author: Kharchenko S.A.
// Description: Collect bool data
// CollectBool()
//========================================================================================
void CollectBool (int _myid, int _nproc, int *_blks, bool *_barr, bool *_barrglob) { // Collect bool data

	const char *funcname = "CollectBool";

// Allocate temporary integer array

	int ntotal = _blks[_nproc+1];
	int *iarr;

	iarr = new int [ntotal];
	if (!iarr) MemoryFail (funcname);

	int i;

	for (i=0;i<ntotal;i++) iarr[i] = 0;

// Init array

	if (_myid != 0) {
		for (i=_blks[_myid];i<_blks[_myid+1];i++) {
			int i1 = i-_blks[_myid];
			if (_barr[i1]) {
				iarr[i] = 1;
			} else {
				iarr[i] = 0;
			};
		};
	} else {
		for (i=_blks[0];i<_blks[1];i++) {
			int i1 = i;
			if (_barr[i1]) {
				iarr[i] = 1;
			} else {
				iarr[i] = 0;
			};
		};
		for (i=_blks[_nproc];i<_blks[_nproc+1];i++) {
			int i1 = i-_blks[_nproc]+_blks[1];
			if (_barr[i1]) {
				iarr[i] = 1;
			} else {
				iarr[i] = 0;
			};
		};
	};

// Combine local integer arrays into one

	AddSparseIArray (_myid, _nproc, ntotal, iarr);

// Store the data in one global boolean array

	for (i=0;i<ntotal;i++) {
		if (iarr[i] != 0) {
			_barrglob[i] = true;
		} else {
			_barrglob[i] = false;
		};
	};

// Delete work array

	delete [] iarr;

};

// Author: Kharchenko S.A.
// Description: Compare two integer values
// compint()
//========================================================================================
int compint (const void *arg1, const void *arg2) { // Compare two integer values
	int *iarg1= (int *) arg1;
	int *iarg2= (int *) arg2;
	if (*iarg1 <  *iarg2) {
		return -1;
	} else if (*iarg1 == *iarg2) {
		return 0;
	} else {
		return 1;
	};
};

// Author: Kharchenko S.A.
// Description: Print integer array
// OutArr()
//========================================================================================
void OutArr (ostream &_stream, const char *_name, int _isize, int *_iarr) { // Print integer array
//	cout << " In OutIArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (Addr = " << _iarr << " Size = " << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
		if ((i%5 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(8) << _iarr[i];
		_stream << setw(12) << _iarr[i];
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print float array
// OutArr()
//========================================================================================
void OutArr (ostream &_stream, const char *_name, int _isize, float *_farr) { // Print float array
//	cout << " In OutDArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (Addr = " << _farr << " Size = " << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
		if ((i%6 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(23) << setprecision(16) << _darr[i] << " ";
		_stream << setw(13) << setprecision(8) << _farr[i] << " ";
//		if ((i%5 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(16) << setprecision(9) << _darr[i];
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print double array
// OutArr()
//========================================================================================
void OutArr (ostream &_stream, const char *_name, int _isize, double *_darr) { // Print double array
//	cout << " In OutDArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (Addr = " << _darr << " Size = " << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
		if ((i%3 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(23) << setprecision(16) << _darr[i] << " ";
		_stream << setw(23) << setprecision(16) << _darr[i] << " ";
//		if ((i%5 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(16) << setprecision(9) << _darr[i];
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print double array with small value threshold
// OutArr()
//========================================================================================
void OutArrThresh (ostream &_stream, const char *_name, int _isize, double *_darr, double _thresh) { // Print double array with small value threshold
//	cout << " In OutDArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (" << _isize << ")" << endl;
	double aux1, aux2;
	for (int i=0;i<_isize;i++) {
		if ((i%3 == 0) && (i > 0)) _stream << endl;
		aux1 = _darr[i];
		aux2 = aux1;
		if (aux2 < 0.0e0) aux2 = -aux2;
		if (aux2 < _thresh) aux1 = 0.0e0;
		_stream << setw(23) << setprecision(16) << aux1 << " ";
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print double complex array
// OutArr()
//========================================================================================
void OutArr (ostream &_stream, const char *_name, int _isize, dcmplx *_carr) { // Print double complex array
//	cout << " In OutDArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (" << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
		if ((i%2 == 0) && (i > 0)) _stream << endl;
//		_stream << left << setw(20) << setprecision(14) << _carr[i] << " ";
//		_stream << setw(20) << setprecision(14) << _carr[i] << " ";
		_stream << setw(20) << setprecision(12) << _carr[i] << " ";
//		if ((i%5 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(16) << setprecision(9) << _carr[i];
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print double complex array
// OutArr()
//========================================================================================
void OutArrThresh (ostream &_stream, const char *_name, int _isize, dcmplx *_carr, double _thresh) { // Print double complex array with small value threshold
//	cout << " In OutDArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (" << _isize << ")" << endl;
	double aux1, aux2;
	for (int i=0;i<_isize;i++) {
		if ((i%2 == 0) && (i > 0)) _stream << endl;
		aux1 = _carr[i].x;
		aux2 = _carr[i].y;
		double daux1 = aux1;
		double daux2 = aux2;
		if (daux1 < 0.0e0) daux1 = -daux1;
		if (daux2 < 0.0e0) daux2 = -daux2;
		if (daux1 < _thresh) aux1 = 0.0e0;
		if (daux2 < _thresh) aux2 = 0.0e0;
		_stream << "(" << setw(22) << setprecision(14) << aux1 << ",";
		_stream << setw(22) << setprecision(14) << aux2 << ") ";
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print boolean array
// OutArr()
//========================================================================================
void OutArr (ostream &_stream, const char *_name, int _isize, bool *_barr) { // Print boolean array
//	cout << " In OutBArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (" << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
		if ((i%20 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(4) << _barr[i];
		_stream << setw(4) << _barr[i];
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print char array
// OutArr()
//========================================================================================
void OutArr (ostream &_stream, const char *_name, int _isize, char *_charr) { // Print char array
//	cout << " In OutDArr, name = " << _name << " isize = " << _isize << endl;
	_stream << _name << " (" << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
		if ((i%20 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(23) << setprecision(16) << _darr[i] << " ";
		_stream << _charr[i] << " ";
//		if ((i%5 == 0) && (i > 0)) _stream << endl;
//		_stream << right << setw(16) << setprecision(9) << _darr[i];
	};
	_stream << endl;
};

// Author: Kharchenko S.A.
// Description: Print integer array
// OutArr()
//========================================================================================
void OutArrAcc (ostream &_stream, const char *_name, int _isize, int *_iarr) { // Store integer array
	_stream << _name << " (" << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
//		_stream << right << setw(8) << _iarr[i] << endl;
		_stream << setw(8) << _iarr[i] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Print double array
// OutArr()
//========================================================================================
void OutArrAcc (ostream &_stream, const char *_name, int _isize, double *_darr) { // Store double array
	_stream << _name << " (" << _isize << ")" << endl;
	for (int i=0;i<_isize;i++) {
//		_stream << right << setw(22) << setprecision(15) << _darr[i] << endl;
		_stream << setw(22) << setprecision(15) << _darr[i] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Store integer array
// OutArr()
//========================================================================================
void OutArrAccV (ostream &_stream, const char *_name, int _isize, int *_iarr) { // Store integer array
	_stream.write ((char*)_iarr,sizeof(int)*_isize);
	if (_isize != 0) {
		cout << _name << " Iarr: First elem = " << _iarr[0] << " Last elem = " << _iarr[_isize-1] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Store double array
// OutArr()
//========================================================================================
void OutArrAccV (ostream &_stream, const char *_name, int _isize, double *_darr) { // Store double array
	_stream.write ((char*)_darr,sizeof(double)*_isize);
	if (_isize != 0) {
		cout << _name << " Darr: First elem = " << _darr[0] << " Last elem = " << _darr[_isize-1] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Store boolean array
// OutArr()
//========================================================================================
void OutArrAccV (ostream &_stream, const char *_name, int _isize, bool *_barr) { // Store boolean array
	_stream.write ((char*)_barr,sizeof(bool)*_isize);
	if (_isize != 0) {
		cout << _name << " Barr: First elem = " << _barr[0] << " Last elem = " << _barr[_isize-1] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Write boolean array to the disk in binary form
// WriteBoolBin()
//========================================================================================
void WriteBoolBin (ostream &_stream, const char *_name, int _n, bool *_barr) { // Read boolean array from disk in binary form

	OutArrAccV (_stream, _name, _n, _barr);

};

// Author: Kharchenko S.A.
// Description: Read integer array
// ReadArrAcc()
//========================================================================================
void ReadArrAcc (istream &_stream, int _isize, int *_iarr) { // Read integer array

	char buffer[20]="";

	_stream.getline(buffer,10);

	cout << buffer << endl;

	for (int i=0;i<_isize;i++) {
		_stream >> _iarr[i];
		_stream.getline(buffer,10);
	};

};

// Author: Kharchenko S.A.
// Description: Read double array
// ReadArrAcc()
//========================================================================================
void ReadArrAcc (istream &_stream, int _isize, double *_darr) { // Read double array

	char buffer[20]="";

	_stream.getline(buffer,10);

	cout << buffer << endl;

	for (int i=0;i<_isize;i++) {
		_stream >> _darr[i];
		_stream.getline(buffer,10);
	};

};

// Author: Kharchenko S.A.
// Description: Read integer array
// ReadArrAcc()
//========================================================================================
void ReadArrAccV (istream &_stream, int _isize, int *_iarr) { // Read integer array
	_stream.read ((char*)_iarr,sizeof(int)*_isize);
	if (_isize != 0) {
		cout << " Read Iarr: First elem = " << _iarr[0] << " Last elem = " << _iarr[_isize-1] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Read double array
// ReadArrAcc()
//========================================================================================
void ReadArrAccV (istream &_stream, int _isize, double *_darr) { // Read double array
	_stream.read ((char*)_darr,sizeof(double)*_isize);
	if (_isize != 0) {
		cout << " Read Darr: First elem = " << _darr[0] << " Last elem = " << _darr[_isize-1] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Read boolean array
// ReadArrAcc()
//========================================================================================
void ReadArrAccV (istream &_stream, int _isize, bool *_barr) { // Read boolean array
	_stream.read ((char*)_barr,sizeof(bool)*_isize);
	if (_isize != 0) {
		cout << " Read Barr: First elem = " << _barr[0] << " Last elem = " << _barr[_isize-1] << endl;
	};
};

// Author: Kharchenko S.A.
// Description: Store double array (direct access)
// FPut()
//========================================================================================
void FPut (fstream &_stream, int _size, double *_darr, int _offset) { // Store double array (direct access)

	int i = _offset;
	i=i*i;

//	_stream.seekg(sizeof(double)*_offset,ios::beg);
//	_stream.seekp(sizeof(double)*_offset);

	_stream.write ((char*)_darr,sizeof(double)*_size);

//	_stream.close ();

//	_stream.open ();

	_stream.flush ();

};

// Author: Kharchenko S.A.
// Description: Read double array (direct access)
// FGet()
//========================================================================================
void FGet (fstream &_stream, int _size, double *_darr, int _offset) { // Read double array (direct access)

//	_stream.seekg(sizeof(double)*_offset,ios::beg);
//	_stream.seekg(sizeof(double)*_offset);

	int i = _offset;
	i = i*i;

	_stream.read ((char*)_darr,sizeof(double)*_size);

};

// Author: Kharchenko S.A.
// Description: Open the set of direct IO files according to the path, name and extention
// OpenDIOFiles()
//========================================================================================
void OpenDIOFiles (int _nfiles, fstream *_farr, // Open the set of direct IO files according to the path, name and extension
					const char *_path, const char *_name, const char *_extention) {

	char strbuff[64];

	for (int ifile=0;ifile<_nfiles;ifile++) {

		sprintf(strbuff, "%s%s%d%s",_path,_name,ifile,_extention);

//		_farr[ifile].open (strbuff, ios::in | ios::out | ios::binary | ios::trunc);
		_farr[ifile].open (strbuff, ios::in | ios::out | ios::binary | ios::app | ios::ate);

		if (!_farr[ifile]) {
			std::cout << " Error when opening direct IO file named " << strbuff << std::endl;
			throw " Error when opening direct IO file";
		};

	};

};

// Author: Kharchenko S.A.
// Description: Close the set of direct IO files
// CloseDIOFiles()
//========================================================================================
void CloseDIOFiles (int _nfiles, fstream *_farr) { // Close the set of direct IO files

	for (int ifile=0;ifile<_nfiles;ifile++) {

		_farr[ifile].close();

	};

};

// Author: Kharchenko S.A.
// Description: Store double array (direct access)
// FPutG()
//========================================================================================
void FPutG (FILE *_file, int _size, void *_varr, int _offset) { // Store double array (direct access)

	size_t i;

	i = fseek(_file,_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fwrite(_varr,1,_size,_file);
	if (i != _size) {
		cout << " Error when writing direct IO file " << endl;
		throw "Error when writing direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Store integer array (direct access)
// FPut()
//========================================================================================
void FPut (FILE *_file, int _size, int *_iarr, int _offset) { // Store integer array (direct access)

	const char *funcname = "FPut_00";

	int i;

	double *darr;

	darr = new double [_size];
	if (!darr) MemoryFail (funcname);

	for (i=0;i<_size;i++) darr[i] = (double) _iarr[i];

	FPut (_file, _size, darr, _offset);

	delete [] darr;

};

// Author: Kharchenko S.A.
// Description: Store double array (direct access)
// FPut()
//========================================================================================
void FPut (FILE *_file, int _size, double *_darr, int _offset) { // Store double array (direct access)

	size_t i;
   
	i = fseek(_file,sizeof(double)*_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fwrite(_darr,sizeof(double),_size,_file);
	if (i != _size) {
		cout << " Error when writing direct IO file " << endl;
		throw "Error when writing direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Store double complex array (direct access)
// FPut()
//========================================================================================
void FPut (FILE *_file, int _size, dcmplx *_darr, int _offset) { // Store double complex array (direct access)

	size_t i;

	i = fseek(_file,sizeof(dcmplx)*_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fwrite(_darr,sizeof(dcmplx),_size,_file);
	if (i != _size) {
		cout << " Error when writing direct IO file " << endl;
		throw "Error when writing direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Read void array (direct access)
// FGetG()
//========================================================================================
void FGetG (FILE *_file, int _size, void *_varr, int _offset) { // Read void array (direct access)

	size_t i;

	i = fseek(_file,_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fread(_varr,1,_size,_file);
	if (i != _size) {
		cout << " Error when reading from direct IO file " << endl;
		throw "Error when reading from direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Read double array (direct access)
// FGet()
//========================================================================================
void FGet (FILE *_file, int _size, double *_darr, int _offset) { // Read double array (direct access)

	size_t i;

	i = fseek(_file,sizeof(double)*_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fread(_darr,sizeof(double),_size,_file);
	if (i != _size) {
		cout << " Error when reading from direct IO file " << endl;
		throw "Error when reading from direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Read double complex array (direct access)
// FGet()
//========================================================================================
void FGet (FILE *_file, int _size, dcmplx *_darr, int _offset) { // Read double complex array (direct access)

	size_t i;

	i = fseek(_file,sizeof(dcmplx)*_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fread(_darr,sizeof(dcmplx),_size,_file);
	if (i != _size) {
		cout << " Error when reading from direct IO file " << endl;
		throw "Error when reading from direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Read char array (direct access)
// FGet()
//========================================================================================
void FGet (FILE *_file, int _size, char *_charr, int _offset) { // Read char array (direct access)

	size_t i;

	i = fseek(_file,sizeof(char)*_offset,SEEK_SET);
	if (i != 0) {
		cout << " Error when moving cursor for direct IO file " << endl;
		throw "Error when moving cursor for direct IO file";
	};

	i = fread(_charr,sizeof(char),_size,_file);
	if (i != _size) {
		cout << " Error when reading from direct IO file " << endl;
		throw "Error when reading from direct IO file";
	};

};

// Author: Kharchenko S.A.
// Description: Open the set of direct IO files according to the path, name and extention
// OpenDIOFiles()
//========================================================================================
void OpenDIOFiles (int _nfiles, FILE **_farr, // Open the set of direct IO files according to the path, name and extension
					const char *_path, const char *_name, const char *_extention) {

	char strbuff[512];

	for (int ifile=0;ifile<_nfiles;ifile++) {

		sprintf(strbuff, "%s%s%d%s",_path,_name,ifile,_extention);

		_farr[ifile] = fopen (strbuff, "r+b");

		if (!_farr[ifile]) {

			_farr[ifile] = fopen (strbuff, "w+b");

			if (!_farr[ifile]) {
				cout << " Error when opening direct IO file named " << strbuff << endl;
				throw " Error when opening direct IO file";
			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Open the set of direct IO files according to the path, name and extention and store file names
// OpenDIOFiles()
//========================================================================================
void OpenDIOFiles (int _nfiles, FILE **_farr, char **_fnames, // Open the set of direct IO files according to the path, name and extension and store file names
					const char *_path, const char *_name, const char *_extention) {

	char strbuff[512];

	for (int ifile=0;ifile<_nfiles;ifile++) {

		sprintf(strbuff, "%s%s%d%s",_path,_name,ifile,_extention);

		strcpy (_fnames[ifile], strbuff);

		_farr[ifile] = fopen (strbuff, "r+b");

		if (!_farr[ifile]) {

			_farr[ifile] = fopen (strbuff, "w+b");

			if (!_farr[ifile]) {
				cout << " Error when opening direct IO file named " << strbuff << endl;
				throw " Error when opening direct IO file";
			};

		};

	};

};

// Author: Kharchenko S.A.
// Description: Flush the set of direct IO files
// FlushDIOFiles()
//========================================================================================
void FlushDIOFiles (int _nfiles, FILE **_farr) { // Flush the set of direct IO files

	for (int ifile=0;ifile<_nfiles;ifile++) {

		if (fflush(_farr[ifile])) {
			cout << " Error when flushing direct IO file" << endl;
			throw " Error when flushing direct IO file";
		};

	};

};

// Author: Kharchenko S.A.
// Description: Close the set of direct IO files
// CloseDIOFiles()
//========================================================================================
void CloseDIOFiles (int _nfiles, FILE **_farr) { // Close the set of direct IO files

	for (int ifile=0;ifile<_nfiles;ifile++) {

		if (fclose(_farr[ifile])) {
			cout << " Error when closing direct IO file" << endl;
			throw " Error when closing direct IO file";
		};

	};

};

// Author: Kharchenko S.A.
// Description: Close the set of direct IO files and remove them
// CloseDIOFiles()
//========================================================================================
void CloseDIOFiles (int _nfiles, FILE **_farr, char **_fnames) { // Close the set of direct IO files and remove them

	for (int ifile=0;ifile<_nfiles;ifile++) {

		if (fclose(_farr[ifile])) {
			cout << " Error when closing direct IO file" << endl;
			throw " Error when closing direct IO file";
		};

		remove(_fnames[ifile]);

//		if (!remove(_fnames[ifile])) {
//			cout << " Error when removing direct IO file" << endl;
//			throw " Error when removing direct IO file";
//		};

	};

};

// Author: Kharchenko S.A.
// Description: Read prescribed symbol from file
// ReadSymbol()
//========================================================================================
int ReadSymbol (char *_file, int &_k, int _fend, char _symbol) { // Read prescribed symbol from file

	char buf = _file[_k];

	while (buf != _symbol && buf != '\x0') {
//		_k++;
//		if (_k < _fend) {
//			buf = _file[_k];
//		} else {
//			buf = '\x0';
//			_k--;
//		};
		if (_k >= _fend-1) {
			_k++;
			break;
		};
		buf = _file[++_k];
	};

	if (buf == '\x0') return(0); 
	else return (1);

};

// Author: Kharchenko S.A.
// Description: Write prescribed symbol from file
// WriteSymbol()
//========================================================================================
void WriteSymbol (char _symbol, char *_file, int &_k, int _fend) { // Write prescribed symbol into file

	if (_k+1<_fend) _file[_k++] = _symbol;

};

// Author: Kharchenko S.A.
// Description: Read current char from file
// ReadChar()
//========================================================================================
char ReadChar (char *_file, int &_k, int _fend) { // Read current char from file
	if (_k >= _fend) return('\xFF');
	char buf = toupper(_file[_k++]);
	if (buf == ';') {
		while (buf != '\n' && _k < _fend) buf = _file[_k++];
		return('\xFF');
	};
	if (isalpha(buf) || isdigit(buf) || buf == '+' || buf == ':' && (_k == 1 || isalpha(_file[_k-2])) || 
		buf == '-'  || buf == '.' || buf == '#' || buf == '_' || 
		buf == '\\' || buf == '!' || buf == '*') {
		return(buf);
	} else return ('\xFF');
};

// Author: Kharchenko S.A.
// Description: Read current char from file
// ReadCharAsIs()
//========================================================================================
char ReadCharAsIs (char *_file, int &_k, int _fend) { // Read current char from file
	if (_k >= _fend) return('\xFF');
	char buf = _file[_k++];
	if (buf == ';') {
		while (buf != '\n' && _k < _fend) buf = _file[_k++];
		return('\xFF');
	};
	if (isalpha(buf) || isdigit(buf) || buf == '+' || buf == ':' && (_k == 1 || isalpha(_file[_k-2])) || 
		buf == '-'  || buf == '.' || buf == '#' || buf == '_' || 
		buf == '\\' || buf == '!') {
		return(buf);
	} else return ('\xFF');
};

// Author: Kharchenko S.A.
// Description: Read current word from file according to length
// ReadWordLen()
//========================================================================================
int ReadWordLen (char *_buf, char *_file, int &_k, int _fend, int _m) { // Read current word from file according to length
	int l = 0;
	for (_m--; _m >= 0; _m--) {
		if (l < CHAR_BUF_SIZE && _k < _fend) {
			_buf[l++] = toupper(_file[_k++]);
		};
	};
	_buf[l] = '\x0'; 
	return (l != 0);
};

// Author: Kharchenko S.A.
// Description: Write current word into file according to length
// WriteWordLen()
//========================================================================================
void WriteWordLen (const char *_buf, char *_file, int &_k, int _fend, int _m) { // Write current word from file according to length
	int lenloc = (int)strlen(_buf);
	int i;
	for (i=0;i<lenloc;i++) {
		if (_k+i<_fend && i < _m) {
			_file[_k++] = _buf[i];
		};
	};
	for (i=lenloc;i<_m;i++) {
		if (_k+i<_fend) {
			_file[_k++] = ' ';
		};
	};
};

// Author: Kharchenko S.A.
// Description: Write int word into file according to length
// WriteIntLen()
//========================================================================================
void WriteIntLen (int _j, char *_file, int &_k, int _fend, int _m) { // Write int word from file according to length
	char buf[CHAR_BUF_SIZE];
	sprintf(buf,"%8i",_j);
	WriteWordLen (buf, _file, _k, _fend, _m);
};

// Author: Kharchenko S.A.
// Description: Write double word into file according to length
// WriteDblLen()
//========================================================================================
void WriteDblLen (double _f, char *_file, int &_k, int _fend, int _m) { // Write int word from file according to length
	char buf[CHAR_BUF_SIZE];
	if (_m == 8) {
//		sprintf(buf,"%8.5g",_f);

		char pb1[64];
		char pb2[64];

		int e = (int) log10(fabs(_f));

		double m = _f / pow(10.0, e);

		sprintf(pb1, "%8.8f", m);

		if (e != 0)
			sprintf(pb2, "%c%d", 'e', e);
		else pb2[0] = 0;

		size_t mlen = _m - strlen(pb2);
		pb1[mlen] = 0;

		sprintf(buf, "%*s%s", (int)mlen, pb1, pb2);

	} else {
		char pb1[64];
		char pb2[64];

		int e = (int) log10(fabs(_f));

		double m = _f / pow(10.0, e);

		sprintf(pb1, "%16.16f", m);

		if (e != 0)
			sprintf(pb2, "%c%d", 'e', e);
		else pb2[0] = 0;

		size_t mlen = _m - strlen(pb2);
		pb1[mlen] = 0;

		sprintf(buf, "%*s%s", (int)mlen, pb1, pb2);

	};
	WriteWordLen (buf, _file, _k, _fend, _m);
};

// Author: Kharchenko S.A.
// Description: Read current word from file
// ReadWord()
//========================================================================================
int ReadWord (char *_buf, char *_file, int &_k, int _fend) { // Read current word from file
	int l = 0;
	while (l < CHAR_BUF_SIZE && _k < _fend) {
		char buf0 = ReadChar (_file, _k, _fend);
		if (buf0 != '\xFF') {
			_buf[l++] = buf0; 
		} else {
			if (l != 0) {
				_buf[l] = '\x0'; 
				return(1);
			};
		};
	};
	_buf[l] = '\x0';
	return (l != 0);
};

// Author: Kharchenko S.A.
// Description: Read current word from file
// ReadWordAsIs()
//========================================================================================
int ReadWordAsIs (char *_buf, char *_file, int &_k, int _fend) { // Read current word from file
	int l = 0;
	while (l < CHAR_BUF_SIZE && _k < _fend) {
		char buf0 = ReadCharAsIs (_file, _k, _fend);
		if (buf0 != '\xFF') {
			_buf[l++] = buf0; 
		} else {
			if (l != 0) {
				_buf[l] = '\x0'; 
				return(1);
			};
		};
	};
	_buf[l] = '\x0';
	return (l != 0);
};

// Author: Kharchenko S.A.
// Description: Read current word from file
// ReadWordInBrakets()
//========================================================================================
int ReadWordInBrakets (char *_buf, char *_file, int &_k, int _fend) { // Read current word in brakets from file
	int fendloc = _fend;
	int i = ReadSymbol (_file, _k, fendloc, '"');
	if (i==0) return i;
	_k++;
	int kloc = _k;
	fendloc = _fend;
	i = ReadSymbol (_file, _k, fendloc, '"');
	if (i==0) return i;
	_k++;
	int l = 0;
	while (l < CHAR_BUF_SIZE && kloc < _k-1) {
		_buf[l++] = _file[kloc++];
	};
	_buf[l] = '\x0';
	return (l != 0);
};

// Author: Kharchenko S.A.
// Description: Determine the length of ascii file
// LengthOfAsciiFile()
//========================================================================================
int LengthOfAsciiFile (char *_filename) { // Determine the length of ascii file
	int length = 0;
	FILE * BLN;
	if ((BLN = fopen(_filename, "rb")) != NULL) {
		fseek (BLN, 0L, SEEK_END);
		length = (unsigned long ) ftell(BLN);
		fclose(BLN);
	};
	return(length);
}

// Author: Kharchenko S.A.
// Description: Read current file into char variable
// ReadAsciiFile()
//========================================================================================
char * ReadAsciiFile (char *_filename, int &_length) { // Read current file into char variable

	const char *funcname = "ReadAsciiFile";

	_length = 0;
	_length = LengthOfAsciiFile (_filename);

	cout << " Read from Ascii file {" << _filename << "}" << " Length = " << _length << endl;

	char *filechar = 0;
	if (_length) {
		filechar = new char [_length+1];
		if (!filechar) MemoryFail (funcname);
		if (filechar != NULL) {
			strcpy (filechar,"");
			FILE *file = fopen(_filename, "rb");
			if (ferror(file)) {
				char str[256] = "";
				perror (str);
				cout << " Error in open = " << str << endl;
				throw " ReadAsciiFile: error in open";
			};
			FGet (file, _length, filechar, 0);
			fclose(file);
		};
	};
	cout << " Reading completed successfully" << endl;
	return (filechar);
};

// Author: Kharchenko S.A.
// Description: Remove comments
// RemoveComments()
//========================================================================================
char * RemoveComments (char *_file, int &_length) { // Remove comments

	const char *funcname = "RemoveComments";

// Scan current file and count the maximal number of commented rows

	int ncomments = 0;

	int fbeg = 0;

	int fend = _length;

	while (fbeg < _length && ReadSymbol (_file, fbeg, fend, '/')) {
		if (fbeg+1 < _length) {
			fbeg++;
			if (_file[fbeg] == '/') {
				ncomments++;
			};
			fbeg++;
		};
	};

	cout << " The maximal number of commented lines =" << ncomments << endl;

// Store indices of commented rows

	int *comment2ind;

	comment2ind = new int [2*ncomments];
	if (!comment2ind) MemoryFail (funcname);

	fbeg = 0;
	fend = _length;

	char buf[CHAR_BUF_SIZE+1];

	ncomments = 0;

	while (fbeg < _length && ReadSymbol (_file, fbeg, fend, '/')) {
		if (fbeg+1 < _length) {
			fbeg++;
			if (_file[fbeg] == '/') {
				comment2ind[ncomments*2] = fbeg-1;
				char cbuf = _file[fbeg];
				int k = fbeg;
				while (cbuf != '\r' && cbuf != '\n' && buf != '\x0') {
					cbuf = _file[++k];
				};
				comment2ind[ncomments*2+1] = k;
				ncomments++;
				fbeg = k-1;
			};
			fbeg++;
		};
	};

// Remove comments from current file

	char *filechar;

	filechar = new char [_length+1];
	if (!filechar) MemoryFail (funcname);

	int i, icomment, nz=0;

	if (ncomments == 0) {
		for (i=0;i<_length;i++) filechar[nz++] = _file[i];
	} else {
		for (i=0;i<comment2ind[0];i++) filechar[nz++] = _file[i];
		for (icomment=0;icomment<ncomments-1;icomment++) {
			for (i=comment2ind[icomment*2+1];i<comment2ind[icomment*2+2];i++) filechar[nz++] = _file[i];
		};
		for (i=comment2ind[ncomments*2-1];i<_length;i++) filechar[nz++] = _file[i];
	};

	filechar[nz] = '\x0';

	_length = nz-1;

// Free work arrays

	delete [] comment2ind;

	cout << " Comments are removed from file" << endl;

	return (filechar);

};

// Author: Kharchenko S.A.
// Description: Add two sorted lists
// AddLists()
//========================================================================================
void AddLists (int _nlist1, int *_iarr1, // Add two sorted lists
				int _nlist2, int *_iarr2,
				int &_nlist, int *_iarr) { 

	int i1, i2, j1, j2;

	_nlist = 0;

	i1 = 0;
	i2 = 0;

	while (i1 < _nlist1 || i2 < _nlist2) {
		if (i1 < _nlist1 && i2 < _nlist2) {
			j1 = _iarr1[i1];
			j2 = _iarr2[i2];
			if (j1 == j2) {
				_iarr[_nlist] = j1;
				i1++;
				i2++;
				_nlist++;
			} else if (j1 < j2) {
				_iarr[_nlist] = j1;
				i1++;
				_nlist++;
			} else if (j1 > j2) {
				_iarr[_nlist] = j2;
				i2++;
				_nlist++;
			};
		} else if (i1 < _nlist1) {
			j1 = _iarr1[i1];
			_iarr[_nlist] = j1;
			i1++;
			_nlist++;
		} else if (i2 < _nlist2) {
			j2 = _iarr2[i2];
			_iarr[_nlist] = j2;
			i2++;
			_nlist++;
		};
	};

};

// Author: Kharchenko S.A.
// Description: Add two sorted 2 index lists
// AddLists()
//========================================================================================
void AddLists (int _nlist1, int *_iarr1, int *_iarrblk1, // Add two sorted 2 index lists
				int _nlist2, int *_iarr2, int *_iarrblk2,
				int &_nlist, int *_iarr, int *_iarrblk) { 

	int i1, i2, j1, j2, jblk1, jblk2;

	_nlist = 0;

	i1 = 0;
	i2 = 0;

	while (i1 < _nlist1 || i2 < _nlist2) {
		if (i1 < _nlist1 && i2 < _nlist2) {
			j1 = _iarr1[i1];
			jblk1 = _iarrblk1[i1];
			j2 = _iarr2[i2];
			jblk2 = _iarrblk2[i2];
			if (jblk1 == jblk2 && j1 == j2) {
				_iarr[_nlist] = j1;
				_iarrblk[_nlist] = jblk1;
				i1++;
				i2++;
				_nlist++;
			} else if (jblk1 < jblk2 || (jblk1 == jblk2 && j1 < j2)) {
				_iarr[_nlist] = j1;
				_iarrblk[_nlist] = jblk1;
				i1++;
				_nlist++;
			} else if (jblk2 < jblk1 || (jblk1 == jblk2 && j2 < j1)) {
				_iarr[_nlist] = j2;
				_iarrblk[_nlist] = jblk2;
				i2++;
				_nlist++;
			};
		} else if (i1 < _nlist1) {
			j1 = _iarr1[i1];
			jblk1 = _iarrblk1[i1];
			_iarr[_nlist] = j1;
			_iarrblk[_nlist] = jblk1;
			i1++;
			_nlist++;
		} else if (i2 < _nlist2) {
			j2 = _iarr2[i2];
			jblk2 = _iarrblk2[i2];
			_iarr[_nlist] = j2;
			_iarrblk[_nlist] = jblk2;
			i2++;
			_nlist++;
		};
	};

};

// Author: Kharchenko S.A.
// Description: Sort elements of the array
// SortElm()
//========================================================================================
void SortElm (int _n, int *_iarr, int *_iarr2) { // Sort elements of the array
	int k, ii, top, bottom, iaux;
	int aux;

// Init variables

	k=0;
	bottom = 0;
	top = _n-2;

// Main cycle

cycle:

	if (bottom <= top) {
		for (ii=bottom;ii<=top;ii++) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _iarr2[ii];
				_iarr2[ii] = _iarr2[ii+1];
				_iarr2[ii+1] = aux;
				k = ii;
			};
		};
		top = k-1;
		for (ii=top;ii>=bottom;ii--) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _iarr2[ii];
				_iarr2[ii] = _iarr2[ii+1];
				_iarr2[ii+1] = aux;
				k = ii;
			};
		};
		bottom = k+1;
		goto cycle;
	};
};

// Author: Kharchenko S.A.
// Description: Sort elements of the array
// SortElm()
//========================================================================================
void SortElm (int _n, int *_iarr, double *_darr) { // Sort elements of the array
	int k, ii, top, bottom, iaux;
	double aux;

// Init variables

	k=0;
	bottom = 0;
	top = _n-2;

// Main cycle

cycle:

	if (bottom <= top) {
		for (ii=bottom;ii<=top;ii++) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _darr[ii];
				_darr[ii] = _darr[ii+1];
				_darr[ii+1] = aux;
				k = ii;
			};
		};
		top = k-1;
		for (ii=top;ii>=bottom;ii--) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _darr[ii];
				_darr[ii] = _darr[ii+1];
				_darr[ii+1] = aux;
				k = ii;
			};
		};
		bottom = k+1;
		goto cycle;
	};
};

// Author: Kharchenko S.A.
// Description: Sort pointers to the double array
// SortPtr()
//========================================================================================
void SortPtr (int _n, int *_iarr, double **_pdarr) { // Sort pointers to the double array
	int k, ii, top, bottom, iaux;
	double *aux;

// Init variables

	k=0;
	bottom = 0;
	top = _n-2;

// Main cycle

cycle:

	if (bottom <= top) {
		for (ii=bottom;ii<=top;ii++) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _pdarr[ii];
				_pdarr[ii] = _pdarr[ii+1];
				_pdarr[ii+1] = aux;
				k = ii;
			};
		};
		top = k-1;
		for (ii=top;ii>=bottom;ii--) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _pdarr[ii];
				_pdarr[ii] = _pdarr[ii+1];
				_pdarr[ii+1] = aux;
				k = ii;
			};
		};
		bottom = k+1;
		goto cycle;
	};
};

// Author: Kharchenko S.A.
// Description: Sort pointers to the double complex array
// SortPtr()
//========================================================================================
void SortPtr (int _n, int *_iarr, dcmplx **_pdarr) { // Sort pointers to the double complex array
	int k, ii, top, bottom, iaux;
	dcmplx *aux;

// Init variables

	k=0;
	bottom = 0;
	top = _n-2;

// Main cycle

cycle:

	if (bottom <= top) {
		for (ii=bottom;ii<=top;ii++) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _pdarr[ii];
				_pdarr[ii] = _pdarr[ii+1];
				_pdarr[ii+1] = aux;
				k = ii;
			};
		};
		top = k-1;
		for (ii=top;ii>=bottom;ii--) {
			if (_iarr[ii] > _iarr[ii+1]) {
				iaux = _iarr[ii];
				_iarr[ii] = _iarr[ii+1];
				_iarr[ii+1] = iaux;
				aux = _pdarr[ii];
				_pdarr[ii] = _pdarr[ii+1];
				_pdarr[ii+1] = aux;
				k = ii;
			};
		};
		bottom = k+1;
		goto cycle;
	};
};

// Author: Kharchenko S.A.
// Description: Round two double values to three digits
// Round()
//========================================================================================
void Round (double &dx, double &dy) { // Round two double values to three digits
	int i, i0;
	i0 = (int) dx;
	i = (int) ((dx-(double)i0) * 1000);
	dx = (double) i;
	dx = dx * 1.0e-3+(double)i0;
	i0 = (int) dy;
	i = (int) ((dy-(double)i0) * 1000);
	dy = (double) i;
	dy = dy * 1.0e-3+(double)i0;
};

// Author: Kharchenko S.A.
// Description: Setup output manipulator
// SetPw()
//========================================================================================
ostream &SetPw (ostream &stream) { // Setup output manipulator
	stream.unsetf(ios::scientific);
//	stream << ' ' << left << setprecision(5) << setw (9);
	stream << ' ' << setprecision(5) << setw (9);
	return stream;
};

// Author: Kharchenko S.A.
// Description: Compute the list of blocks to be checked
// ListBlocks3D()
//========================================================================================
void ListBlocks3D (int _nelem, double *_xyzrelem, int _nballs, double *_xyzrballs, // Compute the list of blocks to be checked
					int *&_ilinks, int *&_jlinks) {

// Compute the size of the integer 3D mesh

	int ncube;

	double daux, daux1;

	daux = pow ((double) _nelem,1.0e0/3.0e0);

	ncube = (int) daux;

	if (ncube < 1) ncube = 1;

// Compute the x, y and z min/max values

	int ix=0, iy=1, iz=2, ir=3;

	double xmin, xmax, ymin, ymax, zmin, zmax;
	double xminl, xmaxl, yminl, ymaxl, zminl, zmaxl;

	int ielem, iball, i, j, k, ii, jj;

	ielem = 0;

	xmin = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
	xmax = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];

	ymin = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
	ymax = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];

	zmin = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
	zmax = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];

	for (ielem=1;ielem<_nelem;ielem++) {
		daux = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
		if (daux < xmin) xmin = daux;
		daux = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];
		if (daux > xmax) xmax = daux;
		daux = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
		if (daux < ymin) ymin = daux;
		daux = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];
		if (daux > ymax) ymax = daux;
		daux = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
		if (daux < zmin) zmin = daux;
		daux = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];
		if (daux > zmax) zmax = daux;
	};

// Set x, y and z linear interpolation data

	double deltax = (xmax-xmin) / ((double) ncube);
	double deltay = (ymax-ymin) / ((double) ncube);
	double deltaz = (zmax-zmin) / ((double) ncube);

// For each subcube count the number of nodes

	int ncube3 = ncube * ncube * ncube;

	int ncube2 = ncube * ncube;

	int *icube, *jcube, *iptr;

	icube = new int [ncube3+1];
	iptr = new int [ncube3];

	for (i=0;i<ncube3+1;i++) icube[i] = 0;

	int ixmin, ixmax, iymin, iymax, izmin, izmax;

	for (ielem=0;ielem<_nelem;ielem++) {

// For current element determine the set of 3D indices that contains the element

		daux = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmin = (int) daux1;
		ixmin--;
		if (ixmin < 0) ixmin = 0;
		if (ixmin >= ncube) ixmin = ncube-1;

		daux = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmax = (int) daux1;
		ixmax++;
		if (ixmax < 0) ixmax = 0;
		if (ixmax >= ncube) ixmax = ncube-1;

		daux = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymin = (int) daux1;
		iymin--;
		if (iymin < 0) iymin = 0;
		if (iymin >= ncube) iymin = ncube-1;

		daux = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymax = (int) daux1;
		iymax++;
		if (iymax < 0) iymax = 0;
		if (iymax >= ncube) iymax = ncube-1;

		daux = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmin = (int) daux1;
		izmin--;
		if (izmin < 0) izmin = 0;
		if (izmin >= ncube) izmin = ncube-1;

		daux = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmax = (int) daux1;
		izmax++;
		if (izmax < 0) izmax = 0;
		if (izmax >= ncube) izmax = ncube-1;

// Register that all detected subcubes contain current element

		for (i=ixmin;i<=ixmax;i++) {
			for (j=iymin;j<=iymax;j++) {
				for (k=izmin;k<=izmax;k++) {
					icube[i*ncube2+j*ncube+k+1]++;
				};
			};
		};

	};

	for (i=0;i<ncube3;i++) icube[i+1] += icube[i];

// Complete computation of the lists

	int nlist = icube[ncube3];

	jcube = new int [nlist];

	for (i=0;i<ncube3;i++) iptr[i] = icube[i];

	for (ielem=0;ielem<_nelem;ielem++) {

// For current element again determine the set of 3D indices that contains the element

		daux = _xyzrelem[ielem*4+ix] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmin = (int) daux1;
		ixmin--;
		if (ixmin < 0) ixmin = 0;
		if (ixmin >= ncube) ixmin = ncube-1;

		daux = _xyzrelem[ielem*4+ix] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-xmin) / deltax;
		ixmax = (int) daux1;
		ixmax++;
		if (ixmax < 0) ixmax = 0;
		if (ixmax >= ncube) ixmax = ncube-1;

		daux = _xyzrelem[ielem*4+iy] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymin = (int) daux1;
		iymin--;
		if (iymin < 0) iymin = 0;
		if (iymin >= ncube) iymin = ncube-1;

		daux = _xyzrelem[ielem*4+iy] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-ymin) / deltay;
		iymax = (int) daux1;
		iymax++;
		if (iymax < 0) iymax = 0;
		if (iymax >= ncube) iymax = ncube-1;

		daux = _xyzrelem[ielem*4+iz] - _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmin = (int) daux1;
		izmin--;
		if (izmin < 0) izmin = 0;
		if (izmin >= ncube) izmin = ncube-1;

		daux = _xyzrelem[ielem*4+iz] + _xyzrelem[ielem*4+ir];
		daux1 = (daux-zmin) / deltaz;
		izmax = (int) daux1;
		izmax++;
		if (izmax < 0) izmax = 0;
		if (izmax >= ncube) izmax = ncube-1;

// Register that all detected subcubes contain current element

		for (i=ixmin;i<=ixmax;i++) {
			for (j=iymin;j<=iymax;j++) {
				for (k=izmin;k<=izmax;k++) {
					int ip = iptr[i*ncube2+j*ncube+k];
					jcube[ip] = ielem;
					iptr[i*ncube2+j*ncube+k]++;
				};
			};
		};

	};

// Prepare mask arrays

	int *imask, *list;

	imask = new int [_nelem];
	list = new int [_nelem];

	int icycle = -1;

	for (i=0;i<_nelem;i++) imask[i] = icycle;

// For each ball compute the set of elements

	_ilinks = new int [_nballs+1];

	int **jlinksarr;

	jlinksarr = new int * [_nballs];

	_ilinks[0] = 0;

	for (iball=0;iball<_nballs;iball++) {

// For current ball determine the set of 3D indices that contains the ball

		xminl = _xyzrballs[iball*4+ix] - _xyzrballs[iball*4+ir];
		daux = xminl;
		daux1 = (daux-xmin) / deltax;
		ixmin = (int) daux1;
		ixmin--;
		if (ixmin < 0) ixmin = 0;
		if (ixmin >= ncube) ixmin = ncube-1;

		xmaxl = _xyzrballs[iball*4+ix] + _xyzrballs[iball*4+ir];
		daux = xmaxl;
		daux1 = (daux-xmin) / deltax;
		ixmax = (int) daux1;
		ixmax++;
		if (ixmax < 0) ixmax = 0;
		if (ixmax >= ncube) ixmax = ncube-1;

		yminl = _xyzrballs[iball*4+iy] - _xyzrballs[iball*4+ir];
		daux = yminl;
		daux1 = (daux-ymin) / deltay;
		iymin = (int) daux1;
		iymin--;
		if (iymin < 0) iymin = 0;
		if (iymin >= ncube) iymin = ncube-1;

		ymaxl = _xyzrballs[iball*4+iy] + _xyzrballs[iball*4+ir];
		daux = ymaxl;
		daux1 = (daux-ymin) / deltay;
		iymax = (int) daux1;
		iymax++;
		if (iymax < 0) iymax = 0;
		if (iymax >= ncube) iymax = ncube-1;

		zminl = _xyzrballs[iball*4+iz] - _xyzrballs[iball*4+ir];
		daux = zminl;
		daux1 = (daux-zmin) / deltaz;
		izmin = (int) daux1;
		izmin--;
		if (izmin < 0) izmin = 0;
		if (izmin >= ncube) izmin = ncube-1;

		zmaxl = _xyzrballs[iball*4+iz] + _xyzrballs[iball*4+ir];
		daux = zmaxl;
		daux1 = (daux-zmin) / deltaz;
		izmax = (int) daux1;
		izmax++;
		if (izmax < 0) izmax = 0;
		if (izmax >= ncube) izmax = ncube-1;

// Scan the list of elements in each subcube

		icycle++;

		nlist = 0;

		for (i=ixmin;i<=ixmax;i++) {
			for (j=iymin;j<=iymax;j++) {
				for (k=izmin;k<=izmax;k++) {
					int ind = i*ncube2+j*ncube+k;
					for (ii=icube[ind];ii<icube[ind+1];ii++) {
						jj = jcube[ii];
						if (imask[jj] != icycle) {

							imask[jj] = icycle;

							daux = _xyzrelem[jj*4+ix] - _xyzrelem[jj*4+ir];
							if (daux > xmaxl) goto ExitElem;
							daux = _xyzrelem[jj*4+ix] + _xyzrelem[jj*4+ir];
							if (daux < xminl) goto ExitElem;
							daux = _xyzrelem[jj*4+iy] - _xyzrelem[jj*4+ir];
							if (daux > ymaxl) goto ExitElem;
							daux = _xyzrelem[jj*4+iy] + _xyzrelem[jj*4+ir];
							if (daux < yminl) goto ExitElem;
							daux = _xyzrelem[jj*4+iz] - _xyzrelem[jj*4+ir];
							if (daux > zmaxl) goto ExitElem;
							daux = _xyzrelem[jj*4+iz] + _xyzrelem[jj*4+ir];
							if (daux < zminl) goto ExitElem;

							list[nlist] = jj;
							nlist++;

ExitElem:;

						};
					};
				};
			};
		};

		qsort (list, nlist, sizeof(int), compint);

		int *jlinkloc;

		jlinkloc = new int [nlist];

		for (i=0;i<nlist;i++) jlinkloc[i] = list[i];

		jlinksarr[iball] = jlinkloc;

		_ilinks[iball+1] = _ilinks[iball] + nlist;

	};

// Transform data

	int nz = _ilinks[_nballs];

	_jlinks = new int [nz];

	for (iball=0;iball<_nballs;iball++) {
		int nzloc = _ilinks[iball+1]-_ilinks[iball];
		int ibs = _ilinks[iball];
		int *jlinkloc = jlinksarr[iball];
		for (i=0;i<nzloc;i++) _jlinks[ibs+i] = jlinkloc[i];
	};

// Free work arrays

	delete [] icube;
	delete [] iptr;
	delete [] jcube;
	delete [] imask;
	delete [] list;
	for (i=0;i<_nballs;i++) {
		delete [] jlinksarr[i];
	};
	delete [] jlinksarr;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
