#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>

#include "globals.h"

using namespace std;

// Preliminary declaration of the routines

extern "C" {

void dlarfgdum_ (int *n, double *alpha, double *x, 
				int *incx, double *tau);

void dsyevdum_ (int *itype, int *n, double *a,
				int *lda, double *w, double *work, int *lwork, int *info);

void dgesvddum_ (int *itype, int *m, int *n, 
				double *a, int *lda, double *s, double *u, int *ldu, 
				double *vt, int *ldvt, double *work, int *lwork, int *info);

void zlarfgdum_ (int *n, dcmplx *alpha, dcmplx *x, 
				int *incx, dcmplx *tau);

void zheevdum_ (int *itype, int *n, dcmplx *a,
				int *lda, double *w, dcmplx *work, int *lwork, double *rwork, int *info);

void zgesvddum_ (int *itype, int *m, int *n, 
				dcmplx *a, int *lda, double *s, dcmplx *u, 
				int *ldu, dcmplx *vt, int *ldvt, dcmplx *work, 
				int *lwork, double *rwork, int *info);

};

// Author: Kharchenko S.A.
// Test Lapack functions that compute SVD and QR
//========================================================================================
void TestLapack () { // Test Lapack functions that compute SVD and QR

	const char *funcname = "TestLapack";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

	int m = 550;
	int n = 150;

// Allocate real memory 

	int lwork = 10*n;

	double *ar;
	double *ur;
	double *vr;
	double *zr;
	double *sv;
	double *workr;

	ar = new double [m*n];
	if (!ar) MemoryFail (funcname);
	ur = new double [m*n];
	if (!ur) MemoryFail (funcname);
	vr = new double [n*n];
	if (!vr) MemoryFail (funcname);
	zr = new double [m*n];
	if (!zr) MemoryFail (funcname);
	sv = new double [n];
	if (!sv) MemoryFail (funcname);
	workr = new double [lwork];
	if (!workr) MemoryFail (funcname);

// Test the case of double data

// Check QR decomposition

	int kii, kjj, kkk;

	for (kii=0;kii<m;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ar[kjj*m+kii] = 1.0e0 / ((float) (kii-kjj) + 0.5e0);
			zr[kjj*m+kii] = ar[kjj*m+kii];
		};
	};

	QrdBlk (m, n, ar, sv);

	for (kii=0;kii<n*n;kii++) vr[kii] = 0.0e0;
	for (kii=0;kii<n;kii++) vr[kii*n+kii] = 1.0e0;

	MvmQBlk (m, n, n,
				ar, sv,
				vr, ur);

	double daux;

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			daux = 0.0e0;
			for (kkk=0;kkk<m;kkk++) {
				daux += ur[kii*m+kkk] * ur[kjj*m+kkk];
			};
			vr[kii*n+kjj] = daux;
		};
	};

	for (kii=0;kii<n;kii++) vr[kii*n+kii] -= 1.0e0;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		daux += vr[kii]*vr[kii];
	};

	cout << " Orthogonality of q = " << sqrt(daux) << endl;

	GetRPartQrd (m, n, ar, vr);

	MvmQBlk (m, n, n,
				ar, sv,
				vr, ur);

	for (kii=0;kii<m*n;kii++) ur[kii] -= zr[kii];

	daux = 0.0e0;

	for (kii=0;kii<m*n;kii++) {
		daux += ur[kii]*ur[kii];
	};

	cout << " Diff || A - Q R || = " << sqrt(daux) << endl;

// Symmetric eigenproblem

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ar[kjj*n+kii] = 1.0e0 / ((float) (kii-kjj) + 0.5e0);
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			zr[kjj*n+kii] = ar[kjj*n+kii]+ar[kii*n+kjj];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ar[kjj*n+kii] = zr[kjj*n+kii];
		};
	};

	int info;

	dsyev_ ("V", "U", &n,
				zr, &n, sv, 
				workr, &lwork, &info);

	cout << " Eig[0] = " << sv[0] << " Eig[n-1] = " << sv[n-1] << endl;

// Check orthogonality

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			daux = 0.0e0;
			for (kkk=0;kkk<n;kkk++) {
				daux += zr[kii*n+kkk] * zr[kjj*n+kkk];
			};
			ur[kii*n+kjj] = daux;
		};
	};

	for (kii=0;kii<n;kii++) ur[kii*n+kii] -= 1.0e0;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		daux += ur[kii]*ur[kii];
	};

	cout << " Orthogonality of ur = " << sqrt(daux) << endl;

// Check the product X \Lambda X^h

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ur[kjj*n+kii] = zr[kjj*n+kii];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			vr[kjj*n+kii] = zr[kii*n+kjj];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ur[kjj*n+kii] *= sv[kjj];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			daux = 0.0e0;
			for (kkk=0;kkk<n;kkk++) {
				daux += ur[kkk*n+kii] * vr[kjj*n+kkk];
			};
			zr[kjj*n+kii] = daux;
		};
	};

	for (kii=0;kii<n*n;kii++) zr[kii] -= ar[kii];

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		daux += zr[kii]*zr[kii];
	};

	cout << " Diff || A - X Lambda X^h || = " << sqrt(daux) << endl;

// Svd

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ar[kjj*n+kii] = 1.0e0 / ((float) (kii-kjj) + 0.5e0);
			zr[kjj*n+kii] = ar[kjj*n+kii];
		};
	};

	dgesvd_ ("A", "A", &n, &n,
				zr, &n, sv, ur, &n, 
				vr, &n, workr, &lwork, &info);

	cout << " Sv[0] = " << sv[0] << " Sv[n-1] = " << sv[n-1] << endl;

// Check orthogonality

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			daux = 0.0e0;
			for (kkk=0;kkk<n;kkk++) {
				daux += ur[kii*n+kkk] * ur[kjj*n+kkk];
			};
			zr[kii*n+kjj] = daux;
		};
	};

	for (kii=0;kii<n;kii++) zr[kii*n+kii] -= 1.0e0;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		daux += zr[kii]*zr[kii];
	};

	cout << " Orthogonality of ur = " << sqrt(daux) << endl;

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			daux = 0.0e0;
			for (kkk=0;kkk<n;kkk++) {
				daux += vr[kii*n+kkk] * vr[kjj*n+kkk];
			};
			zr[kii*n+kjj] = daux;
		};
	};

	for (kii=0;kii<n;kii++) zr[kii*n+kii] -= 1.0e0;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		daux += zr[kii]*zr[kii];
	};

	cout << " Orthogonality of vr = " << sqrt(daux) << endl;

// Check the product U \Sigma V^h

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ur[kjj*n+kii] *= sv[kjj];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			daux = 0.0e0;
			for (kkk=0;kkk<n;kkk++) {
				daux += ur[kkk*n+kii] * vr[kjj*n+kkk];
			};
			zr[kjj*n+kii] = daux;
		};
	};

	for (kii=0;kii<n*n;kii++) zr[kii] -= ar[kii];

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		daux += zr[kii]*zr[kii];
	};

	cout << " Diff || A - U Sigma V^h || = " << sqrt(daux) << endl;

// Allocate complex data

	dcmplx *ac;
	dcmplx *uc;
	dcmplx *vc;
	dcmplx *zc;
	dcmplx *workc;

	ac = new dcmplx [m*n];
	if (!ac) MemoryFail (funcname);
	uc = new dcmplx [m*n];
	if (!uc) MemoryFail (funcname);
	vc = new dcmplx [n*n];
	if (!vc) MemoryFail (funcname);
	zc = new dcmplx [m*n];
	if (!zc) MemoryFail (funcname);
	workc = new dcmplx [lwork];
	if (!workc) MemoryFail (funcname);

// Test the case of complex data

// Check QR decomposition

	dcmplx caux;

	for (kii=0;kii<m;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ac[kjj*m+kii].real (1.0e0 / ((float) (kii-kjj) + 0.5e0));
			ac[kjj*m+kii].imag (-1.0e0 / ((float) (kjj-kii) + 0.43e0));
			zc[kjj*m+kii] = ac[kjj*m+kii];
		};
	};

	QrdBlk (m, n, ac, m, workc);

	for (kii=0;kii<n*n;kii++) vc[kii] = czero;
	for (kii=0;kii<n;kii++) vc[kii*n+kii] = cone;

	MvmQBlk (m, n, n,
				ac, m, workc,
				vc, n, uc, m);

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			caux = czero;
			for (kkk=0;kkk<m;kkk++) {
				caux.PlEq (uc[kii*m+kkk] * conj(uc[kjj*m+kkk]));
			};
			vc[kii*n+kjj] = caux;
		};
	};

	for (kii=0;kii<n;kii++) vc[kii*n+kii] -= cone;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		double daux1 = abs (vc[kii]);
		daux += daux1*daux1;
	};

	cout << " Orthogonality of q = " << sqrt(daux) << endl;

	GetRPartQrd (m, 0, n, ac, m, vc, n);

	MvmQBlk (m, n, n,
				ac, m, workc,
				vc, n, uc, m);

	for (kii=0;kii<m*n;kii++) uc[kii] -= zc[kii];

	for (kii=0;kii<m*n;kii++) {
		double daux1 = abs (uc[kii]);
		daux += daux1*daux1;
	};

	cout << " Diff || Ac - Q R || = " << sqrt(daux) << endl;

// Complex hermitian eigenproblem

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ac[kjj*n+kii].real (1.0e0 / ((float) (kii-kjj) + 0.5e0));
			ac[kjj*n+kii].imag (-1.0e0 / ((float) (kjj-kii) + 0.43e0));
			zc[kjj*n+kii] = ac[kjj*m+kii];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			zc[kjj*n+kii] = ac[kjj*n+kii]+conj(ac[kii*n+kjj]);
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ac[kjj*n+kii] = zc[kjj*n+kii];
		};
	};

	zheev_ ("V", "U", &n,
				zc, &n, sv, 
				workc, &lwork, workr, &info);
//	int MyZheev (char *_jobz, char *_uplo, int *_n, dcmplx *_a, 
//				int *_lda, double *_w, dcmplx *_work, int *_lwork, double *_rwork, int *_info);
//	MyZheev ("V", "U", &n,
//				zc, &n, sv, 
//				workc, &lwork, workr, &info);

	cout << " CEig[0] = " << sv[0] << " CEig[n-1] = " << sv[n-1] << endl;

// Check unitarity

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			caux = czero;
			for (kkk=0;kkk<n;kkk++) {
				caux += zc[kii*n+kkk] * conj(zc[kjj*n+kkk]);
			};
			uc[kii*n+kjj] = caux;
		};
	};

	for (kii=0;kii<n;kii++) uc[kii*n+kii] -= cone;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		double daux1 = abs (uc[kii]);
		daux += daux1*daux1;
	};

	cout << " Orthogonality of uc = " << sqrt(daux) << endl;

// Check the product X \Lambda X^h

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			uc[kjj*n+kii] = zc[kjj*n+kii];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			vc[kjj*n+kii] = conj(zc[kii*n+kjj]);
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			uc[kjj*n+kii] = uc[kjj*n+kii] * sv[kjj];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			caux = czero;
			for (kkk=0;kkk<n;kkk++) {
				caux += uc[kkk*n+kii] * vc[kjj*n+kkk];
			};
			zc[kjj*n+kii] = caux;
		};
	};

	for (kii=0;kii<n*n;kii++) zc[kii] = zc[kii] - ac[kii];

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		double daux1 = abs (zc[kii]);
		daux += daux1*daux1;
	};

	cout << " DiffC || A - X Lambda X^h || = " << sqrt(daux) << endl;

// Svd

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			ac[kjj*n+kii].real (1.0e0 / ((float) (kii-kjj) + 0.5e0));
			ac[kjj*n+kii].imag (-1.0e0 / ((float) (kjj-kii) + 0.43e0));
			zc[kjj*n+kii] = ac[kjj*n+kii];
		};
	};

	zgesvd_ ("A", "A", &n, &n,
				zc, &n, sv, uc, &n, 
				vc, &n, workc, &lwork, workr, &info);

	cout << " Sv[0] = " << sv[0] << " Sv[n-1] = " << sv[n-1] << endl;

// Check orthogonality

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			caux = czero;
			for (kkk=0;kkk<n;kkk++) {
				caux.PlEq (uc[kjj*n+kkk] * conj(uc[kii*n+kkk]));
			};
			zc[kii*n+kjj] = caux;
		};
	};

	for (kii=0;kii<n;kii++) zc[kii*n+kii] -= cone;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		double daux1 = abs (zc[kii]);
		daux += daux1*daux1;
	};

	cout << " Orthogonality of uc = " << sqrt(daux) << endl;

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			caux = czero;
			for (kkk=0;kkk<n;kkk++) {
				caux.PlEq (vc[kjj*n+kkk] * conj(vc[kii*n+kkk]));
			};
			zc[kii*n+kjj] = caux;
		};
	};

	for (kii=0;kii<n;kii++) zc[kii*n+kii] -= cone;

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		double daux1 = abs (zc[kii]);
		daux += daux1*daux1;
	};

	cout << " Orthogonality of vc = " << sqrt(daux) << endl;

// Check the product U \Sigma V^h

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			uc[kjj*n+kii] *= sv[kjj];
		};
	};

	for (kii=0;kii<n;kii++) {
		for (kjj=0;kjj<n;kjj++) {
			caux = czero;
			for (kkk=0;kkk<n;kkk++) {
				caux.PlEq (uc[kkk*n+kii] * vc[kjj*n+kkk]);
			};
			zc[kjj*n+kii] = caux;
		};
	};

	for (kii=0;kii<n*n;kii++) zc[kii] -= ac[kii];

	daux = 0.0e0;

	for (kii=0;kii<n*n;kii++) {
		double daux1 = abs (zc[kii]);
		daux += daux1*daux1;
	};

	cout << " Diff || A - U Sigma V^h || = " << sqrt(daux) << endl;

// Free work arrays

	delete [] ar;
	delete [] ur;
	delete [] vr;
	delete [] sv;
	delete [] workr;
	delete [] ac;
	delete [] uc;
	delete [] vc;
	delete [] workc;

};

// Author: Kharchenko S.A.
// Lapack function zheev implemented via double functions
//========================================================================================
int MyZheevOld (const char *_jobz, const char *_uplo, int *_n, dcmplx *_a, 
			int *_lda, double *_w, dcmplx *_work, int *_lwork, double *_rwork, int *_info) {

	const char *funcname = "MyZheevOld";

	dcmplx czero (0.0e0,0.0e0);
	double machep = 1.0e-13;

	int itype = 0;

	if (_jobz != "V" || _uplo != "U") throw " MyZheev: wrong character parameters ";

// Allocate work memory

	int nloc = *_n;
	int n2loc = 2*nloc;
	int lwork = 10 * n2loc;

	double *ar, *eig, *sv, *workr;
	dcmplx *ac, *uc, *vc, *workc;

	ar = new double [n2loc*n2loc];
	if (!ar) MemoryFail (funcname);
	eig = new double [n2loc];
	if (!eig) MemoryFail (funcname);
	sv = new double [n2loc];
	if (!sv) MemoryFail (funcname);
	workr = new double [lwork];
	if (!workr) MemoryFail (funcname);
	ac = new dcmplx [n2loc*n2loc];
	if (!ac) MemoryFail (funcname);
	uc = new dcmplx [n2loc*n2loc];
	if (!uc) MemoryFail (funcname);
	vc = new dcmplx [n2loc*n2loc];
	if (!vc) MemoryFail (funcname);
	workc = new dcmplx [lwork];
	if (!workc) MemoryFail (funcname);

// Create double matrix data

	int kii, kjj, i, j, info;

	for (kii=0;kii<nloc;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			ar[n2loc*kjj+kii] = _a[*_lda*kjj+kii].x;
			ar[n2loc*(nloc+kjj)+nloc+kii] = _a[*_lda*kjj+kii].x;
			ar[n2loc*kjj+nloc+kii] = _a[*_lda*kjj+kii].y;
			ar[n2loc*(nloc+kjj)+kii] = -_a[*_lda*kjj+kii].y;
		};
	};

// Solve symmetric eigenvalue problem

	dsyev_ (_jobz, _uplo, &n2loc, ar, &n2loc, eig, workr, &lwork, &info);

	if (info != 0) {
		throw " Error in the Lapack routine dsyev";
	};

// For each cluster of eigenvalues compute the complex invariant subspace

	double delta = fabs(eig[0]);
	if (delta<fabs(eig[n2loc-1])) delta = fabs(eig[n2loc-1]);
	delta = delta * machep;

	int ntot = 0;
	int n2tot = 0;
	int n2ini = 0;

	while (n2ini < n2loc) {
		while (true) {
			if ((n2tot)%2 == 0 || 
				(n2tot < n2loc-1 && eig[n2tot+1] < eig[n2ini]+machep)) {
				n2tot++;
			} else {
				break;
			};
		};
		for (i=n2ini;i<=n2tot;i++) {
			for (j=0;j<nloc;j++) {
				ac[(i-n2ini)*nloc+j].x = ar[n2loc*i+j];
				ac[(i-n2ini)*nloc+j].y = ar[n2loc*i+nloc+j];
			};
		};

		int nnloc = n2tot-n2ini+1;

		if (nloc >= nnloc) {
			zgesvd_ ("A", "A", &nloc, &nnloc,
						ac, &nloc, sv, 
						uc, &nloc, vc, &nloc, workc, &lwork, workr, &info);
			if (info != 0) {
				throw " Error in the Lapack routine zgesvd";
			};
			int nnnew = nnloc / 2;
			for (i=0;i<nnnew;i++) {
				for (j=0;j<nloc;j++) {
					_a[ntot*nloc+j] = uc[i*nloc+j];
				};
				_w[ntot] = eig[ntot*2];
				ntot++;
			};
		} else {
			for (i=0;i<nnloc;i++) {
				for (j=0;j<nloc;j++) {
					uc[j*nnloc+i] = conj(ac[n2loc*i+j]);
				};
			};
			zgesvd_ ("A", "A", &nnloc, &nloc,
						uc, &nnloc, sv, 
						ac, &nnloc, vc, &nnloc, workc, &lwork, workr, &info);
			if (info != 0) {
				throw " Error in the Lapack routine zgesvd";
			};
			int nnnew = nnloc / 2;
			for (i=0;i<nnnew;i++) {
				for (j=0;j<nloc;j++) {
					_a[ntot*nloc+j] = vc[j*nnloc+i];
				};
				_w[ntot] = eig[ntot*2];
				ntot++;
			};
		};
		n2ini = n2tot+1;
	};

// Free work memory

	delete [] ar;
	delete [] eig;
	delete [] sv;
	delete [] workr;
	delete [] ac;
	delete [] uc;
	delete [] vc;
	delete [] workc;

	_info = 0;

	return 0;

};

// Author: Kharchenko S.A.
// Lapack function zheev implemented via double functions
//========================================================================================
int MyZheev (const char *_jobz, const char *_uplo, int *_n, dcmplx *_a, 
			int *_lda, double *_w, dcmplx *_work, int *_lwork, double *_rwork, int *_info) {

	const char *funcname = "MyZheev";

	dcmplx czero (0.0e0,0.0e0);

	if (_jobz != "V" || _uplo != "U") throw " MyZheev: wrong character parameters ";

// Allocate work memory

	int nloc = *_n;
	int ldaloc = *_lda;
	int lwork = 10 * nloc;

	double *sv, *workr;
	dcmplx *ac, *uc, *vc, *workc;

	sv = new double [nloc];
	if (!sv) MemoryFail (funcname);
	workr = new double [lwork];
	if (!workr) MemoryFail (funcname);
	ac = new dcmplx [nloc*nloc];
	if (!ac) MemoryFail (funcname);
	uc = new dcmplx [nloc*nloc];
	if (!uc) MemoryFail (funcname);
	vc = new dcmplx [nloc*nloc];
	if (!vc) MemoryFail (funcname);
	workc = new dcmplx [lwork];
	if (!workc) MemoryFail (funcname);

// Create matrix data

	int kii, kjj, info;

	for (kii=0;kii<nloc;kii++) {
		for (kjj=0;kjj<nloc;kjj++) {
			ac[nloc*kjj+kii] = _a[ldaloc*kjj+kii];
		};
	};

// Compute singular value decomposition

	zgesvd_ ("A", "A", &nloc, &nloc,
				ac, &nloc, sv, 
				uc, &nloc, vc, &nloc, 
				workc, &lwork, workr, &info);

	if (info != 0) throw " Error in the Lapack routine ZGESVD";

// For each singular value compute its sign

	dcmplx caux;

	for (kii=0;kii<nloc;kii++) {
		caux = czero;
		for (kjj=0;kjj<nloc;kjj++) {
			caux += uc[kii*nloc+kjj] * vc[kjj*nloc+kii];
		};
		if (caux.x < 0.0e0) sv[kii] = -sv[kii];
	};

// Sort the eigenvalues

	CDoubleInt *srtarr;

	srtarr = new CDoubleInt [nloc];
	if (!srtarr) MemoryFail (funcname);

	for (kii=0;kii<nloc;kii++) {
		srtarr[kii].dvalue = sv[kii];
		srtarr[kii].intvalue = kii;
	};

	qsort (srtarr, nloc, sizeof(CDoubleInt), CDoubleInt::CompareDoubleInt);

// Return the result

	int kold;

	for (kii=0;kii<nloc;kii++) {
		kold = srtarr[kii].intvalue;
		_w[kii] = sv[kold];
		for (kjj=0;kjj<nloc;kjj++) {
			_a[kii*ldaloc+kjj] = uc[kold*nloc+kjj];
		};
	};

// Free work memory

	delete [] sv;
	delete [] workr;
	delete [] ac;
	delete [] uc;
	delete [] vc;
	delete [] workc;
	delete [] srtarr;

	_info = 0;

	return 0;

};

#ifndef _WINDOWS

// Author: Kharchenko S.A.
// Lapack function dlarfg
//========================================================================================
int dlarfg_ (int *n, double *alpha, double *x, // Lapack function dlarfg
				int *incx, double *tau) {

	dlarfgdum_ (n, alpha, x, 
				incx, tau);

	return 0;

};

// Author: Kharchenko S.A.
// Lapack function dsyev
//========================================================================================
int dsyev_ (const char *jobz, const char *uplo, int *n, double *a, // Lapack function dsyev
			int *lda, double *w, double *work, int *lwork, int *info) {

	int itype = 0;

	if (jobz == "V" && uplo == "U") {
		itype = 0;
	};

	dsyevdum_ (&itype, n, a,
				lda, w, work, lwork, info);

	return 0;

};

// Author: Kharchenko S.A.
// Lapack function dgesvd
//========================================================================================
int dgesvd_ (const char *jobu, const char *jobvt, int *m, int *n, // Lapack function dgesvd
				double *a, int *lda, double *s, double *u, int *ldu, 
				double *vt, int *ldvt, double *work, int *lwork, int *info) {

	int itype = 0;

	if (jobu == "N" && jobvt == "N") {
		itype = 0;
	} else if (jobu == "A" && jobvt == "A") {
		itype = 1;
	} else if (jobu == "O" && jobvt == "A") {
		itype = 2;
	};

	dgesvddum_ (&itype, m, n,
				a, lda, s, u, ldu, 
				vt, ldvt, work, lwork, info);

	return 0;

};

// Author: Kharchenko S.A.
// Lapack function zlarfg
//========================================================================================
int zlarfg_ (int *n, dcmplx *alpha, dcmplx *x, // Lapack function zlarfg
			int *incx, dcmplx *tau) {

	zlarfgdum_ (n, alpha, x, 
				incx, tau);

	return 0;

};

// Author: Kharchenko S.A.
// Lapack function zheev
//========================================================================================
int zheev_ (const char *jobz, const char *uplo, int *n, dcmplx *a, 
			int *lda, double *w, dcmplx *work, int *lwork, double *rwork, int *info) {

	int itype = 0;

	if (jobz == "V" && uplo == "U") {
		itype = 0;
	};

	MyZheev (jobz, uplo, n, a,
				lda, w, work, lwork, rwork, info);

//	zheevdum_ (&itype, n, a,
//				lda, w, work, lwork, rwork, info);

	return 0;

};

// Author: Kharchenko S.A.
// Lapack function zgesvd
//========================================================================================
int zgesvd_ (const char *jobu, const char *jobvt, int *m, int *n, // Lapack function zgesvd
				dcmplx *a, int *lda, double *s, dcmplx *u, int *ldu, 
				dcmplx *vt, int *ldvt, dcmplx *work, int *lwork, double *rwork, int *info) {

	int itype = 0;

	if (jobu == "N" && jobvt == "N") {
		itype = 0;
	} else if (jobu == "A" && jobvt == "A") {
		itype = 1;
	} else if (jobu == "O" && jobvt == "A") {
		itype = 2;
	};

	zgesvddum_ (&itype, m, n,
				a, lda, s, u, ldu, 
				vt, ldvt, work, lwork, rwork, info);

	return 0;

};
#endif
