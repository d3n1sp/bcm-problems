#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <iosfwd>
#include <cstring>
#include <cmath>
#include <iomanip>

#include "fourier.h"
#include "globals.h"

using namespace std;


// Author: Kharchenko S.A.
// Description: Implement fast Toeplitz 3D matrix by vector multiplication
// CFourier::FastToeplitz3DMvm()
//========================================================================================
void CFourier::FastToeplitz3DMvm (int _n1, int _n2, int _n3, int _ncols, // Implement fast Toeplitz 3D matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "FastToeplitz3DMvm";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int n12 = _n1*_n2;
	int nold = n12*_n3;
	int nloc = 4*nold;

	int n1_2 = _n1 * 2;
	int n12_2 = n12*2;
	int n12_4 = n12*4;
	int ldwork = nloc*_ncols;
	int ld1 = _n1*2-1;
	int ld2 = _n2*2-1;
	int ld3 = _n3*2-1;
	int ld12 = ld1*ld2;
	int ldcoefs = nloc*2-1;

	dcmplx *xwork, *fxwork;
	dcmplx *fcoefs;

	xwork = new dcmplx [ldwork];
	if (!xwork) MemoryFail (funcname);
	fxwork = new dcmplx [ldwork];
	if (!fxwork) MemoryFail (funcname);
	fcoefs = new dcmplx [ldcoefs];
	if (!fcoefs) MemoryFail (funcname);

// Fill local multiplication arrays

	int icol, i, j, k, ishftold, ishft, indold, ind;

	for (i=0;i<ldwork;i++) xwork[i] = czero;
	for (i=0;i<ldwork;i++) fxwork[i] = czero;
	for (i=0;i<ldcoefs;i++) fcoefs[i] = czero;

	for (icol=0;icol<_ncols;icol++) {
		ishftold = icol*nold;
		ishft = icol*nloc;
		for (i=0;i<_n1;i++) {
			for (j=0;j<_n2;j++) {
				for (k=0;k<_n3;k++) {
					indold = k*n12+j*_n1+i;
					ind = k*n12_4+j*n1_2+i;
					xwork[ishft+ind] = _x[ishftold+indold];
				};
			};
		};
	};

	for (k=0;k<ld3;k++) {
		for (j=0;j<ld2;j++) {
			for (i=0;i<ld1;i++) {
				indold = k*ld12+j*ld1+i;
				ind = n12_2+k*n12_4+_n1+j*n1_2+i;
				fcoefs[ind] = _coefs[indold];
			};
		};
	};

// Perform fast Toeplitz multiplication

	FastToeplitzMvm (nloc, _ncols,
							fcoefs,
							xwork, fxwork);

// Store data back in initial format

	for (icol=0;icol<_ncols;icol++) {
		ishftold = icol*nold;
		ishft = icol*nloc;
		for (i=0;i<_n1;i++) {
			for (j=0;j<_n2;j++) {
				for (k=0;k<_n3;k++) {
					indold = k*n12+j*_n1+i;
					ind = k*n12_4+j*n1_2+i;
					_fx[ishftold+indold] = fxwork[ishft+ind];
				};
			};
		};
	};

// Free work arrays

	delete [] xwork;
	delete [] fxwork;
	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement fast Toeplitz 3D hermitian matrix by vector multiplication
// CFourier::FastToeplitz3DMvmH()
//========================================================================================
void CFourier::FastToeplitz3DMvmH (int _n1, int _n2, int _n3, int _ncols, // Implement fast Toeplitz 3D hermitian matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "FastToeplitz3DMvmH";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ld1 = _n1*2-1;
	int ld2 = _n2*2-1;
	int ld3 = _n3*2-1;
	int ldcoefs = ld1*ld2*ld3;
	int ld12 = ld1*ld2;
	int ld23 = ld2*ld3;

	dcmplx *fcoefs;

	fcoefs = new dcmplx [ldcoefs];
	if (!fcoefs) MemoryFail (funcname);

// Fill local multiplication arrays

	int i, j, k, ind, indold;

	for (i=0;i<ld1;i++) {
		for (j=0;j<ld2;j++) {
			for (k=0;k<ld3;k++) {
				indold = k*ld12+j*ld1+i;
				ind = (ld3-k-1)*ld12+(ld2-j-1)*ld1+ld1-i-1;
				fcoefs[ind] = conj(_coefs[indold]);
			};
		};
	};

// Perform fast Toeplitz multiplication

	FastToeplitz3DMvm (_n1, _n2, _n3, _ncols,
								fcoefs,
								_x, _fx);

// Free work arrays

	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement general Toeplitz 3D matrix by vector multiplication
// CFourier::ExplicitToeplitz3DMvm()
//========================================================================================
void CFourier::ExplicitToeplitz3DMvm (int _n1, int _n2, int _n3, int _ncols, // Implement general Toeplitz 3D matrix by vector multiplication
													dcmplx *_coefs,
													dcmplx *_x,
													dcmplx *_fx) {

	const char *funcname = "ExplicitToeplitz3DMvm";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int icol, i1, j1, k1, i2, j2, k2, ishft, ind, ind1, ind2;

	dcmplx caux;

	int nloc = _n1*_n2*_n3;
	int n12 = _n1*_n2;
	int ld1 = _n1*2-1;
	int ld2 = _n2*2-1;
	int ld12 = ld1*ld2;
	int ishft1, ishft2, ishft3;

	for (icol=0;icol<_ncols;icol++) {
		ishft = icol*nloc;
		for (i1=0;i1<_n1;i1++) {
			for (j1=0;j1<_n2;j1++) {
				for (k1=0;k1<_n3;k1++) {
					ind1 = k1*n12+j1*_n1+i1;
					caux = czero;
					for (i2=0;i2<_n1;i2++) {
						for (j2=0;j2<_n2;j2++) {
							for (k2=0;k2<_n3;k2++) {
								ind2 = k2*n12+j2*_n1+i2;
								ishft1 = i2-i1+_n1-1;
								ishft2 = j2-j1+_n2-1;
								ishft3 = k2-k1+_n3-1;
								ind = ishft3*ld12+ishft2*ld1+ishft1;
								caux += _coefs[ind] * _x[ishft+ind2];
							};
						};
					};
					_fx[ishft+ind1] = caux;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Expand the Toeplitz 3D matrix
// CFourier::ExpandToeplitz3D()
//========================================================================================
void CFourier::ExpandToeplitz3D (int _n1, int _n2, int _n3, // Expand the Toeplitz 3D matrix
											dcmplx *_coefs,
											dcmplx *_amatr) {

	const char *funcname = "ExpandToeplitz3D";

// Multiply

	int nloc = _n1*_n2*_n3;
	int n12 = _n1*_n2;
	int ld1 = _n1*2-1;
	int ld2 = _n2*2-1;
	int ld12 = ld1*ld2;

	int ishft1, ishft2, ishft3;

	int i1, j1, k1, i2, j2, k2, ind, ind1, ind2, indnew;

	for (i1=0;i1<_n1;i1++) {
		for (j1=0;j1<_n2;j1++) {
			for (k1=0;k1<_n3;k1++) {
				ind1 = k1*n12+j1*_n1+i1;
				for (i2=0;i2<_n1;i2++) {
					for (j2=0;j2<_n2;j2++) {
						for (k2=0;k2<_n3;k2++) {
							ind2 = k2*n12+j2*_n1+i2;
							ishft1 = i2-i1+_n1-1;
							ishft2 = j2-j1+_n2-1;
							ishft3 = k2-k1+_n3-1;
							ind = ishft3*ld12+ishft2*ld1+ishft1;
							indnew = nloc*ind2+ind1;
							_amatr[indnew] = _coefs[ind];
						};
					};
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Implement fast Toeplitz 2D matrix by vector multiplication
// CFourier::FastToeplitz2DMvm()
//========================================================================================
void CFourier::FastToeplitz2DMvm (int _n1, int _n2, int _ncols, // Implement fast Toeplitz 2D matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "FastToeplitz2DMvm";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int nold = _n1*_n2;
	int nloc = 2*nold;

	int n1_2 = _n1 * 2;
	int ldwork = nloc*_ncols;
	int ld1 = _n1*2-1;
	int ld2 = _n2*2-1;
	int ld12 = ld1*ld2;
	int ldcoefs = nloc*2-1;

	dcmplx *xwork, *fxwork;
	dcmplx *fcoefs;

	xwork = new dcmplx [ldwork];
	if (!xwork) MemoryFail (funcname);
	fxwork = new dcmplx [ldwork];
	if (!fxwork) MemoryFail (funcname);
	fcoefs = new dcmplx [ldcoefs];
	if (!fcoefs) MemoryFail (funcname);

// Fill local multiplication arrays

	int icol, i, j, ishftold, ishft, indold, ind;

	for (i=0;i<ldwork;i++) xwork[i] = czero;
	for (i=0;i<ldwork;i++) fxwork[i] = czero;
	for (i=0;i<ldcoefs;i++) fcoefs[i] = czero;

	for (icol=0;icol<_ncols;icol++) {
		ishftold = icol*nold;
		ishft = icol*nloc;
		for (i=0;i<_n1;i++) {
			for (j=0;j<_n2;j++) {
				indold = j*_n1+i;
				ind = j*n1_2+i;
				xwork[ishft+ind] = _x[ishftold+indold];
			};
		};
	};

	for (j=0;j<ld2;j++) {
		for (i=0;i<ld1;i++) {
			indold = j*ld1+i;
			ind = _n1+j*n1_2+i;
			fcoefs[ind] = _coefs[indold];
		};
	};

// Perform fast Toeplitz multiplication

	FastToeplitzMvm (nloc, _ncols,
							fcoefs,
							xwork, fxwork);

// Store data back in initial format

	for (icol=0;icol<_ncols;icol++) {
		ishftold = icol*nold;
		ishft = icol*nloc;
		for (i=0;i<_n1;i++) {
			for (j=0;j<_n2;j++) {
				indold = j*_n1+i;
				ind = j*n1_2+i;
				_fx[ishftold+indold] = fxwork[ishft+ind];
			};
		};
	};

// Free work arrays

	delete [] xwork;
	delete [] fxwork;
	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement fast Toeplitz 2D matrix by vector multiplication
// CFourier::FastToeplitz2DMvmH()
//========================================================================================
void CFourier::FastToeplitz2DMvmH (int _n1, int _n2, int _ncols, // Implement fast Toeplitz 2D matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "FastToeplitz2DMvmH";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ld1 = _n1*2-1;
	int ld2 = _n2*2-1;
	int ldcoefs = ld1*ld2;

	dcmplx *fcoefs;

	fcoefs = new dcmplx [ldcoefs];
	if (!fcoefs) MemoryFail (funcname);

// Fill local multiplication arrays

	int i, j, ind, indold;

	for (i=0;i<ld1;i++) {
		for (j=0;j<ld2;j++) {
			indold = j*ld1+i;
			ind = (ld2-j-1)*ld1+ld1-i-1;
			fcoefs[ind] = conj(_coefs[indold]);
		};
	};

// Perform fast Toeplitz multiplication

	FastToeplitz2DMvm (_n1, _n2, _ncols,
								fcoefs,
								_x, _fx);

// Free work arrays

	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement general Toeplitz 2D matrix by vector multiplication
// CFourier::ExplicitToeplitz2DMvm()
//========================================================================================
void CFourier::ExplicitToeplitz2DMvm (int _n1, int _n2, int _ncols, // Implement general Toeplitz 2D matrix by vector multiplication
													dcmplx *_coefs,
													dcmplx *_x,
													dcmplx *_fx) {

	const char *funcname = "ExplicitToeplitz2DMvm";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int icol, i1, j1, i2, j2, ishft, ind, ind1, ind2;

	dcmplx caux;

	int nloc = _n1*_n2;
	int ld1 = _n1*2-1;
	int ishft1, ishft2;

	for (icol=0;icol<_ncols;icol++) {
		ishft = icol*nloc;
		for (i1=0;i1<_n1;i1++) {
			for (j1=0;j1<_n2;j1++) {
				ind1 = j1*_n1+i1;
				caux = czero;
				for (i2=0;i2<_n1;i2++) {
					for (j2=0;j2<_n2;j2++) {
						ind2 = j2*_n1+i2;
						ishft1 = i2-i1+_n1-1;
						ishft2 = j2-j1+_n2-1;
						ind = ishft2*ld1+ishft1;
						caux += _coefs[ind] * _x[ishft+ind2];
					};
				};
				_fx[ishft+ind1] = caux;
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Expand the Toeplitz 2D matrix
// CFourier::ExpandToeplitz2D()
//========================================================================================
void CFourier::ExpandToeplitz2D (int _n1, int _n2, // Expand the Toeplitz 2D matrix
											dcmplx *_coefs,
											dcmplx *_amatr) {

	const char *funcname = "ExpandToeplitz2D";

// Multiply

	int nloc = _n1*_n2;
	int ld1 = _n1*2-1;
	int ishft1, ishft2;

	int i1, j1, i2, j2, ind, ind1, ind2, indnew;

	for (i1=0;i1<_n1;i1++) {
		for (j1=0;j1<_n2;j1++) {
			ind1 = j1*_n1+i1;
			for (i2=0;i2<_n1;i2++) {
				for (j2=0;j2<_n2;j2++) {
					ind2 = j2*_n1+i2;
					ishft1 = i2-i1+_n1-1;
					ishft2 = j2-j1+_n2-1;
					ind = ishft2*ld1+ishft1;
					indnew = nloc*ind2+ind1;
					_amatr[indnew] = _coefs[ind];
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Implement fast Toeplitz matrix by vector multiplication
// CFourier::FastToeplitzMvm()
//========================================================================================
void CFourier::FastToeplitzMvm (int _n, int _ncols, // Implement fast Toeplitz matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "FastToeplitzMvm";

	dcmplx czero (0.0e0,0.0e0);

// Find appropriate degree (power of two) for circulant matrix

	int nloc = 1;

	while (nloc < _n) {
		nloc *= 2;
	};

	nloc *= 2;

// Allocate work memory

	int ldwork = nloc*_ncols;

	dcmplx *xwork, *fxwork;
	dcmplx *fcoefs;

	xwork = new dcmplx [ldwork];
	if (!xwork) MemoryFail (funcname);
	fxwork = new dcmplx [ldwork];
	if (!fxwork) MemoryFail (funcname);
	fcoefs = new dcmplx [nloc];
	if (!fcoefs) MemoryFail (funcname);

// Fill local multiplication arrays

	int i, j;

	for (i=0;i<ldwork;i++) xwork[i] = czero;
	for (i=0;i<ldwork;i++) fxwork[i] = czero;

	for (i=0;i<_ncols;i++) {
		for (j=0;j<_n;j++) {
			xwork[i*nloc+j] = _x[i*_n+j];
		};
	};

	for (i=0;i<nloc;i++) fcoefs[i] = czero;

	for (i=0;i<_n;i++) fcoefs[i] = _coefs[_n-1-i];
	for (i=0;i<_n-1;i++) fcoefs[nloc-i-1] = _coefs[_n+i];

// Perform fast circulant multiplication

	FastCirculantMvm (nloc, _ncols,
							fcoefs,
							xwork, fxwork);

// Store data back in initial format

	for (i=0;i<_ncols;i++) {
		for (j=0;j<_n;j++) {
			_fx[i*_n+j] = fxwork[i*nloc+j];
		};
	};

// Free work arrays

	delete [] xwork;
	delete [] fxwork;
	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement fast Toeplitz hermitian matrix by vector multiplication
// CFourier::FastToeplitzMvm()
//========================================================================================
void CFourier::FastToeplitzMvmH (int _n, int _ncols, // Implement fast Toeplitz hermitian matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "FastToeplitzMvmH";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work memory

	int ldcoef = 2*_n-1;

	dcmplx *fcoefs;

	fcoefs = new dcmplx [ldcoef];
	if (!fcoefs) MemoryFail (funcname);

// Fill local multiplication arrays

	int i;

	for (i=0;i<ldcoef;i++) {
		fcoefs[i] = conj(_coefs[ldcoef-i-1]);
	};

// Perform fast circulant multiplication

	FastToeplitzMvm (_n, _ncols,
							fcoefs,
							_x, _fx);

// Free work arrays

	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement general Toeplitz matrix by vector multiplication
// CFourier::ExplicitToeplitzMvm()
//========================================================================================
void CFourier::ExplicitToeplitzMvm (int _n, int _ncols, // Implement general Toeplitz matrix by vector multiplication
												dcmplx *_coefs,
												dcmplx *_x,
												dcmplx *_fx) {

	const char *funcname = "ExplicitToeplitzMvm";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int i, j, k, ishft, ind;

	dcmplx caux;

	for (i=0;i<_ncols;i++) {
		ishft = i*_n;
		for (j=0;j<_n;j++) {
			caux = czero;
			for (k=0;k<_n;k++) {
				ind = k-j+_n-1;
				caux += _coefs[ind] * _x[ishft+k];
			};
			_fx[ishft+j] = caux;
		};
	};

};

// Author: Kharchenko S.A.
// Description: Expand the Toeplitz matrix
// CFourier::ExpandToeplitz()
//========================================================================================
void CFourier::ExpandToeplitz (int _n, // Expand the Toeplitz matrix
											dcmplx *_coefs,
											dcmplx *_amatr) {

	const char *funcname = "ExpandToeplitz";

// Multiply

	int i, j, ind, ind1;

	for (i=0;i<_n;i++) {
		for (j=0;j<_n;j++) {
			ind = j-i+_n-1;
			ind1 = _n*j+i;
			_amatr[ind1] = _coefs[ind];
		};
	};

};

// Author: Kharchenko S.A.
// Description: Implement fast circulant matrix by vector multiplication
// CFourier::FastCirculantMvm()
//========================================================================================
void CFourier::FastCirculantMvm (int _n, int _ncols, // Implement fast circulant matrix by vector multiplication
													dcmplx *_coefs,
													dcmplx *_x,
													dcmplx *_fx) {

	const char *funcname = "FastCirculantMvm";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work arrays

	int ldwork = _n*_ncols;

	dcmplx *fxwork;
	dcmplx *fcoefs;

	fxwork = new dcmplx [ldwork];
	if (!fxwork) MemoryFail (funcname);
	fcoefs = new dcmplx [_n];
	if (!fcoefs) MemoryFail (funcname);

// Multiply

	FFTPower2 (_n, 1,
					_coefs, fcoefs);

	FFTPower2 (_n, _ncols,
					_x, fxwork);

	int icol, kii, indi;
	dcmplx caux;

	for (icol=0;icol<_ncols;icol++) {
		for (kii=0;kii<_n;kii++) {
			indi = icol*_n+kii;
			caux = fxwork[indi] * fcoefs[kii];
			fxwork[indi] = conj(caux);
		};
	};

	FFTPower2 (_n, _ncols,
					fxwork, _fx);

	double dmone = -1.0e0;
	double dn_inv = 1.0e0 / ((double) _n);

	for (kii=0;kii<_n*_ncols;kii++) {
		_fx[kii].y *= dmone;
		caux = _fx[kii] * dn_inv;
		_fx[kii] = caux;
	};

// Free work arrays

	delete [] fxwork;
	delete [] fcoefs;

};

// Author: Kharchenko S.A.
// Description: Implement Fast Fourier transform for the power of 2
// CFourier::FFTPower2()
//========================================================================================
void CFourier::FFTPower2 (int _n, int _ncols, // Implement Fast Fourier transform for the power of 2
							dcmplx *_x,
							dcmplx *_fx) {

	const char *funcname = "FFTPower2";

//	CFourier::FFTPower2Impl (_n, _ncols,
//										_x, _fx);

	int i, j, ind;
	double *pfxre, *pfxim;

	pfxre = new double [_n];
	if (!pfxre) MemoryFail (funcname);
	pfxim = new double [_n];
	if (!pfxim) MemoryFail (funcname);

	for (i=0;i<_ncols;i++) {

		for (j=0;j<_n;j++) {
			ind = j+i*_n;
			pfxre[j] = _x[ind].x;
			pfxim[j] = _x[ind].y;
		};

		CFourier::FTF1C (pfxre, pfxim, 
								_n, 1, 1, -1.0e0);

		for (j=0;j<_n;j++) {
			ind = j+i*_n;
			_fx[ind].x = pfxre[j];
			_fx[ind].y = pfxim[j];
		};

	};

	delete [] pfxre;
	delete [] pfxim;

};

// Author: Kharchenko S.A.
// Description: Implement Fast Fourier transform for the power of 2
// CFourier::FFTPower2Impl()
//========================================================================================
void CFourier::FFTPower2Impl (int _n, int _ncols, // Implement Fast Fourier transform for the power of 2
							dcmplx *_x,
							dcmplx *_fx) {

	const char *funcname = "FFTPower2Impl";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

// Check degree on entry

	int degree = 0;
	int nloc = 1;

	while (nloc < _n) {
		nloc *= 2;
		degree++;
	};

	if (nloc != _n) throw " CFourier::FFTPower2: array size is not a power of 2";

// Set local arrays that control multilevel computations

	const int degreemax = 26;

	static const int nvalues [degreemax+1] = {1,
																2,    4,    8,   16,   32,   64,   128,   256,   512,   1024,
															2048, 4096, 8192,16384,32768,65536,131072,262144,524288,1048576,
															2097152, 4194304, 8388608, 16777216, 33554432, 67108864};

	static int ilevstate  [degreemax+1];
	static int ilevsize   [degreemax+1];
	static int ilevsize2  [degreemax+1];
	static int ilevdegree [degreemax+1];
	static int iacoefs    [degreemax+2];

	if (degree > degreemax) throw " CFourier::FFTPower2: array size is too large";

	ilevsize   [0] = _n / 2;
	ilevdegree [0] = degree;

	int i;

	for (i=1;i<degree;i++) {
		ilevsize   [i] = ilevsize   [i-1] / 2;
		ilevdegree [i] = ilevdegree [i-1] - 1;
	};

	for (i=0;i<degree;i++) {
		ilevsize2[i] = _n / (2 * ilevsize[i]);
	};

// Allocate work arrays that store the temporary data

	int n_2 = _n / 2;
	int nzcoefs = 4 * _n;
	int ldlev = _n*_ncols;

	dcmplx *fxlev;
	dcmplx *fxwork;
	dcmplx *coefs;

	fxlev = new dcmplx [ldlev];
	if (!fxlev) MemoryFail (funcname);
	fxwork = new dcmplx [ldlev];
	if (!fxwork) MemoryFail (funcname);
	coefs = new dcmplx [nzcoefs];
	if (!coefs) MemoryFail (funcname);

// Init coefficients

	int idegree;

	int nz = 0;
	iacoefs[0] = 0;

	for (idegree=0;idegree<=degree;idegree++) {

		nloc = nvalues[idegree];

		if (nz+nloc > nzcoefs) throw " CFourier::FFTPower2: insufficient size of coefs array";

		double phi = 2 * DPI / ((double) nloc);

		double dre = cos(phi);
		double dim = sin(phi);

		dcmplx ceps (dre,dim);

		coefs[nz] = cone;

		for (i=1;i<nloc;i++) coefs[nz+i] = coefs[nz+i-1] * ceps;

		nz += nloc;

		iacoefs[idegree+1] = nz;

	};

// Init data of the main cycle

	for (i=0;i<ldlev;i++) fxlev[i] = _x[i];

// Main cycle over the levels

	enum FFT_STATE {INI_STATE = 0,
					COLS_FFT_STATE,
					OK_STATE};

	int ilev = 0;

	ilevstate[ilev] = INI_STATE;

	while (ilevstate[0] != OK_STATE) {

// Initial FFT by rows (n1=2)

		if (ilevstate[ilev] == INI_STATE) {

			ExplicitFourier ('R', 2, ilevsize[ilev], _ncols*ilevsize2[ilev],
									fxlev, fxwork);

// Perform the scaling

			idegree = ilevdegree[ilev];

			ScalingFourier (2, ilevsize[ilev], _ncols*ilevsize2[ilev],
									coefs+iacoefs[idegree],
									fxwork);

// Final FFT by columns (recursive call if necessary)

			ilevstate[ilev] = COLS_FFT_STATE;
			if (ilevsize[ilev] == 2) {

				ExplicitFourier ('C', 2, ilevsize[ilev], _ncols*ilevsize2[ilev],
										fxwork,
										fxlev);

			} else {
				for (i=0;i<ldlev;i++) fxlev[i] = fxwork[i];
				ilevstate[ilev+1] = INI_STATE;
				ilev++;
			};
		};

// Transpose the data

		if (ilevstate[ilev] == COLS_FFT_STATE) {

			ilevstate[ilev] = OK_STATE;

			TransposeFourier (2, ilevsize[ilev], _ncols*ilevsize2[ilev],
									fxlev,
									fxwork);

			if (ilev == 0) {
				for (i=0;i<ldlev;i++) _fx[i] = fxwork[i];
			} else {
				for (i=0;i<ldlev;i++) fxlev[i] = fxwork[i];
			};
			ilev--;
		};

// End point of the cycle

	};

// Free work arrays

	delete [] fxlev;
	delete [] fxwork;
	delete [] coefs;

};

// Author: Kharchenko S.A.
// Description: Implement implicit two level Fourier transform by columns
// CFourier::ImplicitTwoLevelFourier()
//========================================================================================
void CFourier::ImplicitTwoLevelFourier (int _n, int _ncols, int _ndivide, // Implement implicit two level Fourier transform by columns
										dcmplx *_x,
										dcmplx *_fx) {

	const char *funcname = "ImplicitTwoLevelFourier";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

// Compute Fourier parameter

	int ndivide2 = _n / _ndivide;

	double phi = 2 * DPI / ((double) _n);

	double dre = cos(phi);
	double dim = sin(phi);

	dcmplx ceps (dre,dim);

// Compute set of degrees

	dcmplx *cepsarr;

	cepsarr = new dcmplx [_n];
	if (!cepsarr) MemoryFail (funcname);

	cepsarr[0] = cone;

	int i;

	for (i=1;i<_n;i++) cepsarr[i] = cepsarr[i-1] * ceps;

// Perform computations

	dcmplx *fxloc;

	fxloc = new dcmplx [_n*_ncols];
	if (!fxloc) MemoryFail (funcname);

// Rows Fourier

	ExplicitFourier ('R', _ndivide, ndivide2, _ncols,
							_x, _fx);

// Scaling by rotation factors

	ScalingFourier (_ndivide, ndivide2, _ncols,
							cepsarr, 
							_fx);

// Columns Fourier

	ExplicitFourier ('C', ndivide2, _ndivide, _ncols,
							_fx, fxloc);

// Transpose

	TransposeFourier (_ndivide, ndivide2, _ncols,
								fxloc, 
								_fx);

// Free work array

	delete [] cepsarr;
	delete [] fxloc;

};

// Author: Kharchenko S.A.
// Description: Perform scaling of the matrix
// CFourier::ScalingFourier()
//========================================================================================
void CFourier::ScalingFourier (int _n1, int _n2, int _ncols, // Perform scaling of the matrix
						dcmplx *_cepsarr, 
						dcmplx *_x) {

//	const char *funcname = "ScalingFourier";

	int nloc = _n1 * _n2;

	int kii, kjj, kkk, i, imod, indi;
	dcmplx caux;

	for (kii=0;kii<_n2;kii++) {
		for (kjj=0;kjj<_n1;kjj++) {
			kkk = kjj*_n2+kii;
			imod = (kii*kjj) % nloc;
			for (i=0;i<_ncols;i++) {
				indi = i*nloc+kkk;
				caux = _x[indi] * _cepsarr[imod];
				_x[indi] = caux;
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Perform transposition of the matrix
// CFourier::TransposeFourier()
//========================================================================================
void CFourier::TransposeFourier (int _n1, int _n2, int _ncols, // Perform transposition of the matrix
						dcmplx *_x, 
						dcmplx *_fx) {

//	const char *funcname = "TransposeFourier";

	int nloc = _n1 * _n2;

	int kii, kjj, kkk, kkk1, i;
	dcmplx caux;

	for (kii=0;kii<_n2;kii++) {
		for (kjj=0;kjj<_n1;kjj++) {
			kkk = kjj*_n2+kii;
			kkk1 = kii*_n1+kjj;
			for (i=0;i<_ncols;i++) {
				_fx[i*nloc+kkk1] = _x[i*nloc+kkk];
			};
		};
	};

};

// Author: Kharchenko S.A.
// Description: Implement explicit Fourier transform by columns/rows
// CFourier::ExplicitFourier()
//========================================================================================
void CFourier::ExplicitFourier (char _descr, int _n, int _n2, int _ncols, // Implement explicit Fourier transform by columns/rows
								dcmplx *_x,
								dcmplx *_fx) {

	const char *funcname = "ExplicitFourier";

	dcmplx czero (0.0e0,0.0e0);
	dcmplx cone (1.0e0,0.0e0);

// Compute Fourier parameter

	double phi = 2 * DPI / ((double) _n);

	double dre = cos(phi);
	double dim = sin(phi);

	dcmplx ceps (dre,dim);

// Compute set of degrees

	dcmplx *cepsarr;

	cepsarr = new dcmplx [_n];
	if (!cepsarr) MemoryFail (funcname);

	cepsarr[0] = cone;

	int i;

	for (i=1;i<_n;i++) cepsarr[i] = cepsarr[i-1] * ceps;

// Perform computations

	int kii, kjj, imod, indi, indj, icol;

	for (i=0;i<_n*_n2*_ncols;i++) _fx[i] = czero;

	if (_descr == 'C') {
		for (kii=0;kii<_n;kii++) {
			for (kjj=0;kjj<_n;kjj++) {
				imod = (kii*kjj) % _n;
				for (i=0;i<_n2*_ncols;i++) {
					indi = i*_n+kii;
					indj = i*_n+kjj;
					_fx[indi] += cepsarr[imod] * _x[indj];
				};
			};
		};
	} else {
		int ldmatr = _n*_n2;
		for (kii=0;kii<_n;kii++) {
			for (kjj=0;kjj<_n;kjj++) {
				imod = (kii*kjj) % _n;
				for (icol=0;icol<_ncols;icol++) {
					for (i=0;i<_n2;i++) {
						indi = icol*ldmatr+kii*_n2+i;
						indj = icol*ldmatr+kjj*_n2+i;
						_fx[indi] += cepsarr[imod] * _x[indj];
					};
				};
			};
		};
	};

// Free work array

	delete [] cepsarr;

};

// Author: Kharchenko S.A.
// Description: Implement explicit circulant matrix by vector multiplication
// CFourier::ExplicitCirculantMvm()
//========================================================================================
void CFourier::ExplicitCirculantMvm (int _n, int _ncols, // Implement explicit circulant matrix by vector multiplication
													dcmplx *_coefs,
													dcmplx *_x,
													dcmplx *_fx) {

	const char *funcname = "ExplicitCirculantMvm";

	dcmplx czero (0.0e0,0.0e0);

// Main cycle

	int kii, kjj, icol;
	int ind, indi, indj;

	for (kii=0;kii<_n*_ncols;kii++) _fx[kii] = czero;

	for (kii=0;kii<_n;kii++) {
		for (kjj=0;kjj<_n;kjj++) {
			ind = kii-kjj+_n;
			if (ind >= _n) ind -= _n;
			for (icol=0;icol<_ncols;icol++) {
				indi = icol*_n+kii;
				indj = icol*_n+kjj;
				_fx[indi] += _coefs[ind] * _x[indj];
			};
		};
	};

};

// Author: Unknown
// Description: Implement Fast Fourier Transform for the power of 2
// CFourier::FTF1C()
//========================================================================================
void CFourier::FTF1C (double *_ar0, double *_ai0, // Implement Fast Fourier transform for the power of 2
								int _n, int _in, int _k, double _p)
{

// _ar0 - real part
// _ai0 - image part
// _n   - dimension of arrays ar0, ai0

	double *ar, *ai;
	int i, j, ii, jj, m, mm, n1, n2;
	double r, t, c, s, co, si, a, b;

	ar=_ar0-1;
	ai=_ai0-1;

	n1=_n/2; mm=n1/2;
	n2=n1+_k;
	j=_in; 
	jj=j;

	do {
		j+=_k;
		if (j>n1) break;
		ii=jj+n1;
		r=ar[j];  ar[j]=ar[ii];  ar[ii]=r;
		r=ai[j];  ai[j]=ai[ii];  ai[ii]=r;
		j+=_k;
		m=mm;
		while (jj>m) {
			jj-=m;
			m=m/2;
		}
		jj+=m;
		if (jj<=j) continue;
		r=ar[j];  ar[j]=ar[jj];  ar[jj]=r;
		r=ai[j];  ai[j]=ai[jj];  ai[jj]=r;
		i=j+n2;
		ii=jj+n2;
		r=ar[i];  ar[i]=ar[ii];  ar[ii]=r;
		r=ai[i];  ai[i]=ai[ii];  ai[ii]=r;
	} while (1);

	i=_k;
	t=DPI;
	if (_p==0) return;
	if (_p>0) t=-t;
	_p=-t;

	do {
		si=0;      co=1.;
		s=sin(t);  c=cos(t);
		t=0.5*t;
		ii=i;
		i*=2;
		for (m=_in; m<=ii; m+=_k) {
			for (j=m; j<=_n; j+=i) {
				jj=j+ii;
				a=ar[jj];
				b=ai[jj];
				r=a*co-b*si;
				ar[jj]=ar[j]-r;
				ar[j] =ar[j]+r;
				r=b*co+a*si;
				ai[jj]=ai[j]-r;
				ai[j] =ai[j]+r;
			}
			r= c*co-s*si;
			si=c*si+s*co;
			co=r;
		}
	} while (i<_n);
	return;
}

// Author: Kharchenko S.A.
// Description: Implement fast diagonal extraction of the Toeplitz products
// CFourier::ExplicitDiagonalPartT2TMultiply()
//========================================================================================
void CFourier::FastDiagonalPartT2TMultiply (int _n, int _idiag, // Implement fast diagonal extraction of the Toeplitz products
												dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB) {

	const char *funcname = "FastDiagonalPartT2TMultiply";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int i, i1, j, k, inda, indb, ibeg, iend;
	dcmplx caux;

	if (_idiag >= 0) {
		ibeg = 0;
		iend = _n-_idiag;
	} else {
		ibeg = -_idiag;
		iend = _n;
	};

	i = ibeg;
	j = i+_idiag;
	caux = czero;
	for (k=0;k<_n;k++) {
		inda = k-i+_n-1;
		indb = j-k+_n-1;
		caux += _coefsA[inda]*_coefsB[indb];
	};
	_diagAB[i] = caux;
	for (i=ibeg+1;i<iend;i++) {
		caux = _diagAB[i-1];
		i1 = i-1;
		inda = -i1-1+_n-1;
		indb = _idiag+i1+_n;
		caux += _coefsA[inda]*_coefsB[indb];
		inda = _n-i1-1+_n-1;
		indb = _idiag+i1;
		caux -= _coefsA[inda]*_coefsB[indb];
		_diagAB[i] = caux;
	};

};

// Author: Kharchenko S.A.
// Description: Implement explicit diagonal extraction of the Toeplitz products
// CFourier::ExplicitDiagonalPartT2TMultiply()
//========================================================================================
void CFourier::ExplicitDiagonalPartT2TMultiply (int _n, int _idiag, // Implement explicit diagonal extraction of the Toeplitz products
												dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB) {

	const char *funcname = "ExplicitDiagonalPartT2TMultiply";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int i, j, k, inda, indb, ibeg, iend;
	dcmplx caux;

	if (_idiag >= 0) {
		ibeg = 0;
		iend = _n-_idiag;
	} else {
		ibeg = -_idiag;
		iend = _n;
	};

	for (i=ibeg;i<iend;i++) {
		j = i+_idiag;
		caux = czero;
		for (k=0;k<_n;k++) {
			inda = k-i+_n-1;
			indb = j-k+_n-1;
			caux += _coefsA[inda]*_coefsB[indb];
		};
		_diagAB[i] = caux;
	};

};

// Author: Kharchenko S.A.
// Description: Implement fast diagonal extraction of the 2D Toeplitz products
// CFourier::FastDiagonalPartT2T2DMultiply()
//========================================================================================
void CFourier::FastDiagonalPartT2T2DMultiply (int _n1, int _n2, int _idiag1, int _idiag2, // Implement fast diagonal extraction of the 2D Toeplitz products
												dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB, dcmplx *_work) {

	const char *funcname = "FastDiagonalPartT2T2DMultiply";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int ld1 = 2*_n1-1;
	dcmplx *pwork, *pwork1;
	int i2, i2prev, j2, i, k, inda, indb, ibeg2, iend2;
	dcmplx caux;

	pwork = _work;
	pwork1 = _work+_n1;

	if (_idiag2 >= 0) {
		ibeg2 = 0;
		iend2 = _n2-_idiag2;
	} else {
		ibeg2 = -_idiag2;
		iend2 = _n2;
	};

	i2 = ibeg2;
	j2 = i2+_idiag2;
	for (i=0;i<_n1;i++) pwork[i] = czero;
	for (k=0;k<_n2;k++) {
		inda = k-i2+_n2-1;
		indb = j2-k+_n2-1;
		for (i=0;i<_n1;i++) pwork1[i] = czero;
		CFourier::FastDiagonalPartT2TMultiply (_n1, _idiag1,
															_coefsA+inda*ld1, _coefsB+indb*ld1, pwork1);
		for (i=0;i<_n1;i++) pwork[i] += pwork1[i];
	};
	for (i=0;i<_n1;i++) _diagAB[ibeg2*_n1+i] = pwork[i];
	for (i2=ibeg2+1;i2<iend2;i2++) {
		i2prev = i2-1;
		for (i=0;i<_n1;i++) pwork[i] = _diagAB[i2prev*_n1+i];
		inda = -i2prev-1+_n2-1;
		indb = _idiag2+i2prev+_n2;
		for (i=0;i<_n1;i++) pwork1[i] = czero;
		CFourier::FastDiagonalPartT2TMultiply (_n1, _idiag1,
															_coefsA+inda*ld1, _coefsB+indb*ld1, pwork1);
		for (i=0;i<_n1;i++) pwork[i] += pwork1[i];
		inda = _n2-i2prev-1+_n2-1;
		indb = _idiag2+i2prev;
		for (i=0;i<_n1;i++) pwork1[i] = czero;
		CFourier::FastDiagonalPartT2TMultiply (_n1, _idiag1,
															_coefsA+inda*ld1, _coefsB+indb*ld1, pwork1);
		for (i=0;i<_n1;i++) pwork[i] -= pwork1[i];
		for (i=0;i<_n1;i++) _diagAB[i2*_n1+i] = pwork[i];
	};

};

// Author: Kharchenko S.A.
// Description: Implement fast diagonal extraction of the 3D Toeplitz products
// CFourier::FastDiagonalPartT2T3DMultiply()
//========================================================================================
void CFourier::FastDiagonalPartT2T3DMultiply (int _n1, int _n2, int _n3, // Implement fast diagonal extraction of the 3D Toeplitz products
																int _idiag1, int _idiag2, int _idiag3,
																dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB, dcmplx *_work) {

	const char *funcname = "FastDiagonalPartT2T3DMultiply";

	dcmplx czero (0.0e0,0.0e0);

// Multiply

	int ld1 = 2*_n1-1;
	int ld2 = 2*_n2-1;
	int ld12 = ld1*ld2;
	int n12 = _n1*_n2;
	dcmplx *pwork, *pwork1, *pwork2;
	int i3, i3prev, j3, i, k, inda, indb, ibeg3, iend3;
	dcmplx caux;

	pwork = _work;
	pwork1 = _work+n12;
	pwork2 = _work+2*n12;

	if (_idiag3 >= 0) {
		ibeg3 = 0;
		iend3 = _n3-_idiag3;
	} else {
		ibeg3 = -_idiag3;
		iend3 = _n3;
	};

	i3 = ibeg3;
	j3 = i3+_idiag3;
	for (i=0;i<n12;i++) pwork[i] = czero;
	for (k=0;k<_n3;k++) {
		inda = k-i3+_n3-1;
		indb = j3-k+_n3-1;
		for (i=0;i<n12;i++) pwork1[i] = czero;
		CFourier::FastDiagonalPartT2T2DMultiply (_n1, _n2, _idiag1, _idiag2,
															_coefsA+inda*ld12, _coefsB+indb*ld12, pwork1, pwork2);
		for (i=0;i<n12;i++) pwork[i] += pwork1[i];
	};
	for (i=0;i<n12;i++) _diagAB[ibeg3*n12+i] = pwork[i];
	for (i3=ibeg3+1;i3<iend3;i3++) {
		i3prev = i3-1;
		for (i=0;i<n12;i++) pwork[i] = _diagAB[i3prev*n12+i];
		inda = -i3prev-1+_n3-1;
		indb = _idiag3+i3prev+_n3;
		for (i=0;i<n12;i++) pwork1[i] = czero;
		CFourier::FastDiagonalPartT2T2DMultiply (_n1, _n2, _idiag1, _idiag2, 
															_coefsA+inda*ld12, _coefsB+indb*ld12, pwork1, pwork2);
		for (i=0;i<n12;i++) pwork[i] += pwork1[i];
		inda = _n3-i3prev-1+_n3-1;
		indb = _idiag3+i3prev;
		for (i=0;i<n12;i++) pwork1[i] = czero;
		CFourier::FastDiagonalPartT2T2DMultiply (_n1, _n2, _idiag1, _idiag2,
															_coefsA+inda*ld12, _coefsB+indb*ld12, pwork1, pwork2);
		for (i=0;i<n12;i++) pwork[i] -= pwork1[i];
		for (i=0;i<n12;i++) _diagAB[i3*n12+i] = pwork[i];
	};

};

// Author: Kharchenko S.A.
// Description: Compute coefficients of the complex conjugate matrix for the Toeplitz matrix
// CFourier::ToeplitzComplexConjugate()
//========================================================================================
void CFourier::ToeplitzComplexConjugate (int _n, // Compute coefficients of the complex conjugate matrix for the Toeplitz matrix
														dcmplx *_coefsA, dcmplx *_coefsAH) {

	const char *funcname = "ToeplitzComplexConjugate";

	int ishift = 2*_n-1;
	int i, ip;
	dcmplx caux;

	for (i=0;i<ishift;i++) {
		ip = ishift-i-1;
		caux = _coefsA[ip];
		_coefsAH[i] = conj(caux);
	};

};

// Author: Kharchenko S.A.
// Description: Compute coefficients of the complex conjugate matrix for the Toeplitz 2D matrix
// CFourier::Toeplitz2DComplexConjugate()
//========================================================================================
void CFourier::Toeplitz2DComplexConjugate (int _n1, int _n2, // Compute coefficients of the complex conjugate matrix for the Toeplitz 2D matrix
															dcmplx *_coefsA, dcmplx *_coefsAH) {

	const char *funcname = "Toeplitz2DComplexConjugate";

	int ishift1 = 2*_n1-1;
	int ishift2 = 2*_n2-1;
	int i1, i2, i1p, i2p, ind, indp;
	dcmplx caux;

	for (i1=0;i1<ishift1;i1++) {
		for (i2=0;i2<ishift2;i2++) {
			i1p = ishift1-i1-1;
			i2p = ishift2-i2-1;
			ind = i2*ishift1+i1;
			indp = i2p*ishift1+i1p;
			caux = _coefsA[indp];
			_coefsAH[ind] = conj(caux);
		};
	};

};

// Author: Kharchenko S.A.
// Description: Compute coefficients of the complex conjugate matrix for the Toeplitz 3D matrix
// CFourier::Toeplitz3DComplexConjugate()
//========================================================================================
void CFourier::Toeplitz3DComplexConjugate (int _n1, int _n2, int _n3, // Compute coefficients of the complex conjugate matrix for the Toeplitz 3D matrix
															dcmplx *_coefsA, dcmplx *_coefsAH) {

	const char *funcname = "Toeplitz3DComplexConjugate";

	int ishift1 = 2*_n1-1;
	int ishift2 = 2*_n2-1;
	int ishift3 = 2*_n3-1;
	int i1, i2, i3, i1p, i2p, i3p, ind, indp;
	dcmplx caux;

	for (i1=0;i1<ishift1;i1++) {
		for (i2=0;i2<ishift2;i2++) {
			for (i3=0;i3<ishift3;i3++) {
				i1p = ishift1-i1-1;
				i2p = ishift2-i2-1;
				i3p = ishift3-i3-1;
				ind = i3*ishift2*ishift1+i2*ishift1+i1;
				indp = i3p*ishift2*ishift1+i2p*ishift1+i1p;
				caux = _coefsA[indp];
				_coefsAH[ind] = conj(caux);
			};
		};
	};

};
