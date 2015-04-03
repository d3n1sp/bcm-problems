// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "globals.h"

// Fourier.h: Description of the Fourier transform
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Fourier
#define __Fourier

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFourier
{
public:
// Functions
// Constructors and destructor
	CFourier () {}; // Memory allocation zero data constructor
	~CFourier () { // Destructor
//		cout << " On entry to CFourier destructor " << endl;
//		cout << " On return from CFourier destructor " << endl;
	};
	CFourier (const CFourier &_func) {}; // Copy constructor
// Operator functions
	CFourier &operator= (const CFourier &_func) {return *this;}; // Equality operator
// Toeplitz functions
	static void FastToeplitz3DMvm (int _n1, int _n2, int _n3, int _ncols, // Implement fast Toeplitz 3D matrix by vector multiplication
												dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void FastToeplitz3DMvmH (int _n1, int _n2, int _n3, int _ncols, // Implement fast Toeplitz 3D hermitian matrix by vector multiplication
												dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExplicitToeplitz3DMvm (int _n1, int _n2, int _n3, int _ncols, // Implement explicit Toeplitz 3D matrix by vector multiplication
													dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExpandToeplitz3D (int _n1, int _n2, int _n3, // Expand the Toeplitz 3D matrix
											dcmplx *_coefs, dcmplx *_amatr);
	static void FastToeplitz2DMvm (int _n1, int _n2, int _ncols, // Implement fast Toeplitz 2D matrix by vector multiplication
												dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void FastToeplitz2DMvmH (int _n1, int _n2, int _ncols, // Implement fast Toeplitz 2D hermitian matrix by vector multiplication
												dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExplicitToeplitz2DMvm (int _n1, int _n2, int _ncols, // Implement explicit Toeplitz 2D matrix by vector multiplication
													dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExpandToeplitz2D (int _n1, int _n2, // Expand the Toeplitz 2D matrix
											dcmplx *_coefs, dcmplx *_amatr);
	static void FastToeplitzMvm (int _n, int _ncols, // Implement fast Toeplitz matrix by vector multiplication
											dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void FastToeplitzMvmH (int _n, int _ncols, // Implement fast Toeplitz hermitian matrix by vector multiplication
											dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExplicitToeplitzMvm (int _n, int _ncols, // Implement explicit Toeplitz matrix by vector multiplication
												dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExpandToeplitz (int _n, // Expand the Toeplitz matrix
											dcmplx *_coefs, dcmplx *_amatr);
// Extract diagonal part of the Toeplitz type matrices products
	static void FastDiagonalPartT2T3DMultiply (int _n1, int _n2, int _n3, // Implement fast diagonal extraction of the 3D Toeplitz products
																int _idiag1, int _idiag2, int _idiag3,
																dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB, dcmplx *_work);
	static void FastDiagonalPartT2T2DMultiply (int _n1, int _n2, int _idiag1, int _idiag2, // Implement fast diagonal extraction of the 2D Toeplitz products
																dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB, dcmplx *_work);
	static void FastDiagonalPartT2TMultiply (int _n, int _idiag, // Implement fast diagonal extraction of the Toeplitz products
															dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB);
	static void ExplicitDiagonalPartT2TMultiply (int _n, int _idiag, // Implement explicit diagonal extraction of the Toeplitz products
																dcmplx *_coefsA, dcmplx *_coefsB, dcmplx *_diagAB);
// Compute complex conjugate matrix for the Toeplitz type matrices
	static void Toeplitz3DComplexConjugate (int _n1, int _n2, int _n3, // Compute coefficients of the complex conjugate matrix for the Toeplitz 3D matrix
															dcmplx *_coefsA, dcmplx *_coefsAH);
	static void Toeplitz2DComplexConjugate (int _n1, int _n2, // Compute coefficients of the complex conjugate matrix for the Toeplitz 2D matrix
															dcmplx *_coefsA, dcmplx *_coefsAH);
	static void ToeplitzComplexConjugate (int _n, // Compute coefficients of the complex conjugate matrix for the Toeplitz matrix
														dcmplx *_coefsA, dcmplx *_coefsAH);
// Circulant multiplication functions
	static void FastCirculantMvm (int _n, int _ncols, // Implement fast circulant matrix by vector multiplication
												dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
	static void ExplicitCirculantMvm (int _n, int _ncols, // Implement explicit circulant matrix by vector multiplication
													dcmplx *_coefs, dcmplx *_x, dcmplx *_fx);
// Fourier functions
	static void FFTPower2 (int _n, int _ncols, // Use the fastest Fast Fourier transform algorithm for the power of 2
									dcmplx *_x, dcmplx *_fx);
	static void FFTPower2Impl (int _n, int _ncols, // Implement Fast Fourier transform for the power of 2
									dcmplx *_x, dcmplx *_fx);
	static void ExplicitFourier (char _descr, int _n, int _n2, int _ncols, // Implement explicit Fourier transform by columns/rows
											dcmplx *_x, dcmplx *_fx);
	static void FTF1C (double *_ar0, double *_ai0, // Implement Fast Fourier transform for the power of 2
								int _n, int _in, int _k, double _p);
private:
	static void ScalingFourier (int _n1, int _n2, int _ncols, // Perform scaling of the matrix
							dcmplx *_cepsarr, 
							dcmplx *_x);
	static void TransposeFourier (int _n1, int _n2, int _ncols, // Perform transposition of the matrix
							dcmplx *_x, 
							dcmplx *_fx);
	static void ImplicitTwoLevelFourier (int _n, int _ncols, int _ndivide, // Implement implicit two level Fourier transform by columns
										dcmplx *_x,
										dcmplx *_fx);
public:
// Input/Output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CFourier &_func) {return _stream;}; // Output function description
// Friend classes
};

#endif

