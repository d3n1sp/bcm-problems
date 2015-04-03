#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "globals.h"

using namespace std;

//========================================================================
// Standard 4x4 matrix by matrix complex multiplication
void MtrXMtrStd4C (int _nmatr, dcmplx *_x, dcmplx *_y, dcmplx *_xy) { // Standard 4x4 matrix by matrix complex multiplication

//	const char *funcname = "MtrXMtrStd4C";

	dcmplx czero (0.0e0,0.0e0);

	int bla = 4;
	int bla_2 = 16;

	int imatr, kii, kjj, kkk;

	dcmplx caux;

	for (imatr=0;imatr<_nmatr;imatr++) {

		dcmplx *xloc = _x + imatr*bla_2;
		dcmplx *yloc = _y + imatr*bla_2;

		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				caux = czero;
				for (kkk=0;kkk<bla;kkk++) {
					caux += xloc[kii*bla+kkk] * yloc[kjj*bla+kkk];
				};
				_xy[kii*bla+kjj] += caux;
			};
		};

	};

};

//========================================================================
// Loop unroll 4x4 matrix by matrix complex multiplication
void MtrXMtrLoop4C (int _nmatr, dcmplx *_x, dcmplx *_y, dcmplx *_xy) { // Loop unroll 4x4 matrix by matrix complex multiplication

//	const char *funcname = "MtrXMtrLoop4C";

	dcmplx czero (0.0e0,0.0e0);

	int bla = 4;
	int bla_2 = 16;

	int imatr, kii, kjj, kkk, kiib/*, kjjb*/;

	dcmplx *xloc;
	dcmplx *yloc;

	dcmplx *xxloc, *yyloc, *xyloc;

	dcmplx caux;

	for (imatr=0;imatr<_nmatr;imatr++) {

		int ishft = imatr*bla_2;

		xyloc = _xy;

		xloc = _x + ishft;

		kiib = 0;
		for (kii=0;kii<bla;kii++) {
			yloc = _y + ishft;
			yyloc = yloc;
			for (kjj=0;kjj<bla;kjj++) {
				xxloc = xloc+kiib;
//				caux = czero;
				for (kkk=0;kkk<bla;kkk++) {
					*xyloc += *xxloc++ * *yyloc++;
//					caux += *xxloc++ * *yyloc++;
				};
//				*xyloc = caux;
				xyloc++;
			};
			kiib += bla;
		};

	};

};

//========================================================================
// Standard 8x8 matrix by matrix complex multiplication
void MtrXMtrStd8C (int _nmatr, dcmplx *_x, dcmplx *_y, dcmplx *_xy) { // Standard 8x8 matrix by matrix complex multiplication

//	const char *funcname = "MtrXMtrStd8C";

	dcmplx czero (0.0e0,0.0e0);

	int bla = 8;
	int bla_2 = 64;

	int imatr, kii, kjj, kkk;

	dcmplx caux;

	for (imatr=0;imatr<_nmatr;imatr++) {

		dcmplx *xloc = _x + imatr*bla_2;
		dcmplx *yloc = _y + imatr*bla_2;

		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				caux = czero;
				for (kkk=0;kkk<bla;kkk++) {
					caux += xloc[kii*bla+kkk] * yloc[kjj*bla+kkk];
				};
				_xy[kii*bla+kjj] += caux;
			};
		};

	};

};

//========================================================================
// Loop unroll 8x8 matrix by matrix complex multiplication
void MtrXMtrLoop8C (int _nmatr, dcmplx *_x, dcmplx *_y, dcmplx *_xy) { // Loop unroll 8x8 matrix by matrix complex multiplication

//	const char *funcname = "MtrXMtrLoop8C";

	dcmplx czero (0.0e0,0.0e0);

	int bla = 8;
	int bla_2 = 64;

	int imatr, kii, kjj, kkk, kiib, kjjb;

	dcmplx *xloc;
	dcmplx *yloc;

	dcmplx *xxloc, *yyloc, *xyloc;

	dcmplx caux;

	for (imatr=0;imatr<_nmatr;imatr++) {

		int ishft = imatr*bla_2;

		xyloc = _xy;

		xloc = _x + ishft;

		kiib = 0;
		for (kii=0;kii<bla;kii++) {
			yloc = _y + ishft;
			kjjb = 0;
			for (kjj=0;kjj<bla;kjj++) {
				xxloc = xloc+kiib;
				yyloc = yloc+kjjb;
				caux = czero;
				for (kkk=0;kkk<bla;kkk++) {
					caux += *xxloc++ * *yyloc++;
				};
				*xyloc += caux;
				xyloc++;
				kjjb += bla;
			};
			kiib += bla;
		};

	};

};

//========================================================================
// Standard 8x8 matrix by matrix complex multiplication
void MtrXMtrStd8R (int _nmatr, double *_x, double *_y, double *_xy) { // Standard 8x8 matrix by matrix real multiplication

//	const char *funcname = "MtrXMtrStd8R";

	double dzero = 0.0e0;

	int bla = 8;
	int bla_2 = 64;

	int imatr, kii, kjj, kkk;

	double daux;

	for (imatr=0;imatr<_nmatr;imatr++) {

		double *xloc = _x + imatr*bla_2;
		double *yloc = _y + imatr*bla_2;

		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				daux = dzero;
				for (kkk=0;kkk<bla;kkk++) {
					daux += xloc[kii*bla+kkk] * yloc[kjj*bla+kkk];
				};
				_xy[kii*bla+kjj] += daux;
			};
		};

	};

};

//========================================================================
// Loop unroll 8x8 matrix by matrix real multiplication
void MtrXMtrLoop8R (int _nmatr, double *_x, double *_y, double *_xy) { // Loop unroll 8x8 matrix by matrix real multiplication

//	const char *funcname = "MtrXMtrLoop8R";

	double dzero = 0.0e0;

	int bla = 8;
	int bla_2 = 64;

	int imatr, kii, kjj, kkk, kiib, kjjb;

	double *xloc;
	double *yloc;

	double daux;

	double *xxloc, *yyloc, *xyloc;

	for (imatr=0;imatr<_nmatr;imatr++) {

		int ishft = imatr*bla_2;

		xyloc = _xy;

		xloc = _x + ishft;

		kiib = 0;
		for (kii=0;kii<bla;kii++) {
			yloc = _y + ishft;
			kjjb = 0;
			for (kjj=0;kjj<bla;kjj++) {
				xxloc = xloc+kiib;
				yyloc = yloc+kjjb;
				daux = dzero;
				for (kkk=0;kkk<bla;kkk++) {
//					*xyloc += *xxloc++ * *yyloc++;
					daux += *xxloc++ * *yyloc++;
				};
				*xyloc += daux;
				xyloc++;
				kjjb += bla;
			};
			kiib += bla;
		};

	};

};

//========================================================================
// Full unroll 4x4 matrix by matrix complex multiplication
void MtrXMtrUnr4C (int _nmatr, dcmplx *_x, dcmplx *_y, dcmplx *_xy) { // Full unroll 4x4 matrix by matrix complex multiplication

//	const char *funcname = "MtrXMtrUnr4C";

	int bla = 4;
	int bla_2 = 16;

	int imatr;

	dcmplx *xloc;
	dcmplx *yloc;

	dcmplx *xxloc, *yyloc;

	dcmplx caux0, caux1, caux2, caux3;

	for (imatr=0;imatr<_nmatr;imatr++) {

		int ishft = imatr*bla_2;

		xloc = _x + ishft;
		yloc = _y + ishft;

		xxloc = xloc;
		yyloc = yloc;

		caux0 = xxloc[ 0]*yyloc[ 0] + xxloc[ 1]*yyloc[ 1] + xxloc[ 2]*yyloc[ 2] + xxloc[ 3]*yyloc[ 3];
		caux1 = xxloc[ 0]*yyloc[ 4] + xxloc[ 1]*yyloc[ 5] + xxloc[ 2]*yyloc[ 6] + xxloc[ 3]*yyloc[ 7];
		caux2 = xxloc[ 0]*yyloc[ 8] + xxloc[ 1]*yyloc[ 9] + xxloc[ 2]*yyloc[10] + xxloc[ 3]*yyloc[11];
		caux3 = xxloc[ 0]*yyloc[12] + xxloc[ 1]*yyloc[13] + xxloc[ 2]*yyloc[14] + xxloc[ 3]*yyloc[15];

		_xy[ 0] += caux0;
		_xy[ 1] += caux1;
		_xy[ 2] += caux2;
		_xy[ 3] += caux3;

		caux0 = xxloc[ 4]*yyloc[ 0] + xxloc[ 5]*yyloc[ 1] + xxloc[ 6]*yyloc[ 2] + xxloc[ 7]*yyloc[ 3];
		caux1 = xxloc[ 4]*yyloc[ 4] + xxloc[ 5]*yyloc[ 5] + xxloc[ 6]*yyloc[ 6] + xxloc[ 7]*yyloc[ 7];
		caux2 = xxloc[ 4]*yyloc[ 8] + xxloc[ 5]*yyloc[ 9] + xxloc[ 6]*yyloc[10] + xxloc[ 7]*yyloc[11];
		caux3 = xxloc[ 4]*yyloc[12] + xxloc[ 5]*yyloc[13] + xxloc[ 6]*yyloc[14] + xxloc[ 7]*yyloc[15];

		_xy[ 4] += caux0;
		_xy[ 5] += caux1;
		_xy[ 6] += caux2;
		_xy[ 7] += caux3;

		caux0 = xxloc[ 8]*yyloc[ 0] + xxloc[ 9]*yyloc[ 1] + xxloc[10]*yyloc[ 2] + xxloc[11]*yyloc[ 3];
		caux1 = xxloc[ 8]*yyloc[ 4] + xxloc[ 9]*yyloc[ 5] + xxloc[10]*yyloc[ 6] + xxloc[11]*yyloc[ 7];
		caux2 = xxloc[ 8]*yyloc[ 8] + xxloc[ 9]*yyloc[ 9] + xxloc[10]*yyloc[10] + xxloc[11]*yyloc[11];
		caux3 = xxloc[ 8]*yyloc[12] + xxloc[ 9]*yyloc[13] + xxloc[10]*yyloc[14] + xxloc[11]*yyloc[15];

		_xy[ 8] += caux0;
		_xy[ 9] += caux1;
		_xy[10] += caux2;
		_xy[11] += caux3;

		caux0 = xxloc[12]*yyloc[ 0] + xxloc[13]*yyloc[ 1] + xxloc[14]*yyloc[ 2] + xxloc[15]*yyloc[ 3];
		caux1 = xxloc[12]*yyloc[ 4] + xxloc[13]*yyloc[ 5] + xxloc[14]*yyloc[ 6] + xxloc[15]*yyloc[ 7];
		caux2 = xxloc[12]*yyloc[ 8] + xxloc[13]*yyloc[ 9] + xxloc[14]*yyloc[10] + xxloc[15]*yyloc[11];
		caux3 = xxloc[12]*yyloc[12] + xxloc[13]*yyloc[13] + xxloc[14]*yyloc[14] + xxloc[15]*yyloc[15];

		_xy[12] += caux0;
		_xy[13] += caux1;
		_xy[14] += caux2;
		_xy[15] += caux3;

	};

};

//========================================================================
// Test unroll computations
void TestUnroll () { // Test unroll computations

	const char *funcname = "TestUnroll";

	dcmplx czero (0.0e0,0.0e0);
	double dzero = 0.0e0;

	double tottim;
	clock_t time0, time1;

	int nmatr = 300000;
	int bla = 4;

	int bla_2 = bla*bla;

	dcmplx *x;
	dcmplx *y;
	dcmplx *xy;

	x = new dcmplx [nmatr*bla_2];
	if (!x) MemoryFail(funcname);
	y = new dcmplx [nmatr*bla_2];
	if (!y) MemoryFail(funcname);
	xy = new dcmplx [bla_2];
	if (!xy) MemoryFail(funcname);

	double *xr;
	double *yr;
	double *xyr;

	xr = new double [nmatr*bla_2];
	if (!x) MemoryFail(funcname);
	yr = new double [nmatr*bla_2];
	if (!y) MemoryFail(funcname);
	xyr = new double [bla_2];
	if (!xy) MemoryFail(funcname);

// Init x and y

	int imatr, kii, kjj;

	for (imatr=0;imatr<nmatr;imatr++) {
		for (kii=0;kii<bla;kii++) {
			for (kjj=0;kjj<bla;kjj++) {
				x[imatr*bla_2+kii+kjj*bla].real(1.0e0/((double) (kii-kjj+imatr) + 0.43e0));
				x[imatr*bla_2+kii+kjj*bla].imag(-1.0e0/((double) (kjj-kii+imatr) + 0.57e0));
				y[imatr*bla_2+kii+kjj*bla].real(1.0e0/((double) (kii-kjj-imatr) + 0.47e0));
				y[imatr*bla_2+kii+kjj*bla].imag(-1.0e0/((double) (kjj-kii-imatr) + 0.53e0));
				xr[imatr*bla_2+kii+kjj*bla] = 1.0e0/((double) (kii-kjj+imatr) + 0.43e0);
				yr[imatr*bla_2+kii+kjj*bla] = 1.0e0/((double) (kii-kjj-imatr) + 0.57e0);
			};
		};
	};

// Perform tree types of computations with timing and compare the results

	ofstream fout ("Results.dat");

	for (kii=0;kii<bla_2;kii++) xy[kii] = czero;
	for (kii=0;kii<bla_2;kii++) xyr[kii] = dzero;

	time0 = clock ();

	MtrXMtrStd4C (nmatr, x, y, xy);
//	MtrXMtrStd8R (nmatr, xr, yr, xyr);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	cout  << "  Time0  = " << tottim << " sec." << endl;

//	OutArr (fout," Std ",bla_2,xy);

	for (kii=0;kii<bla_2;kii++) xy[kii] = czero;
	for (kii=0;kii<bla_2;kii++) xyr[kii] = dzero;

	time0 = clock ();

	MtrXMtrLoop4C (nmatr, x, y, xy);
//	MtrXMtrLoop8R (nmatr, xr, yr, xyr);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	cout  << "  Time1  = " << tottim << " sec." << endl;

	OutArr (fout," Loop ",bla_2,xy);

	for (kii=0;kii<bla_2;kii++) xy[kii] = czero;

	time0 = clock ();

	MtrXMtrUnr4C (nmatr, x, y, xy);

	time1 = clock ();

	tottim = (double) (time1-time0) / (double) CLOCKS_PER_SEC;

	cout  << "  Time2  = " << tottim << " sec." << endl;

	OutArr (fout," Unr ",bla_2,xy);

};
