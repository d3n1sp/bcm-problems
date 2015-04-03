//------------------------------------------------------------------------------------------------
// File: complex.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
/*=========================================*/
/*                COMPLEX                  */
/*=========================================*/
#ifndef ___DCMPLXMY___
#define ___DCMPLXMY___

#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

//////////////////////////////////////////////////////////////////
//.class of complex numbers (with member operators);
struct dcmplx_my {
// Data
	double x, y;
// Functions:
	dcmplx_my () { // Zero constructor
		x = 0.0e0;
		y = 0.0e0;
	};
	dcmplx_my (double _x, double _y) { // Data constructor
		x = _x;
		y = _y;
	};
	inline void real (double _x) {x=_x;}; // Init x part
	inline void imag (double _y) {y=_y;}; // Init y part
	inline double real () { return x; }; // Return x part
	inline double imag () { return y; }; // Return y part
	inline dcmplx_my conj () { // Return complex conjugate
		dcmplx_my temp;
		temp.x =  x;
		temp.y = -y;
		return temp;
	};
	inline void MlEq (dcmplx_my __z2) { // *= operator with complex value
		long double temp = x*__z2.y+y*__z2.x;
		x = x*__z2.x-y*__z2.y;
		y = temp;
	};
	inline void PlEq (dcmplx_my __z2) { // += operator with complex value
		x += __z2.x;
		y += __z2.y;
	};
	inline void MnEq (dcmplx_my __z2) { // -= operator with complex value
		x -= __z2.x;
		y -= __z2.y;
	};
	static void OutArrAccV  (std::ostream &_stream, char *_name, int _isize, dcmplx_my *_carr); // Store complex array
	static void ReadArrAccV (std::istream &_stream, int _isize, dcmplx_my *_carr); // Read complex array
};

dcmplx_my   conj (dcmplx_my&);
double      norm (dcmplx_my&);
double      abs  (dcmplx_my&);
//dcmplx_my   comp (double x, double y = 0.);
//dcmplx_my   polar(double _mag, double _angle = 0.);

int         operator==(dcmplx_my&, dcmplx_my&);
int         operator!=(dcmplx_my&, dcmplx_my&);
dcmplx_my&  operator+(dcmplx_my&);
dcmplx_my&  operator-(dcmplx_my&);

dcmplx_my&  operator+=(dcmplx_my&, dcmplx_my);
dcmplx_my&  operator+=(dcmplx_my&, double);
dcmplx_my&  operator-=(dcmplx_my&, dcmplx_my);
dcmplx_my&  operator-=(dcmplx_my&, double);
dcmplx_my&  operator*=(dcmplx_my&, dcmplx_my);
dcmplx_my&  operator*=(dcmplx_my&, double);
dcmplx_my&  operator/=(dcmplx_my&, dcmplx_my);
dcmplx_my&  operator/=(dcmplx_my&, double);

dcmplx_my   operator+(dcmplx_my, dcmplx_my);
dcmplx_my   operator+(double, dcmplx_my);
dcmplx_my   operator+(dcmplx_my, double);
dcmplx_my   operator-(dcmplx_my, dcmplx_my);
dcmplx_my   operator-(double, dcmplx_my);
dcmplx_my   operator-(dcmplx_my, double);
dcmplx_my   operator*(dcmplx_my, dcmplx_my);
//dcmplx_my   operator*(dcmplx_my, const dcmplx_my);
dcmplx_my   operator*(dcmplx_my, double);
dcmplx_my   operator*(double, dcmplx_my);
dcmplx_my   operator/(dcmplx_my, dcmplx_my);
dcmplx_my   operator/(dcmplx_my, double);
dcmplx_my   operator/(double, dcmplx_my);

std::ostream &operator<< (std::ostream &_stream, const dcmplx_my &_z);

inline dcmplx_my conj(dcmplx_my& z) {
	dcmplx_my temp;
	temp.x =  z.x;
	temp.y = -z.y;
	return temp;
};

inline double norm(dcmplx_my& z) {
	return z.x*z.x+z.y*z.y;
};

inline double abs(dcmplx_my& z) {
	return sqrt(z.x*z.x+z.y*z.y);
};

//inline dcmplx_my comp (double x, double y) {
//	dcmplx_my zz; 
//	zz.x = x; 
//	zz.y = y;
//	return zz;
//};

//inline dcmplx_my polar(double _mag, double _angle) {
//	dcmplx_my zz; 
//	zz.x = _mag*cos(_angle); 
//	zz.y = _mag*sin(_angle);
//	return zz;
//};

inline int operator==(dcmplx_my& _z1, dcmplx_my& _z2) {
	return _z1.x == _z2.x&& _z1.y == _z2.y;
};

inline int operator!=(dcmplx_my& _z1, dcmplx_my& _z2) {
	return _z1.x != _z2.x || _z1.y != _z2.y;
};

inline dcmplx_my& operator+(dcmplx_my& z) {
	return z;
};

inline dcmplx_my& operator-(dcmplx_my& z) {
	z.x = -z.x; 
	z.y = -z.y;
	return z;
};

inline dcmplx_my& operator+=(dcmplx_my& z, dcmplx_my __z2) {
	z.x += __z2.x;
	z.y += __z2.y;
	return z;
};

inline dcmplx_my& operator+=(dcmplx_my& z, double __re_val2) {
	z.x += __re_val2;
	return z;
}

inline dcmplx_my& operator-=(dcmplx_my& z, dcmplx_my __z2) {
	z.x -= __z2.x;
	z.y -= __z2.y;
	return z;
};

inline dcmplx_my& operator-=(dcmplx_my& z, double __re_val2) {
	z.x -= __re_val2;
	return z;
};

inline dcmplx_my& operator*=(dcmplx_my& z, dcmplx_my __z2) {
	long double temp = z.x*__z2.y+z.y*__z2.x;
	z.x = z.x*__z2.x-z.y*__z2.y;
	z.y = temp;
	return z;
};

inline dcmplx_my& operator*=(dcmplx_my& z, double __re_val2) {
	z.x *= __re_val2;
	z.y *= __re_val2;
	return z;
};

inline dcmplx_my& operator/=(dcmplx_my& z, dcmplx_my __z2) {
	long double sum_sqrs = norm(__z2);
	dcmplx_my zz = conj(__z2);
	z *= zz;
	z.x /= sum_sqrs;
	z.y /= sum_sqrs;
	return z;
};

inline dcmplx_my& operator/=(dcmplx_my& z, double __re_val2) {
	z.x /= __re_val2;
	z.y /= __re_val2;
	return z;
};

inline dcmplx_my operator+(dcmplx_my __z1, dcmplx_my __z2) {
	dcmplx_my zz; 
	zz.x = __z1.x + __z2.x;
	zz.y = __z1.y + __z2.y;
	return zz;
};

inline dcmplx_my operator+(double __re_val1, dcmplx_my __z2) {
	dcmplx_my zz;
	zz.x = __re_val1 + __z2.x;
	zz.y = __z2.y;
	return zz;
};

inline dcmplx_my operator+(dcmplx_my __z1, double __re_val2) {
	dcmplx_my zz;
	zz.x = __z1.x + __re_val2;
	zz.y = __z1.y;
	return zz;
};

inline dcmplx_my operator-(dcmplx_my __z1, dcmplx_my __z2) {
	dcmplx_my zz;
	zz.x = __z1.x - __z2.x;
	zz.y = __z1.y - __z2.y;
	return zz;
};

inline dcmplx_my operator-(double __re_val1, dcmplx_my __z2) {
	dcmplx_my zz;
	zz.x = __re_val1 - __z2.x;
	zz.y = -__z2.y;
	return  zz;
};

inline dcmplx_my operator-(dcmplx_my __z1, double __re_val2) {
	dcmplx_my zz;
	zz.x = __z1.x - __re_val2;
	zz.y = __z1.y;
	return  zz;
};

inline dcmplx_my operator*(dcmplx_my __z1, dcmplx_my __z2) {
	dcmplx_my zz;
	zz.x = __z1.x*__z2.x - __z1.y*__z2.y;
	zz.y = __z1.x*__z2.y + __z1.y*__z2.x;
	return  zz;
};

//inline dcmplx_my operator*(dcmplx_my __z1, const dcmplx_my __z2) {
//	dcmplx_my zz;
//	zz.x = __z1.x*__z2.x - __z1.y*__z2.y;
//	zz.y = __z1.x*__z2.y + __z1.y*__z2.x;
//	return  zz;
//};

inline dcmplx_my operator*(dcmplx_my __z1, double __re_val2) {
	dcmplx_my zz;
	zz.x = __z1.x*__re_val2;
	zz.y = __z1.y*__re_val2;
	return  zz;
};

inline dcmplx_my operator*(double __re_val1, dcmplx_my __z2) {
	dcmplx_my zz;
	zz.x = __z2.x*__re_val1;
	zz.y = __z2.y*__re_val1;
	return  zz;
};

inline dcmplx_my operator/(dcmplx_my __z1, dcmplx_my __z2) {
	dcmplx_my zz = conj(__z2);
	zz = zz / norm(__z2);
	zz = __z1 * zz;
	return zz;
};

inline dcmplx_my operator/(dcmplx_my __z1, double __re_val2) {
	dcmplx_my zz;
	zz.x = __z1.x/__re_val2;
	zz.y = __z1.y/__re_val2;
	return  zz;
};

inline dcmplx_my operator/(double __re_val1, dcmplx_my __z2) {
	dcmplx_my zz = conj(__z2);
	zz = zz / norm(__z2);
	zz = __re_val1 * zz;
	return zz;
};

inline std::ostream &operator<< (std::ostream &_stream, const dcmplx_my &_z) { // Output complex number
//	_stream << std::left << std::setw(1) << '(' << std::setw(20) << std::setprecision(14) << _z.x;
//	_stream << std::setw(1) << ',' << std::right << std::setw(20) << std::setprecision(14) << _z.y << std::setw(1) << ')';
	_stream << std::setw(1) << '(' << std::setw(20) << std::setprecision(14) << _z.x;
	_stream << std::setw(1) << ',' << std::setw(20) << std::setprecision(14) << _z.y << std::setw(1) << ')';
	return _stream;
};

//========================================================================
// Store complex array
inline void dcmplx_my::OutArrAccV (std::ostream &_stream, char *_name, int _isize, dcmplx_my *_carr) { // Store complex array

	double *darr;

	int i;

	darr = new double [_isize*2];

	for (i=0;i<_isize;i++) {
		darr[i*2] = _carr[i].x;
		darr[i*2+1] = _carr[i].y;
	};

	_stream.write ((char*)darr,sizeof(double)*_isize*2);
	if (_isize != 0) {
		std::cout << _name << " Write Carr: First elem = " << _carr[0] << " Last elem = " << _carr[_isize-1] << std::endl;
	};

	delete [] darr;

};

//========================================================================
// Read complex array
inline void dcmplx_my::ReadArrAccV (std::istream &_stream, int _isize, dcmplx_my *_carr) { // Store complex array

	double *darr;

	int i;

	darr = new double [_isize*2];

	_stream.read ((char*)darr,sizeof(double)*_isize*2);

	for (i=0;i<_isize;i++) {
		_carr[i].x = darr[i*2];
		_carr[i].y = darr[i*2+1];
	};

	if (_isize != 0) {
		std::cout << " Read Carr: First elem = " << _carr[0] << " Last elem = " << _carr[_isize-1] << std::endl;
	};

	delete [] darr;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif

#endif
