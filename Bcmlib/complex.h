/*=========================================*/
/*                COMPLEX                  */
/*=========================================*/
#ifndef ___COMPLEX___
#define ___COMPLEX___

#include <math.h>

/////////////////////////////////////////////////////////////////////////////
//...constants rounded for 21 decimals for use with class of complex numbers;
#undef  M_E
#define M_E         2.71828182845904523536
#undef  M_LOG2E
#define M_LOG2E     1.44269504088896340736
#undef  M_LOG10E
#define M_LOG10E    0.434294481903251827651
#undef  M_LN2
#define M_LN2       0.693147180559945309417
#undef  M_LN10
#define M_LN10      2.30258509299404568402
#undef  M_PI
#define M_PI        3.14159265358979323846
#undef  M_PI_2
#define M_PI_2      1.57079632679489661923
#undef  M_PI_4
#define M_PI_4      0.785398163397448309616
#undef  M_1_PI
#define M_1_PI      0.318309886183790671538
#undef  M_2_PI
#define M_2_PI      0.636619772367581343076
#undef  M_1_SQRTPI
#define M_1_SQRTPI  0.564189583547756286948
#undef  M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257390
#undef  M_SQRT2
#define M_SQRT2     1.41421356237309504880
#undef  M_SQRT_2
#define M_SQRT_2    0.707106781186547524401
#undef  M_G
#define M_G         9.81

//////////////////////////////////////////////////////////////////
//...class for support of complex numbers (with member operators);
struct  complex {
        double x, y;
        complex () {
            x = y = 0.;
        };
        complex (double __x, double __y = 0.) {
            x = __x;
            y = __y;
        };
        inline void real (double __x) {x = __x;};
        inline void imag (double __y) {y = __y;};
        inline double real () { return x; }; 
        inline double imag () { return y; };
        
        inline complex conj() {
            complex __zz; __zz.x =  x; __zz.y = -y;
            return  __zz;
        };

        inline void Mult (complex __z2) {
            long double __zz_imag = x*__z2.y+y*__z2.x;
				x = x*__z2.x-y*__z2.y;
				y = (double)__zz_imag;
        };
        inline void Plus (complex __z2) {
            x += __z2.x;
            y += __z2.y;
        };
        inline void Minus (complex __z2) {
            x -= __z2.x;
            y -= __z2.y;
        };
};

////////////////////////////////////////////////////////////////
//...auxilliary function for preparing class of complex numbers;
#undef sqr
inline long double  sqr(long double x) { return x*x;}
inline      double  sqr(     double x) { return x*x;}
inline      int	  sqr(		  int x) { return x*x;}

////////////////////////////////////////////////////////////////////////////////////
//...different realization another functions for preparing class of complex numbers;
#ifdef __STANDART_REALIZATION___
double  real(complex& __z);
double  imag(complex& __z);
double  norm(complex& __z);
double  abs (complex& __z);
double  fabs(complex& __z);
int operator==(complex& __z1, complex& __z2);
int operator!=(complex& __z1, complex& __z2);

inline double  real(complex& __z) 
{ 
    return __z.x;
}
inline double  imag(complex& __z) 
{ 
    return __z.y;
}
inline double  norm(complex& __z) 
{ 
    return sqr(__z.x)+sqr(__z.y);
}
inline double  abs (complex& __z) 
{ 
    return (double)sqrt(sqr(__z.x)+sqr(__z.y));
}
inline double fabs (complex& __z) 
{ 
    return (double)sqrt(sqr(__z.x)+sqr(__z.y));
}
inline int operator==(complex& __z1, complex& __z2) 
{ 
    return __z1.x == __z2.x && __z1.y == __z2.y;
}
inline int operator!=(complex& __z1, complex& __z2) 
{ 
    return __z1.x != __z2.x || __z1.y != __z2.y;
}
complex  comp (double __x,   double __y = 0.);
complex  polar(double __mag, double __angle = 0.);
complex  conj (complex&);

inline complex comp(double __x, double __y)
{
    complex __zz; __zz.x = __x; __zz.y = __y;
    return __zz;
}
inline complex polar(double __mag, double __angle)
{
    complex __zz; __zz.x = __mag*cos(__angle); __zz.y = __mag*sin(__angle);
    return __zz;
}
inline complex conj(complex& __z)
{
    complex __zz; __zz.x = __z.x; __zz.y = -__z.y;
    return __zz;
}
complex  operator-(complex&);
complex  operator+(complex&, complex&);
complex  operator+(double, complex&);
complex  operator+(complex&, double);
complex  operator-(complex&, complex&);
complex  operator-(double, complex&);
complex  operator-(complex&, double);
complex  operator*(complex&, complex&);
complex  operator*(complex&, double);
complex  operator*(double, complex&);
complex  operator/(complex&, complex&);
complex  operator/(complex&, double);
complex  operator/(double, complex&);

inline complex operator-(complex& __z)
{
    complex __zz; __zz.x = -__z.x; __zz.y = -__z.y;
    return __zz;
}
inline complex operator+(complex& __z1, complex& __z2)
{
    complex __zz; __zz.x = __z1.x + __z2.x; __zz.y = __z1.y + __z2.y;
    return __zz;
}
inline complex operator+(double __re_val1, complex& __z2)
{
    complex __zz; __zz.x = __re_val1 + __z2.x; __zz.y = __z2.y;
    return __zz;
}
inline complex operator+(complex& __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x + __re_val2; __zz.y = __z1.y;
    return __zz;
}
inline complex operator-(complex& __z1, complex& __z2)
{
    complex __zz; __zz.x = __z1.x - __z2.x; __zz.y = __z1.y - __z2.y;
    return __zz;
}
inline complex operator-(double __re_val1, complex& __z2)
{
    complex __zz; __zz.x = __re_val1 - __z2.x; __zz.y = -__z2.y;
    return __zz;
}
inline complex operator-(complex& __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x - __re_val2; __zz.y = __z1.y;
    return __zz;
}
inline complex operator*(complex& __z1, complex& __z2)
{
    complex __zz; __zz.x = __z1.x*__z2.x - __z1.y*__z2.y; __zz.y = __z1.x*__z2.y + __z1.y*__z2.x;
    return __zz;
}
inline complex operator*(complex& __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x*__re_val2; __zz.y = __z1.y*__re_val2;
    return __zz;
}
inline complex operator*(double __re_val1, complex& __z2)
{
    complex __zz; __zz.x = __z2.x*__re_val1; __zz.y = __z2.y*__re_val1;
    return __zz;
}
inline complex operator/(complex& __z1, complex& __z2)
{
    complex __zz = conj(__z2);
    __zz = __zz / norm(__z2);
    __zz = __z1 * __zz;
    return __zz;
}
inline complex operator/(complex& __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x/__re_val2; __zz.y = __z1.y/__re_val2;
    return __zz;
}
inline complex operator/(double __re_val1, complex& __z2)
{
    complex __zz = conj(__z2);
    __zz = __zz / norm(__z2);
    __zz = __re_val1 * __zz;
    return __zz;
}
complex& operator+ (complex&);
complex& operator+=(complex&, complex&);
complex& operator+=(complex&, double);
complex& operator-=(complex&, complex&);
complex& operator-=(complex&, double);
complex& operator*=(complex&, complex&);
complex& operator*=(complex&, double);
complex& operator/=(complex&, complex&);
complex& operator/=(complex&, double);

inline complex& operator +(complex& __z)
{ 
    return __z; 
}
inline complex& operator+=(complex& z, complex& __z2)
{
    z.x += __z2.x;
    z.y += __z2.y;
    return z;
}
inline complex& operator+=(complex& z, double __re_val2)
{
    z.x += __re_val2;
    return z;
}
inline complex& operator-=(complex& z, complex& __z2)
{
    z.x -= __z2.x;
    z.y -= __z2.y;
    return z;
}
inline complex& operator-=(complex& z, double __re_val2)
{
    z.x -= __re_val2;
    return z;
}
inline complex& operator*=(complex& z, complex& __z2)
{
    long double temp = z.x*__z2.y+z.y*__z2.x;
    z.x = z.x*__z2.x-z.y*__z2.y;
    z.y = temp;
    return z;
}
inline complex& operator*=(complex& z, double __re_val2)
{
    z.x *= __re_val2;
    z.y *= __re_val2;
    return z;
}
inline complex& operator/=(complex& z, complex& __z2)
{
    long double sum_sqrs = norm(__z2);
    complex __zz = conj(__z2);
    z *= __zz;
    z.x /= sum_sqrs;
    z.y /= sum_sqrs;
    return z;
}
inline complex& operator/=(complex& z, double __re_val2)
{
    z.x /= __re_val2;
    z.y /= __re_val2;
    return z;
}
#else
double  real(complex __z);
double  imag(complex __z);
double  norm(complex __z);
double  abs (complex __z);
double  fabs(complex __z);
int operator==(complex __z1, complex __z2);
int operator!=(complex __z1, complex __z2);

inline double  real(complex __z) 
{ 
    return __z.x;
}
inline double  imag(complex __z) 
{ 
    return __z.y;
}
inline double  norm(complex __z) 
{ 
	 return (double)(sqr(__z.x)+sqr(__z.y));
}
inline double  abs (complex __z) 
{ 
    return sqrt(sqr(__z.x)+sqr(__z.y));
}
inline double fabs (complex __z) 
{ 
    return sqrt(sqr(__z.x)+sqr(__z.y));
}
inline int operator==(complex __z1, complex __z2) 
{ 
    return __z1.x == __z2.x && __z1.y == __z2.y;
}
inline int operator!=(complex __z1, complex __z2) 
{ 
    return __z1.x != __z2.x || __z1.y != __z2.y;
}
complex  comp (double __x,   double __y = 0.);
complex  polar(double __mag, double __angle = 0.);
complex  conj (complex);

inline complex comp(double __x, double __y)
{
    complex __zz; __zz.x = __x; __zz.y = __y;
    return __zz;
}
inline complex polar(double __mag, double __angle)
{
    complex __zz; __zz.x = __mag*cos(__angle); __zz.y = __mag*sin(__angle);
    return __zz;
}
inline complex conj(complex __z)
{
    complex __zz; __zz.x = __z.x; __zz.y = -__z.y;
    return __zz;
}
complex  operator-(complex);
complex  operator+(complex, complex);
complex  operator+(double, complex);
complex  operator+(complex, double);
complex  operator-(complex, complex);
complex  operator-(double, complex);
complex  operator-(complex, double);
complex  operator*(complex, complex);
complex  operator*(complex, double);
complex  operator*(double, complex);
complex  operator/(complex, complex);
complex  operator/(complex, double);
complex  operator/(double, complex);

inline complex operator-(complex __z)
{
    complex __zz; __zz.x = -__z.x; __zz.y = -__z.y;
    return __zz;
}
inline complex operator+(complex __z1, complex __z2)
{
    complex __zz; __zz.x = __z1.x + __z2.x; __zz.y = __z1.y + __z2.y;
    return __zz;
}
inline complex operator+(double __re_val1, complex __z2)
{
    complex __zz; __zz.x = __re_val1 + __z2.x; __zz.y = __z2.y;
    return __zz;
}
inline complex operator+(complex __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x + __re_val2; __zz.y = __z1.y;
    return __zz;
}
inline complex operator-(complex __z1, complex __z2)
{
    complex __zz; __zz.x = __z1.x - __z2.x; __zz.y = __z1.y - __z2.y;
    return __zz;
}
inline complex operator-(double __re_val1, complex __z2)
{
    complex __zz; __zz.x = __re_val1 - __z2.x; __zz.y = -__z2.y;
    return __zz;
}
inline complex operator-(complex __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x - __re_val2; __zz.y = __z1.y;
    return __zz;
}
inline complex operator*(complex __z1, complex __z2)
{
    complex __zz; __zz.x = __z1.x*__z2.x - __z1.y*__z2.y; __zz.y = __z1.x*__z2.y + __z1.y*__z2.x;
    return __zz;
}
inline complex operator*(complex __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x*__re_val2; __zz.y = __z1.y*__re_val2;
    return __zz;
}
inline complex operator*(double __re_val1, complex __z2)
{
    complex __zz; __zz.x = __z2.x*__re_val1; __zz.y = __z2.y*__re_val1;
    return __zz;
}
inline complex operator/(complex __z1, complex __z2)
{
    complex  __zz = conj(__z2);
    __zz = __zz / norm(__z2);
    __zz = __z1 * __zz;
    return __zz;
}
inline complex operator/(complex __z1, double __re_val2)
{
    complex __zz; __zz.x = __z1.x/__re_val2; __zz.y = __z1.y/__re_val2;
    return __zz;
}
inline complex operator/(double __re_val1, complex __z2)
{
    complex  __zz = conj(__z2);
    __zz = __zz / norm(__z2);
    __zz = __re_val1 * __zz;
    return __zz;
}
complex& operator+ (complex&);
complex& operator+=(complex&, complex);
complex& operator+=(complex&, double);
complex& operator-=(complex&, complex);
complex& operator-=(complex&, double);
complex& operator*=(complex&, complex);
complex& operator*=(complex&, double);
complex& operator/=(complex&, complex);
complex& operator/=(complex&, double);

inline complex& operator +(complex& __z)
{ 
    return __z; 
}
inline complex& operator+=(complex& z, complex __z2)
{
    z.x += __z2.x;
    z.y += __z2.y;
    return z;
}
inline complex& operator+=(complex& z, double __re_val2)
{
    z.x += __re_val2;
    return z;
}
inline complex& operator-=(complex& z, complex __z2)
{
    z.x -= __z2.x;
    z.y -= __z2.y;
    return z;
}
inline complex& operator-=(complex& z, double __re_val2)
{
    z.x -= __re_val2;
    return z;
}
inline complex& operator*=(complex& z, complex __z2)
{
    long double temp = z.x*__z2.y+z.y*__z2.x;
    z.x = z.x*__z2.x-z.y*__z2.y;
	 z.y = (double)temp;
	 return z;
}
inline complex& operator*=(complex& z, double __re_val2)
{
    z.x *= __re_val2;
    z.y *= __re_val2;
    return z;
}
inline complex& operator/=(complex& z, complex __z2)
{
    long double sum_sqrs = norm(__z2);
    complex __zz = conj(__z2);
    z *= __zz;
	 z.x = (double)(z.x/sum_sqrs);
	 z.y = (double)(z.y/sum_sqrs);
	 return z;
}
inline complex& operator/=(complex& z, double __re_val2)
{
    z.x /= __re_val2;
    z.y /= __re_val2;
    return z;
}
#endif
//////////////////////////////////////
//...;
inline double to_double(complex a) { return real(a);}
inline int to_int(complex a) { return static_cast<int>(real(a));}
#endif

