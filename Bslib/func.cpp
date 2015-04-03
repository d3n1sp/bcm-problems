#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "func.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFunc::CFunc()
//========================================================================================
CFunc::CFunc () { // Memory allocation zero data constructor

//	const char *funcname = "CFunc_00";

	c0 = 0.0e0;
	c1 = 0.0e0;
	c2 = 0.0e0;
	c3 = 0.0e0;
	c_1 = 0.0e0;
	cadd_1 = 0.0e0;

};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CFunc::CFunc()
//========================================================================================
CFunc::CFunc (const CFunc &_func) { // Copy constructor

//	const char *funcname = "CFunc_copy";

	c0 = _func.c0;
	c1 = _func.c1;
	c2 = _func.c2;
	c3 = _func.c3;
	c_1 = _func.c_1;
	cadd_1 = _func.cadd_1;

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CFunc::operator=()
//========================================================================================
CFunc &CFunc::operator= (const CFunc &_func) { // Equality operator

//	const char *funcname = "CFunc_=";

	c0 = _func.c0;
	c1 = _func.c1;
	c2 = _func.c2;
	c3 = _func.c3;
	c_1 = _func.c_1;
	cadd_1 = _func.cadd_1;

	return *this;

};

// Author: Kharchenko S.A.
// Description: Compute function value
// CFunc::ComputeValue()
//========================================================================================
double CFunc::ComputeValue (double _value) { // Compute function value

//	const char *funcname = "ComputeValue";

	double fvalue = 0.0e0;

	fvalue += c0;

	fvalue += c1 * _value;

	double daux = _value * _value;
	fvalue += c2 * daux;

	daux *= _value;
	fvalue += c3 * daux;

	if (c_1 != 0.0e0) {
		daux = cadd_1+_value;
		fvalue += c_1 / daux;
	};

	return fvalue;

};

// Author: Kharchenko S.A.
// Description: Compute derivative value
// CFunc::ComputeDerivative()
//========================================================================================
double CFunc::ComputeDerivative (double _value) { // Compute derivative value

//	const char *funcname = "ComputeDerivative";

	double fvalue = 0.0e0;

	fvalue += c1;

	fvalue += 2.0e0 * c2 * _value;

	double daux = _value * _value;
	fvalue += 3.0e0 * c3 * daux;

	if (c_1 != 0.0e0) {
		daux = cadd_1+_value;
		daux = daux * daux;
		fvalue -= c_1 / daux;
	};

	return fvalue;

};

// Author: Kharchenko S.A.
// Description: Output the nonlinear function description
// CFunc::operator<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CFunc &_func) { // Output the nonlinear function description

	_stream << " Function: C0 = " << _func.c0 << "; C1 = " << _func.c1;
	_stream << "; C2 = " << _func.c2 << "; C3 = " << _func.c3 << endl;
	_stream << "           C_1 = " << _func.c_1 << "; Cadd_1 = " << _func.cadd_1 << endl;

	return _stream;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CFuncCrValue::CFuncCrValue()
//========================================================================================
CFuncCrValue::CFuncCrValue () { // Memory allocation zero data constructor

//	const char *funcname = "CFuncCrValue_00";

	crvalue = 0.0e0;

};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CFuncCrValue::CFuncCrValue()
//========================================================================================
CFuncCrValue::CFuncCrValue (const CFuncCrValue &_func) { // Copy constructor

//	const char *funcname = "CFuncCrValue_copy";

	crvalue = _func.crvalue;
	func_undercr = _func.func_undercr;
	func_overcr = _func.func_overcr;

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CFuncCrValue::operator=()
//========================================================================================
CFuncCrValue &CFuncCrValue::operator= (const CFuncCrValue &_func) { // Equality operator

//	const char *funcname = "CFuncCrValue_=";

	crvalue = _func.crvalue;
	func_undercr = _func.func_undercr;
	func_overcr = _func.func_overcr;

	return *this;

};

// Author: Kharchenko S.A.
// Description: Compute function value
// CFuncCrValue::ComputeValue()
//========================================================================================
double CFuncCrValue::ComputeValue (double _value) { // Compute function value

//	const char *funcname = "ComputeValue";

	double fvalue = 0.0e0;

	if (_value <= crvalue) {
		fvalue = func_undercr.ComputeValue (_value);
	} else {
		fvalue = func_overcr.ComputeValue (_value);
	};

	return fvalue;

};

// Author: Kharchenko S.A.
// Description: Compute derivative value
// CFuncCrValue::ComputeDerivative()
//========================================================================================
double CFuncCrValue::ComputeDerivative (double _value) { // Compute derivative value

//	const char *funcname = "ComputeDerivative";

	double fvalue = 0.0e0;

	if (_value <= crvalue) {
		fvalue = func_undercr.ComputeDerivative (_value);
	} else {
		fvalue = func_overcr.ComputeDerivative (_value);
	};

	return fvalue;

};

// Author: Kharchenko S.A.
// Description: Output the nonlinear function description
// CFuncCrValue::operator<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CFuncCrValue &_func) { // Output the nonlinear function description

	_stream << " Function with critical value: CrValue = " << _func.crvalue << endl;
	_stream << " Under critical part:" << endl;
	_stream << _func.func_undercr << endl;
	_stream << " Over critical part:" << endl;
	_stream << _func.func_overcr << endl;

	return _stream;

};
