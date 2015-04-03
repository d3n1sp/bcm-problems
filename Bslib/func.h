// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "globals.h"

// Func.h: Description of the nonlinear function
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Func
#define __Func

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CFunc
{
// Data
public:
	double c0;     // c0  specifies the coefficient for degree zero of a value in additive nonlinear function representation
	double c1;     // c1  specifies the coefficient for degree one  of a value in additive nonlinear function representation
	double c2;     // c2  specifies the coefficient for degree two  of a value in additive nonlinear function representation
	double c3;     // c3  specifies the coefficient for degree tree of a value in additive nonlinear function representation
	double c_1;    // c_1 specifies the coefficient for degree minus one of a value in additive nonlinear function representation
	double cadd_1; // cadd_1 specifies the addition coefficient for degree minus one of a value in additive nonlinear function representation
public:
// Functions
// Constructors and destructor
	CFunc (); // Memory allocation zero data constructor
	~CFunc () { // Destructor
//		cout << " On entry to CFunc destructor " << endl;
//		cout << " On return from CFunc destructor " << endl;
	};
	CFunc (const CFunc &_func); // Copy constructor
// Operator functions
	CFunc &operator= (const CFunc &_func); // Equality operator
// Compute values and derivatives
	double ComputeValue (double _value); // Compute function value
	double ComputeDerivative (double _value); // Compute derivative value
// Input/Output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CFunc &_func); // Output function description
// Friend classes
};

class CFuncCrValue
{
// Data
public:
	double crvalue;     // crvalue contains critical value of the function
	CFunc func_undercr; // func_undercr contains under critical part of the function
	CFunc func_overcr;  // func_overcr  contains over  critical part of the function
public:
// Functions
// Constructors and destructor
	CFuncCrValue (); // Memory allocation zero data constructor
	~CFuncCrValue () { // Destructor
//		cout << " On entry to CFuncCrValue destructor " << endl;
//		cout << " On return from CFuncCrValue destructor " << endl;
	};
	CFuncCrValue (const CFuncCrValue &_func); // Copy constructor
// Operator functions
	CFuncCrValue &operator= (const CFuncCrValue &_func); // Equality operator
// Compute values and derivatives
	double ComputeValue (double _value); // Compute function value
	double ComputeDerivative (double _value); // Compute derivative value
	static double FuncValue (void *_obj, int _itype, int _iparam, double _value) {
		double fvalue;
		if (_itype == -1) {
			fvalue = ((CFuncCrValue *)_obj)->crvalue;
		} else if (_itype == 0) {
			fvalue = ((CFuncCrValue *)_obj)->ComputeValue (_value);
		} else if (_itype == 1) {
			fvalue = ((CFuncCrValue *)_obj)->ComputeDerivative (_value);
		};
		return fvalue;
	};
// Input/Output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CFuncCrValue &_func); // Output function description
// Friend classes
};

#endif

