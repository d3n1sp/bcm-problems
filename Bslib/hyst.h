//------------------------------------------------------------------------------------------------
// File: hyst.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
// Preliminary declarations

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "tree.h"
#include "globals.h"
#include "ExchangeMPI.h"

// hyst.h: Description of the hystogram data
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Hyst
#define __Hyst

// Preliminary declarations

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

class CHyst
{
protected:
// Data
	int iflag;         // iflag specifies that the hystogram was initialized or not by some data
	int nhyst;         // nhyst is the size of the hystogram array
	int icenter;       // icenter is the pointer to the central interval of the hystogram
	double elemmax;    // elemmax spacifies the maximal element in absolute value
	double elemmin;    // elemmin spacifies the minimal element in absolute value
	double *bounds;    // bounds[nhyst] describes bounds for each hystogram interval
	int *counts;       // counts[nhyst] describes the number of elements in each interval bound
public:
// Functions
// Constructors and destructor
	CHyst (); // Memory allocation zero data constructor
	~CHyst () { // Destructor
//		cout << " On entry to CHyst destructor " << endl;
		delete [] bounds;
		delete [] counts;
//		cout << " On return from CHyst destructor " << endl;
	};
	CHyst (const CHyst &_hyst) {throw " Copy constructor for Hyst class called";}; // Copy constructor
// Operator functions
	CHyst &operator= (const CHyst &_hyst) {throw " Equality operator for Hyst class called";}; // Equality operator
// Update hystogram functions
	void UpdateHystogram (int _n, double *_values); // Update the hystogram data
	void MakeGlobalHystogram (const CMPIComm &_comm); // Perform parallel update of the hystogram data
	void MakeGlobalHystogram (); // Perform serial update of the hystogram data
// Input/Output functions
	friend std::ostream &operator<< (std::ostream &_stream, const CHyst &_hyst); // Output hystogram
// Friend classes
};

#endif
