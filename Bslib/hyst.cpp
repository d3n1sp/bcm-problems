//------------------------------------------------------------------------------------------------
// File: hyst.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "hyst.h"
#include "globals.h"

using namespace std;

// Author: Kharchenko S.A.
// CHyst: Memory allocation zero data constructor
//========================================================================================
CHyst::CHyst () { // Memory allocation zero data constructor

	const char *funcname = "CHyst_00";

	iflag = 0;
	nhyst = 32;
	icenter = 16;
	elemmax = 0.0e0;
	elemmin = 0.0e0;

	bounds = new double [nhyst];
	if (!bounds) MemoryFail (funcname);
	counts = new int [nhyst];
	if (!counts) MemoryFail (funcname);

	bounds[0] = 1.0e-16;

	int i;

	for (i=1;i<nhyst;i++) bounds[i] = bounds[i-1] * 10.0e0;

	for (i=0;i<nhyst;i++) counts[i] = 0;

};

// Author: Kharchenko S.A.
// CHyst: Update the hystogram data
//========================================================================================
void CHyst::UpdateHystogram (int _n, double *_values) { // Update the hystogram data

//	const char *funcname = "UpdateHystogram";

// Init if necessary maximal and minimal elements

	if (_n == 0) return;

	double aux = _values[0];
	if (aux < 0.0e0) aux = -aux;

	if (iflag == 0) {
		elemmin = aux;
		elemmax = aux;
	};

// Update hystogram

	for (int i=0;i<_n;i++) {
		aux = _values[i];
		if (aux < 0.0e0) aux = -aux;
		if (aux > elemmax) elemmax = aux;
		if (aux < elemmin) elemmin = aux;

		int ipos = icenter;

		while (ipos >= 0) {
			if (aux >= bounds[ipos]) {
				if (aux < bounds[ipos+1]) {
					counts[ipos]++;
					ipos = -1;
				} else {
					if (ipos < nhyst-1) {
						ipos++;
					} else {
						counts[nhyst-1]++;
						ipos = -1;
					};
				};
			} else {
				if (ipos>0) {
					ipos--;
				} else {
					counts[0]++;
					ipos = -1;
				};
			};
		};
	};

	iflag = 1;

};

// Author: Kharchenko S.A.
// CHyst: Perform parallel update of the hystogram data
//========================================================================================
void CHyst::MakeGlobalHystogram (const CMPIComm &_comm) { // Perform parallel update of the hystogram data

	const char *funcname = "MakeGlobalHystogram";

// Update maximal and minimal elements

	int i=0;

	int *iarr;

	iarr = new int [nhyst];
	if (!iarr) MemoryFail (funcname);

	double aux;

	int nproc = _comm.GetNproc ();

	if (nproc == 1) goto Exit;

	CMPIExchange::ExchangeArrayMPI (_comm,
												DOUBLEVALUE, MAXIMUM,
												1, &elemmax, &aux);

	elemmax = aux;

	CMPIExchange::ExchangeArrayMPI (_comm,
												DOUBLEVALUE, MINIMUM,
												1, &elemmin, &aux);

	elemmin = aux;

	CMPIExchange::ExchangeArrayMPI (_comm,
												INTEGERVALUE, ADD,
												nhyst, counts, iarr);

	for (i=0;i<nhyst;i++) counts[i] = iarr[i];

Exit:;

	delete [] iarr;

};

// Author: Kharchenko S.A.
// CHyst: Perform serial update of the hystogram data
//========================================================================================
void CHyst::MakeGlobalHystogram () { // Perform serial update of the hystogram data

//	const char *funcname = "MakeGlobalHystogram";

};

// Author: Kharchenko S.A.
// Hyst: Output hystogram
//========================================================================================
ostream &operator<< (ostream &_stream, const CHyst &_hyst) { // Output hystogram

	int indmin = _hyst.nhyst-1, indmax=0;
	int nloc = 0;

	int i;

	for (i=0;i<_hyst.nhyst;i++) {
		if (_hyst.counts[i] != 0) indmax = i;
		nloc += _hyst.counts[i];
	};

	_stream << " Hystogram: N = " << nloc << "; SMin = " << _hyst.elemmin << "; SMax = " << _hyst.elemmax << ";" << endl;

	for (i=_hyst.nhyst-1;i>=0;i--) {
		if (_hyst.counts[i] != 0) indmin = i;
	};

	OutArr (_stream, " Hyst =", indmax-indmin+1, _hyst.counts+indmin);

	return _stream;

};

