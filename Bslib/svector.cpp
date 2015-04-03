//------------------------------------------------------------------------------------------------
// File: svector.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "globals.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSVector: Zero length vector constructor
//========================================================================================
CSVector::CSVector () { // Zero length vector constructor

	const char *funcname = "CSVector_00";

	nv = 0;
	nrhs = 0;
	vect = new double [nv];
	if (!vect) MemoryFail (funcname);
	for (int i=0;i<nv;i++) {
		vect[i] = 0.0e0;
	};
	nparts = 1;
	blksv = new int [nparts+1];
	if (!blksv) MemoryFail (funcname);
	blksv[0] = 0;
	blksv[1] = nv;
	inda = new int [nparts];
	if (!inda) MemoryFail (funcname);
	inda[0] = 0;
};

// Author: Kharchenko S.A.
// CSVector: Constructor
//========================================================================================
CSVector::CSVector (int _n) { // Constructor

	const char *funcname = "CSVector_01";

	nv = _n;
	nrhs = 1;
	vect = new double [nv*nrhs];
	if (!vect) MemoryFail (funcname);
	for (int i=0;i<nv*nrhs;i++) {
		vect[i] = 0.0e0;
	};
	nparts = 1;
	blksv = new int [nparts+1];
	if (!blksv) MemoryFail (funcname);
	blksv[0] = 0;
	blksv[1] = nv;
	inda = new int [nparts];
	if (!inda) MemoryFail (funcname);
	inda[0] = 0;
};

// Author: Kharchenko S.A.
// CSVector: Constructor
//========================================================================================
CSVector::CSVector (int _n, int _nrhs) { // Constructor

	const char *funcname = "CSVector_02";

	nv = _n;
	nrhs = _nrhs;
	vect = new double [nv*nrhs];
	if (!vect) MemoryFail (funcname);
	for (int i=0;i<nv*nrhs;i++) {
		vect[i] = 0.0e0;
	};
	nparts = 1;
	blksv = new int [nparts+1];
	if (!blksv) MemoryFail (funcname);
	blksv[0] = 0;
	blksv[1] = nv;
	inda = new int [nparts];
	if (!inda) MemoryFail (funcname);
	inda[0] = 0;

};

// Author: Kharchenko S.A.
// CSVector: One more constructor
//========================================================================================
CSVector::CSVector (int _n, int _np, int _nrhs) { // One more constructor

	const char *funcname = "CSVector_03";

	nv = _n;
	nrhs = _nrhs;
	vect = new double [nv*nrhs];
	if (!vect) MemoryFail (funcname);
	int i;
	for (i=0;i<nv*nrhs;i++) {
		vect[i] = 0.0e0;
	};
	nparts = _np;
	blksv = new int [nparts+1];
	if (!blksv) MemoryFail (funcname);
	for (i=0;i<nparts;i++) {
		blksv[i] = 0;
	};
	blksv[nparts] = nv;
	inda = new int [nparts];
	if (!inda) MemoryFail (funcname);
	for (i=0;i<nparts;i++) {
		inda[i] = 0;
	};
};

// Author: Kharchenko S.A.
// CSVector: Constructor for parallel mode
//========================================================================================
CSVector::CSVector (const CGSMatrix &_gmtra, const CSVector &_x) { // Constructor for parallel mode

	const char *funcname = "CSVector_04";

// Count nloc

	int nloc=0;

	int ilist;
	for (ilist=0;ilist<_gmtra.nlist;ilist++) {
		int iblk = _gmtra.listb[ilist];
		int ni = _gmtra.bl2ndr[iblk+1]-_gmtra.bl2ndr[iblk];
		nloc += ni;
	};

// Allocate the data

	nv = nloc;
	nrhs = _x.nrhs;
	vect = new double [nv*nrhs];
	if (!vect) MemoryFail (funcname);
	nparts = _gmtra.nlist;
	blksv = new int [nparts+1];
	if (!blksv) MemoryFail (funcname);
	inda = new int [nparts];
	if (!inda) MemoryFail (funcname);

// Fill the data

	blksv[0] = 0;
	nloc = 0;

	for (ilist=0;ilist<_gmtra.nlist;ilist++) {
		int iblk = _gmtra.listb[ilist];
		int ni = _gmtra.bl2ndr[iblk+1]-_gmtra.bl2ndr[iblk];
		int ibs = _gmtra.bl2ndr[iblk];
		for (int irhs=0;irhs<nrhs;irhs++) {
			for (int kii=0;kii<ni;kii++) {
				vect[nv*irhs+nloc+kii] = _x.vect[_x.nv*irhs+ibs+kii];
			};
		};
		nloc += ni;
		blksv[ilist+1] = nloc;
		inda[ilist] = _gmtra.bl2ndr[iblk];
	};

};

// Author: Kharchenko S.A.
// CSVector: Destructor
//========================================================================================
CSVector::~CSVector () { // Destructor
//	static int icount = 0;
//	icount++;
//	std::cout << " CSVector destructor: " << icount << endl;
	nv = 0;
	nrhs = 0;
	nparts = 0;
	delete [] vect;
	delete [] blksv;
	delete [] inda;
//	std::cout << " On return from CSVector destructor " << endl;
};

// Author: Kharchenko S.A.
// CSVector: Copy constructor
//========================================================================================
CSVector::CSVector (const CSVector &_vv) { // Copy constructor

	const char *funcname = "CSVector_copy";

	nv = _vv.nv;
	nrhs = _vv.nrhs;
	vect = new double [nv*nrhs];
	if (!vect) MemoryFail (funcname);
	int i;
	for (i=0;i<nv*nrhs;i++) {
		vect[i] = _vv.vect[i];
	};
	nparts = _vv.nparts;
	blksv = new int [nparts+1];
	if (!blksv) MemoryFail (funcname);
	for (i=0;i<=nparts;i++) {
		blksv[i] = _vv.blksv[i];
	};
	inda = new int [nparts];
	if (!inda) MemoryFail (funcname);
	for (i=0;i<nparts;i++) {
		inda[i] = _vv.inda[i];
	};
};

// Author: Kharchenko S.A.
// CSVector: Equality operator
//========================================================================================
CSVector &CSVector::operator= (const CSVector &_v2) { // Equality operator

	const char *funcname = "CSVector_=";

	if (nv != _v2.nv || nrhs != _v2.nrhs || nparts != _v2.nparts) {
		delete [] vect;
		delete [] blksv;
		delete [] inda;
		nv = _v2.nv;
		nrhs = _v2.nrhs;
		vect = new double [nv*nrhs];
		if (!vect) MemoryFail (funcname);
		nparts = _v2.nparts;
		blksv = new int [nparts+1];
		if (!blksv) MemoryFail (funcname);
		inda = new int [nparts];
		if (!inda) MemoryFail (funcname);
	};
	int i;
	for (i=0;i<nv*nrhs;i++) {
		vect[i] = _v2.vect[i];
	};
	for (i=0;i<nparts+1;i++) {
		blksv[i] = _v2.blksv[i];
	};
	for (i=0;i<nparts;i++) {
		inda[i] = _v2.inda[i];
	};
	return *this;
};

// Author: Kharchenko S.A.
// CSVector: Set a float value for all components
//========================================================================================
void CSVector::SetSVect (double _value) { // Set a float value for all components
	for (int i=0;i<nv*nrhs;i++) {
		vect[i] = _value;
	};
};

// Author: Kharchenko S.A.
// CSVector: Set a float value for a component
//========================================================================================
void CSVector::SetComponentSVect (int _index, double _value) { // Set a float value for component
	vect[_index] = _value;
};

// Author: Kharchenko S.A.
// CSVector: Add two vectors
//========================================================================================
CSVector CSVector::operator+ (const CSVector &_v2) { // Add two vectors
	CSVector temp(nv,nparts,nrhs);
	int i;
	for (i=0;i<nv*nrhs;i++) {
		temp.vect[i] = vect[i] + _v2.vect[i];
	};
	for (i=0;i<=nparts;i++) {
		temp.blksv[i] = blksv[i];
	};
	for (i=0;i<nparts;i++) {
		temp.inda[i] = inda[i];
	};
	return temp;
};

// Author: Kharchenko S.A.
// CSVector: Substract two vectors
//========================================================================================
CSVector CSVector::operator- (const CSVector &_v2) { // Substract two vectors
	CSVector temp(nv,nparts,nrhs);
	int i;
	for (i=0;i<nv*nrhs;i++) {
		temp.vect[i] = vect[i] - _v2.vect[i];
	};
	for (i=0;i<=nparts;i++) {
		temp.blksv[i] = blksv[i];
	};
	for (i=0;i<nparts;i++) {
		temp.inda[i] = inda[i];
	};
	return temp;
};

// Author: Kharchenko S.A.
// CSVector: Multiply vector by a scalar
//========================================================================================
CSVector CSVector::operator* (double _alpha) { // Multiply vector by a scalar
	CSVector temp(nv,nparts,nrhs);
	int i;
	for (i=0;i<nv*nrhs;i++) {
		temp.vect[i] = vect[i] * _alpha;
	};
	for (i=0;i<=nparts;i++) {
		temp.blksv[i] = blksv[i];
	};
	for (i=0;i<nparts;i++) {
		temp.inda[i] = inda[i];
	};
	return temp;
};

// Author: Kharchenko S.A.
// CSVector: Compute y += alpha*x
//========================================================================================
void CSVector::DaxpySVect (double _alpha, const CSVector &_x) { // Compute y += alpha*x
	for (int i=0;i<nv*nrhs;i++) vect[i] += _alpha * _x.vect[i];
};

// Author: Kharchenko S.A.
// CSVector: Compute y += alpha*x*a
//========================================================================================
void CSVector::DaxpySVectBlk (double _alpha, const double *_a, const CSVector &_x) { // Compute y += alpha*x*a

	int i, j, k;
	double aux;

	for (i=0;i<nrhs;i++) {
		for (j=0;j<nv;j++) {
			aux = 0.0e0;
			for (k=0;k<nrhs;k++) {
				aux += _x.vect[k*nv+j] * _a[i*nrhs+k];
			};
			vect[i*nv+j] += _alpha * aux;
		};
	};

};

// Author: Kharchenko S.A.
// CSVector: Compute x *= alpha
//========================================================================================
void CSVector::ScaleSVect (double _alpha) { // Compute x *= alpha
	for (int i=0;i<nv*nrhs;i++) vect[i] *= _alpha;
};

// Author: Kharchenko S.A.
// CSVector: Component product
//========================================================================================
double CSVector::CompProdSVect (int _i, const CSVector &_v2) { // Component product
	double sc;
	sc = vect[_i] * _v2.vect[_i];
	return sc;
};	

// Author: Kharchenko S.A.
// CSVector: Scalar product
//========================================================================================
double CSVector::ScProdSVect (const CSVector &_v2) { // Scalar product
	double sc;
	sc = 0.0e0;
	for (int i=0;i<nv;i++) {
		sc = sc + vect[i] * _v2.vect[i];
	};
	return sc;
};

// Author: Kharchenko S.A.
// CSVector: Maximum of component products
//========================================================================================
double CSVector::MaxProdSVect (const CSVector &_v2) const { // Maximum of component products
	double sc, aux;
	sc = 0.0e0;
	for (int i=0;i<nv;i++) {
		aux = vect[i] * _v2.vect[i];
		if (aux < 0.0e0) aux = -aux;
		if (sc < aux) sc = aux;
	};
	return sc;
};

// Author: Kharchenko S.A.
// CSVector: Maximum of chosen component products
//========================================================================================
double CSVector::MaxProdSVect (const CSVector &_v2, const bool *_barr) const { // Maximum of chosen component products
	double sc, aux;
	sc = 0.0e0;
	for (int i=0;i<nv;i++) {
		if (_barr[i]) {
			aux = vect[i] * _v2.vect[i];
			if (aux < 0.0e0) aux = -aux;
			if (sc < aux) sc = aux;
		};
	};
	return sc;
};

// Author: Kharchenko S.A.
// CSVector: Reorder the vector
//========================================================================================
void CSVector::OrdVect (char _ordtyp, const int *_order, const CSVector &_x) { // Reorder the vector
	int i;
	if (_ordtyp == 'D') {
		for (i=0;i<nv;i++) vect[_order[i]] = _x.vect[i];
	} else {
		for (i=0;i<nv;i++) vect[i] = _x.vect[_order[i]];
	};
};

// Author: Kharchenko S.A.
// CSVector: Order the vector according to the supernodes partitioning
//========================================================================================
void CSVector::OrdVect (char _ortype, // Order the vector according to the supernodes partitioning
						int _nsupc, int *_sprndc, int *_sprndcnew, int *_order) {

	const char *funcname = "OrdVect";

	double *vloc;

	vloc = new double [nv];
	if (!vloc) MemoryFail (funcname);

	int isupo, isupn, blai, ibsn, ibso, kii, irhs;

	for (irhs=0;irhs<nrhs;irhs++) {

		for (int i=0;i<nv;i++) vloc[i] = vect[i+irhs*nv];

		if (_ortype == 'D') {

			for (isupo=0;isupo<_nsupc;isupo++) {
				isupn = _order[isupo];
				ibso = _sprndc[isupo];
				blai = _sprndcnew[isupn+1]-_sprndcnew[isupn];
				ibsn = _sprndcnew[isupn];
				for (kii=0;kii<blai;kii++) vect[irhs*nv+ibsn+kii] = vloc[ibso+kii];
			};
	
		} else {

			for (isupo=0;isupo<_nsupc;isupo++) {
				isupn = _order[isupo];
				ibso = _sprndc[isupo];
				blai = _sprndcnew[isupn+1]-_sprndcnew[isupn];
				ibsn = _sprndcnew[isupn];
				for (kii=0;kii<blai;kii++) vect[irhs*nv+ibso+kii] = vloc[ibsn+kii];
			};

		};

	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSVector: Read vector from disk
//========================================================================================
CSVector CSVector::ReadVect (istream &_fin, int _nv) { // Read vector from disk

	CSVector temp(_nv);

	int i;
	char buffer[20]="";
 
	_fin.getline(buffer,10);
	cout << buffer << endl;
	_fin.getline(buffer,10);
	cout << buffer << endl;

	for (i=0; i<_nv; i++)
	{
		_fin >> temp.vect[i];
	};

	return temp;
};

// Author: Kharchenko S.A.
// CSVector: Read vector from disk in binary form
//========================================================================================
CSVector CSVector::ReadVectBin (istream &_stream) { // Read vector from disk in binary form

	int npartsl, nrhsl, nvl;

	ReadArrAccV (_stream, 1, &npartsl);
	ReadArrAccV (_stream, 1, &nrhsl);
	ReadArrAccV (_stream, 1, &nvl);

	cout << " CSVector:" << endl;
	cout << " Nparts = " << npartsl << " Nrhs = " << nrhsl << " Nv = " << nvl << endl;

	CSVector temp (nvl,npartsl,nrhsl);

	ReadArrAccV (_stream, nvl*nrhsl,  temp.vect);
	ReadArrAccV (_stream, npartsl+1, temp.blksv);
	ReadArrAccV (_stream, npartsl,    temp.inda);

	return temp;

};

// Author: Kharchenko S.A.
// CSVector: Output vector
//========================================================================================
ostream &operator<< (ostream &_stream, const CSVector &_v) { // Output vector

	_stream << " CSVector:" << endl;
	_stream << " Nv = " << _v.nv << " Nrhs = " << _v.nrhs << " Nparts = " << _v.nparts << endl;

	OutArr (_stream, " V =", _v.nv*_v.nrhs, _v.vect);
	OutArr (_stream, " BlksV =", _v.nparts+1, _v.blksv);
	OutArr (_stream, " IndA  =", _v.nparts,   _v.inda);

	return _stream;
};

// Author: Kharchenko S.A.
// CSVector: Output vector to the disk
//========================================================================================
void CSVector::WriteVect (ostream &_stream) { // Output vector to the disk

	int i;
	for (i=0;i<nv;i++) _stream << "   " << setw(23) << setprecision(16) << vect[i] << endl;

};

// Author: Kharchenko S.A.
// CSVector: Output vector to the disk in a binary form
//========================================================================================
void CSVector::WriteVectBin (ostream &_stream) { // Output vector to the disk in a binary form

	OutArrAccV (_stream, " Nparts =", 1,     &nparts);
	OutArrAccV (_stream, " Nrhs =", 1,         &nrhs);
	OutArrAccV (_stream, " Nv =", 1,             &nv);
	OutArrAccV (_stream, " V =", nv*nrhs,       vect);
	OutArrAccV (_stream, " BlksV =", nparts+1, blksv);
	OutArrAccV (_stream, " IndA  =", nparts,    inda);

};

// Author: Kharchenko S.A.
// CSVector: Output vector to the disk in a binary form
//========================================================================================
void CSVector::WriteVectBin (FILE * _file, int &_offset) { // Output vector to the disk in a binary form

	FPut (_file,1,&nparts,_offset);
	_offset += 1;
	FPut (_file,1,&nrhs,_offset);
	_offset += 1;
	FPut (_file,1,&nv,_offset);
	_offset += 1;
	FPut (_file,nv*nrhs,vect,_offset);
	_offset += nv*nrhs;
	FPut (_file,nparts+1,blksv,_offset);
	_offset += nparts+1;
	FPut (_file,nparts,inda,_offset);
	_offset += nparts;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
