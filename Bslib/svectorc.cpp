//------------------------------------------------------------------------------------------------
// File: svectorc.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include "globals.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CSVectorC: Reorder vector components
//========================================================================================
void CSVectorC::OrderVector (int _ordtype, int _n, int *_order, dcmplx *_x, dcmplx *_ax) { // Reorder vector components

	int kii, kjj;

	if (_ordtype == 1) {
		for (kii=0;kii<_n;kii++) {
			kjj = _order[kii];
			_ax[kjj] = _x[kii];
		};
	} else {
		for (kii=0;kii<_n;kii++) {
			kjj = _order[kii];
			_ax[kii] = _x[kjj];
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Multiply by the rectangular block and add
//========================================================================================
void CSVectorC::MvmBlock  (int _m, int _n, dcmplx *_ablk, dcmplx *_x, dcmplx *_ax) { // Multiply by the rectangular block and add

	int kii, kjj;

	dcmplx *paxloc;
	dcmplx *paloc = _ablk;
	dcmplx caux;

	for (kii=0;kii<_n;kii++) {
		caux = _x[kii];
		paxloc = _ax;
		for (kjj=0;kjj<_m;kjj++) {
			*paxloc += *paloc * caux;
//			*paxloc.x += paloc->x * caux.x - paloc->y * caux.y;
//			*paxloc.y += paloc->x * caux.y + paloc->y * caux.x;
			paloc++;
			paxloc++;
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Multiply by the rectangular block and add
//========================================================================================
void CSVectorC::MvmBlockH (int _m, int _n, dcmplx *_ablk, dcmplx *_x, dcmplx *_ax) { // Multiply by the rectangular block and add

	int kii, kjj;

	dcmplx *pxloc;
	dcmplx *paxloc;
	dcmplx *paloc = _ablk;

	paxloc = _ax;
	for (kii=0;kii<_n;kii++) {
		pxloc = _x;
		for (kjj=0;kjj<_m;kjj++) {
			*paxloc += conj(*paloc) * *pxloc;
//			paxloc->x += paloc->x * *pxloc.x + paloc->y * *pxloc.y;
//			paxloc->y += paloc->x * *pxloc.y - paloc->y * *pxloc.x;
			paloc++;
			pxloc++;
		};
		paxloc++;
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Zero length vector constructor
//========================================================================================
CSVectorC::CSVectorC () { // Zero length vector constructor

	const char *funcname = "CSVectorC_00";

	dcmplx czero (0.0e0,0.0e0);

	nv = 0;
	nrhs = 0;
	vect = new dcmplx [nv+1];
	if (!vect) MemoryFail (funcname);
	for (int i=0;i<nv;i++) {
		vect[i] = czero;
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
// CSVectorC: Constructor
//========================================================================================
CSVectorC::CSVectorC (int _n) { // Constructor

	const char *funcname = "CSVectorC_01";

	dcmplx czero (0.0e0,0.0e0);

	nv = _n;
	nrhs = 1;
	vect = new dcmplx [nv*nrhs+1];
	if (!vect) MemoryFail (funcname);
	for (int i=0;i<nv*nrhs;i++) {
		vect[i] = czero;
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
// CSVectorC: Constructor
//========================================================================================
CSVectorC::CSVectorC (int _n, int _nrhs) { // Constructor

	const char *funcname = "CSVectorC_02";

	dcmplx czero (0.0e0,0.0e0);

	nv = _n;
	nrhs = _nrhs;
	vect = new dcmplx [nv*nrhs+1];
	if (!vect) MemoryFail (funcname);
	for (int i=0;i<nv*nrhs;i++) {
		vect[i] = czero;
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
// CSVectorC: One more constructor
//========================================================================================
CSVectorC::CSVectorC (int _n, int _np, int _nrhs) { // One more constructor

	const char *funcname = "CSVectorC_03";

	dcmplx czero (0.0e0,0.0e0);

	nv = _n;
	nrhs = _nrhs;
	vect = new dcmplx [nv*nrhs+1];
	if (!vect) MemoryFail (funcname);
	int i;
	for (i=0;i<nv*nrhs;i++) {
		vect[i] = czero;
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
// CSVectorC: Constructor for parallel mode
//========================================================================================
CSVectorC::CSVectorC (const CGSMatrixCS &_gmtra, const CSVectorC &_x) { // Constructor for parallel mode

	const char *funcname = "CSVectorC_04";

	dcmplx czero (0.0e0,0.0e0);

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
	vect = new dcmplx [nv*nrhs+1];
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
// CSVectorC: Destructor
//========================================================================================
CSVectorC::~CSVectorC () { // Destructor
//	static int icount = 0;
//	icount++;
//	std::cout << " CSVectorC destructor: " << std::endl;
	nv = 0;
	nrhs = 0;
	nparts = 0;
	delete [] vect;
	delete [] blksv;
	delete [] inda;
//	std::cout << " On return from CSVectorC destructor " << std::endl;
};

// Author: Kharchenko S.A.
// CSVectorC: Copy constructor
//========================================================================================
CSVectorC::CSVectorC (const CSVectorC &_vv) { // Copy constructor

	const char *funcname = "CSVectorC_copy";

	nv = _vv.nv;
	nrhs = _vv.nrhs;
	vect = new dcmplx [nv*nrhs+1];
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
// CSVectorC: Equality operator
//========================================================================================
CSVectorC &CSVectorC::operator= (const CSVectorC &_v2) { // Equality operator

	const char *funcname = "CSVectorC_=";

	if (nv != _v2.nv || nrhs != _v2.nrhs || nparts != _v2.nparts) {
		delete [] vect;
		delete [] blksv;
		delete [] inda;
		nv = _v2.nv;
		nrhs = _v2.nrhs;
		vect = new dcmplx [nv*nrhs+1];
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
// CSVectorC: Set a complex value for all components
//========================================================================================
void CSVectorC::SetSVect (dcmplx _value) { // Set a complex value for all components
	for (int i=0;i<nv*nrhs;i++) {
		vect[i] = _value;
	};
};

// Author: Kharchenko S.A.
// CSVectorC: Set a complex value for a component
//========================================================================================
void CSVectorC::SetComponentSVect (int _index, dcmplx _value) { // Set a complex value for component
	vect[_index] = _value;
};

// Author: Kharchenko S.A.
// CSVectorC: Add two vectors
//========================================================================================
CSVectorC CSVectorC::operator+ (const CSVectorC &_v2) { // Add two vectors
	CSVectorC temp(nv,nparts,nrhs);
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
// CSVectorC: Substract two vectors
//========================================================================================
CSVectorC CSVectorC::operator- (const CSVectorC &_v2) { // Substract two vectors
	CSVectorC temp(nv,nparts,nrhs);
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
// CSVectorC: Multiply vector by a scalar
//========================================================================================
CSVectorC CSVectorC::operator* (dcmplx _alpha) { // Multiply vector by a scalar
	CSVectorC temp(nv,nparts,nrhs);
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
// CSVectorC: Scalar product
//========================================================================================
dcmplx CSVectorC::ScProdSVect (const CSVectorC &_v2) { // Scalar product
	dcmplx sc (0.0e0,0.0e0);
	for (int i=0;i<nv;i++) {
		sc = sc + vect[i] * conj(_v2.vect[i]);
	};
	return sc;
};

// Author: Kharchenko S.A.
// CSVectorC: Compute y += alpha*x
//========================================================================================
void CSVectorC::DaxpySVect (dcmplx _alpha, const CSVectorC &_x) { // Compute y += alpha*x
//	for (int i=0;i<nv*nrhs;i++) vect[i] += _alpha * _x.vect[i];
	for (int i=0;i<nv*nrhs;i++) vect[i].PlEq (_alpha * _x.vect[i]);
};

// Author: Kharchenko S.A.
// CSVectorC: Compute y += alpha*x*a
//========================================================================================
void CSVectorC::DaxpySVectBlk (dcmplx _alpha, const dcmplx *_a, const CSVectorC &_x) { // Compute y += alpha*x*a

	int i, j, k;
	dcmplx aux;

	for (i=0;i<nrhs;i++) {
		for (j=0;j<nv;j++) {
			aux.real(0.0e0);
			aux.imag(0.0e0);
			for (k=0;k<nrhs;k++) {
//				aux += _x.vect[k*nv+j] * _a[i*nrhs+k];
				aux.PlEq (_x.vect[k*nv+j] * _a[i*nrhs+k]);
			};
//			vect[i*nv+j] += _alpha * aux;
			vect[i*nv+j].PlEq (_alpha * aux);
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Compute x *= alpha
//========================================================================================
void CSVectorC::ScaleSVect (dcmplx _alpha) { // Compute x *= alpha
	for (int i=0;i<nv*nrhs;i++) vect[i] *= _alpha;
};

// Author: Kharchenko S.A.
// CSVectorC: Reorder the vector
//========================================================================================
void CSVectorC::OrdVect (char _ordtyp, const int *_order, const CSVectorC &_x) { // Reorder the vector
	int i;
	if (_ordtyp == 'D') {
		for (i=0;i<nv;i++) vect[_order[i]] = _x.vect[i];
	} else {
		for (i=0;i<nv;i++) vect[i] = _x.vect[_order[i]];
	};
};

// Author: Kharchenko S.A.
// CSVectorC: Order the vector according to the supernodes partitioning
//========================================================================================
void CSVectorC::OrdVect (char _ortype, // Order the vector according to the supernodes partitioning
							int _nsupc, int *_sprndc, int *_sprndcnew, int *_order) {

	const char *funcname = "OrdVect";

	dcmplx *vloc;

	vloc = new dcmplx [nv];
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
// CSVectorC: Perform point scaling of the vector according to the blocks partitioning
//========================================================================================
void CSVectorC::ScaleVect (char _sctype, // Perform point scaling of the vector according to the blocks partitioning
							const CSVectorC &_sclpt) {

//	const char *funcname = "ScaleVect";

// Perform scaling operations

	int i, j, irhs, kii, iblkscl, ibeg, ibegn, ni;

	iblkscl = 0;

	for (i=0;i<nparts;i++) {
		ibeg = inda[i];
		int icheck = 0;
		for (j=iblkscl;j<_sclpt.nparts;j++) {
			if (_sclpt.inda[j] == ibeg) {
				icheck = 1;
				iblkscl = j;
				break;
			};
		};
		if (icheck == 0) throw " ScaleVect: required block is not found ";
		ibeg = blksv[i];
		ibegn = _sclpt.blksv[iblkscl];
		ni = blksv[i+1]-blksv[i];
		if (_sctype == 'D') {
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<ni;kii++) {
					vect[irhs*nv+ibeg+kii] *= _sclpt.vect[ibegn+kii];
				};
			};
		} else {
			for (irhs=0;irhs<nrhs;irhs++) {
				for (kii=0;kii<ni;kii++) {
					vect[irhs*nv+ibeg+kii] /= _sclpt.vect[ibegn+kii];
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSVectorC: Perform block scaling of the vector according to the blocks partitioning
//========================================================================================
void CSVectorC::BlockScaleVect (char _transp, const CGSMatrixCS &_mtra, // Perform block scaling of the vector according to the blocks partitioning
								const CSVectorC &_sclblk) {

	const char *funcname = "BlockScaleVect";

	dcmplx czero (0.0e0,0.0e0);

// Allocate work data

	dcmplx *vloc;

	int iblk = _mtra.listb[0];

	vloc = new dcmplx [_mtra.mtrarr[iblk].blamx];
	if (!vloc) MemoryFail (funcname);

// Perform scaling operations

	int i, j, irhs, kii, kjj, iblkscl, ibeg, ni, ibs, isup;

	dcmplx caux;

	iblk = 0;
	iblkscl = 0;

	for (i=0;i<nparts;i++) {
		ibeg = inda[i];
		int icheck = 0;
		for (j=iblk;j<_mtra.nblksc;j++) {
			if (_mtra.bl2ndc[j] == ibeg) {
				icheck = 1;
				iblk = j;
				break;
			};
		};
		if (icheck == 0) throw " BlockScaleVect: required block is not found ";
		icheck = 0;
		for (j=iblkscl;j<_sclblk.nparts;j++) {
			if (_sclblk.inda[j] == ibeg) {
				icheck = 1;
				iblkscl = j;
				break;
			};
		};
		if (icheck == 0) throw " BlockScaleVect: required block is not found ";
		ibeg = blksv[i];
		ibs = _sclblk.blksv[iblkscl];
		for (isup=_mtra.blksc[iblk];isup<_mtra.blksc[iblk+1];isup++) {
			ni = _mtra.sprndc[isup+1]-_mtra.sprndc[isup];
			for (irhs=0;irhs<nrhs;irhs++) {
				if (_transp == 'N') {
					for (kii=0;kii<ni;kii++) {
						caux = czero;
						for (kjj=0;kjj<ni;kjj++) {
							caux += _sclblk.vect[ibs+kjj*ni+kii] * vect[irhs*nv+ibeg+kjj];
						};
						vloc[kii] = caux;
					};
				} else {
					for (kii=0;kii<ni;kii++) {
						caux = czero;
						for (kjj=0;kjj<ni;kjj++) {
							caux += _sclblk.vect[ibs+kii*ni+kjj] * vect[irhs*nv+ibeg+kjj];
						};
						vloc[kii] = caux;
					};
				};
				for (kii=0;kii<ni;kii++) vect[irhs*nv+ibeg+kii] = vloc[kii];
			};
			ibs += ni*ni;
			ibeg += ni;
		};
	};

// Free work array

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSVectorC: Output vector
//========================================================================================
ostream &operator<< (ostream &_stream, const CSVectorC &_v) { // Output vector

	_stream << " CSVectorC:" << endl;
	_stream << " Nv = " << _v.nv << " Nrhs = " << _v.nrhs << " Nparts = " << _v.nparts << endl;

	OutArr (_stream, " V =", _v.nv*_v.nrhs, _v.vect);
	OutArr (_stream, " BlksV =", _v.nparts+1, _v.blksv);
	OutArr (_stream, " IndA  =", _v.nparts,   _v.inda);

	return _stream;
	
};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
