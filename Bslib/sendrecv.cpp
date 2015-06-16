//------------------------------------------------------------------------------------------------
// File: sendrecv.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

//#include "mpi.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "globals.h"
#include "mvm.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "fct.h"
#include "fctdiag.h"
#include "qrd.h"
#include "ExchangeMPI.h"

using namespace std;

// Author: Kharchenko S.A.
// CFctDiagR: Send scaling data
//========================================================================================
void CFctDiagR::SendScaling (bool _is2index, const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Send scaling data
								const CGSMatrix &_gmtra) {

	const char *funcname = "SendScaling";

	double *sendarr;

	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	int nztot = 0;

	int ilist, nlocblk, isupbeg, isupend;
	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		int nzloc = diasize[iblk];

		nztot += nlocblk+4*nzloc;

	};

	sendarr = new double [nztot];
	if (!sendarr) MemoryFail (funcname);

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc  = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		int kii;
		for (kii=0;kii<nlocblk;kii++) sendarr[nztot+kii] = sclptloc[kii];
		nztot += nlocblk;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = scllloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = scluloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = sclinvlloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = sclinvuloc[kii];
		nztot += nzloc;

	};

	CMPIComm commloc = _comm;

	int msgtag = _msgtag;

	if (nztot > 0) {
		CMPIExchange::Send (commloc, _iproc, msgtag,
									nztot*sizeof(double), (char *)sendarr);
	};

	delete [] sendarr;

};

// Author: Kharchenko S.A.
// CFctDiagR: Receive scaling data
//========================================================================================
void CFctDiagR::ReceiveScaling (bool _is2index, const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Receive scaling data
									const CGSMatrix &_gmtra) {

	const char *funcname = "ReceiveScaling";

//	double dzero = 0.0e0;

	double *recvarr;
	double *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	int nztot = 0;

	int ilist, nlocblk, isupbeg, isupend;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		nztot += nlocblk+4*nzloc;

	};

	recvarr = new double [nztot];
	if (!recvarr) MemoryFail (funcname);

	CMPIComm commloc = _comm;
	CMPIStatus status;

	int msgtag = _msgtag;

	if (nztot > 0) {
		CMPIExchange::Recv (commloc, _iproc, msgtag,
								nztot*sizeof(double), (char *)recvarr, status);
	};

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		if (_is2index) {
			nlocblk = isupend+1-isupbeg;
		} else {
			nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];
		};

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		int kii;
		for (kii=0;kii<nlocblk;kii++) sclptloc[kii] = recvarr[nztot+kii];
		nztot += nlocblk;

		for (kii=0;kii<nzloc;kii++) scllloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) scluloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sclinvlloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sclinvuloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

	};

	delete [] recvarr;

};

// Author: Kharchenko S.A.
// CFctDiagR: Send diagg array
//========================================================================================
void CFctDiagR::SendDiag (const CMPIComm &_comm, int _nlistloc, int* _listloc, int _iproc, int _msgtag, // Send diagg array
							const CGSMatrix &_gmtra) {

	const char *funcname = "SendDiag";

// Read the data from disk if necessary

	if (nfiles > 0) {

		for (int ilist=0;ilist<_nlistloc;ilist++) {

			int iblk = _listloc[ilist];
			double *diaggloc;

			int nzloc  = diasize[iblk];

			diaggloc = new double [nzloc];
			if (!diaggloc) MemoryFail (funcname);

			diagg[iblk] = diaggloc;

			int ifileloc  = blk2file[iblk*5+4];
			int ibsloc    = blk2bs  [iblk*5+4];

			FGet (files[ifileloc],nzloc,diaggloc,ibsloc);

		};

	};

// Send

	int nztot = 0;

	int ilist;
	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		nztot += nzloc;

	};

	double *sendarr;

	sendarr = new double [nztot];
	if (!sendarr) MemoryFail (funcname);

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		double *diaggloc = diagg[iblk];

		for (int kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = diaggloc[kii];

		nztot += nzloc;

	};

	CMPIComm commloc = _comm;

	int tag = _msgtag;

	if (nztot > 0) {
		CMPIExchange::Send (commloc, _iproc, tag,
									nztot*sizeof(double), (char *)sendarr);
	};

	delete [] sendarr;

// Free the memory

	if (nfiles > 0) {

		for (int ilist=0;ilist<_nlistloc;ilist++) {

			int iblk = _listloc[ilist];
			double *diaggloc = diagg[iblk];

			delete [] diaggloc;

			diaggloc = new double [0];
			if (!diaggloc) MemoryFail (funcname);

			diagg[iblk] = diaggloc;

		};

	};

};

// Author: Kharchenko S.A.
// CFctDiagR: Receive diagg array
//========================================================================================
void CFctDiagR::ReceiveDiag (const CMPIComm &_comm, char _optype, int _nlistloc, int* _listloc, int _iproc, int _msgtag) { // Receive diagg array

	const char *funcname = "ReceiveDiag";

// Allocate the diagg data for receive

	int nztot = 0;

	int ilist;
	for (ilist=0;ilist<_nlistloc;ilist++) {
		int iblk = _listloc[ilist];
		int nzloc = diasize[iblk];
		nztot += nzloc;
	};

	double *diaggrcv;

	diaggrcv = new double [nztot];
	if (!diaggrcv) MemoryFail (funcname);

	int nz = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		nz += nzloc;

	};

	int tag = _msgtag;

	CMPIComm commloc = _comm;
	CMPIStatus status;

	if (nz > 0) {
		CMPIExchange::Recv (commloc, _iproc, tag,
								nz*sizeof(double), (char *)diaggrcv, status);
	};

// Update

	nz = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		double *diaggloc;

		if (nfiles == 0) {
			diaggloc = diagg[iblk];
		} else {

			diaggloc = new double [nzloc];
			if (!diaggloc) MemoryFail (funcname);

			int ifileloc = blk2file[iblk*5+4];
			int ibsloc   = blk2bs  [iblk*5+4];

			FGet (files[ifileloc],nzloc,diaggloc,ibsloc);

		};

		if (_optype == '+') {
			for (int kii=0;kii<nzloc;kii++) diaggloc[kii] += diaggrcv[nz+kii];
		} else {
			for (int kii=0;kii<nzloc;kii++) diaggloc[kii] = diaggrcv[nz+kii];
		};

		if (nfiles > 0) {

			int ifileloc = blk2file[iblk*5+4];
			int ibsloc   = blk2bs  [iblk*5+4];

			FPut (files[ifileloc],nzloc,diaggloc,ibsloc);

			delete [] diaggloc;

		};

		nz += nzloc;

	};

	delete [] diaggrcv;

};

// Author: Kharchenko S.A.
// CFctDiagC: Send scaling data
//========================================================================================
void CFctDiagC::SendScaling (const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Send scaling data
								const CGSMatrix &_gmtra) {

	const char *funcname = "SendScaling";

	dcmplx *sendarr;

	dcmplx *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	int nztot = 0;

	int ilist, nlocblk, isupbeg, isupend;
	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];

		int nzloc = diasize[iblk];

		nztot += nlocblk+4*nzloc;

	};

	sendarr = new dcmplx [nztot];
	if (!sendarr) MemoryFail (funcname);

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc  = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		int kii;
		for (kii=0;kii<nlocblk;kii++) sendarr[nztot+kii] = sclptloc[kii];
		nztot += nlocblk;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = scllloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = scluloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = sclinvlloc[kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = sclinvuloc[kii];
		nztot += nzloc;

	};

	CMPIComm commloc = _comm;
	int msgtag = _msgtag;

	if (nztot > 0) {
		CMPIExchange::Send (commloc, _iproc, msgtag,
									nztot*sizeof(dcmplx), (char *)sendarr);
	};

	delete [] sendarr;

};

// Author: Kharchenko S.A.
// CFctDiagC: Receive scaling data
//========================================================================================
void CFctDiagC::ReceiveScaling (const CMPIComm &_comm, int _nlistloc, int *_listloc, int _iproc, int _msgtag, // Receive scaling data
									const CGSMatrix &_gmtra) {

	const char *funcname = "ReceiveScaling";

	dcmplx czero (0.0e0,0.0e0);

	dcmplx *recvarr;
	dcmplx *sclptloc, *scllloc, *scluloc, *sclinvlloc, *sclinvuloc;

	int nztot = 0;

	int ilist, nlocblk, isupbeg, isupend;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];

		nztot += nlocblk+4*nzloc;

	};

	recvarr = new dcmplx [nztot];
	if (!recvarr) MemoryFail (funcname);

	CMPIComm commloc = _comm;
	CMPIStatus status;

	int msgtag = _msgtag;

	if (nztot > 0) {
		CMPIExchange::Recv (commloc, _iproc, msgtag,
								nztot*sizeof(dcmplx), (char *)recvarr, status);
	};

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		isupbeg = _gmtra.blksr[iblk];
		isupend = _gmtra.blksr[iblk+1]-1;
		nlocblk = _gmtra.sprndr[isupend+1]-_gmtra.sprndr[isupbeg];

		sclptloc   = sclpt  [iblk];
		scllloc    = scll   [iblk];
		scluloc    = sclu   [iblk];
		sclinvlloc = sclinvl[iblk];
		sclinvuloc = sclinvu[iblk];

		int kii;
		for (kii=0;kii<nlocblk;kii++) sclptloc[kii] = recvarr[nztot+kii];
		nztot += nlocblk;

		for (kii=0;kii<nzloc;kii++) scllloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) scluloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sclinvlloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

		for (kii=0;kii<nzloc;kii++) sclinvuloc[kii] = recvarr[nztot+kii];
		nztot += nzloc;

	};

	delete [] recvarr;

};

// Author: Kharchenko S.A.
// CFctDiagC: Send diagg array
//========================================================================================
void CFctDiagC::SendDiag (const CMPIComm &_comm, int _nlistloc, int* _listloc, int _iproc, int _msgtag, // Send diagg array
							const CGSMatrix &_gmtra) {

	const char *funcname = "SendDiag";

// Read the data from disk if necessary

	if (nfiles > 0) {

		for (int ilist=0;ilist<_nlistloc;ilist++) {

			int iblk = _listloc[ilist];
			dcmplx *diaggloc;

			int nzloc  = diasize[iblk];

			diaggloc = new dcmplx [nzloc];
			if (!diaggloc) MemoryFail (funcname);

			diagg[iblk] = diaggloc;

			int ifileloc  = blk2file[iblk*5+4];
			int ibsloc    = blk2bs  [iblk*5+4];

			FGet (files[ifileloc],nzloc,diaggloc,ibsloc);

		};

	};

// Send

	int nztot = 0;

	int ilist;
	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		nztot += nzloc;

	};

	dcmplx *sendarr;

	sendarr = new dcmplx [nztot];
	if (!sendarr) MemoryFail (funcname);

	nztot = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		dcmplx *diaggloc = diagg[iblk];

		for (int kii=0;kii<nzloc;kii++) sendarr[nztot+kii] = diaggloc[kii];

		nztot += nzloc;

	};

	CMPIComm commloc = _comm;
	int tag = _msgtag;

	if (nztot > 0) {
		CMPIExchange::Send (commloc, _iproc, tag,
									nztot*sizeof(dcmplx), (char *)sendarr);
	};

	delete [] sendarr;

// Free the memory

	if (nfiles > 0) {

		for (int ilist=0;ilist<_nlistloc;ilist++) {

			int iblk = _listloc[ilist];
			dcmplx *diaggloc = diagg[iblk];

			delete [] diaggloc;

			diaggloc = new dcmplx [0];
			if (!diaggloc) MemoryFail (funcname);

			diagg[iblk] = diaggloc;

		};

	};

};

// Author: Kharchenko S.A.
// CFctDiagC: Receive diagg array
//========================================================================================
void CFctDiagC::ReceiveDiag (const CMPIComm &_comm, char _optype, int _nlistloc, int* _listloc, int _iproc, int _msgtag) { // Receive diagg array

	const char *funcname = "ReceiveDiag";

// Allocate the diagg data for receive

	int nztot = 0;

	int ilist;
	for (ilist=0;ilist<_nlistloc;ilist++) {
		int iblk = _listloc[ilist];
		int nzloc = diasize[iblk];
		nztot += nzloc;
	};

	dcmplx *diaggrcv;

	diaggrcv = new dcmplx [nztot];
	if (!diaggrcv) MemoryFail (funcname);

	int nz = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		nz += nzloc;

	};

	int tag = _msgtag;

	CMPIComm commloc = _comm;
	CMPIStatus status;

	if (nz > 0) {
		CMPIExchange::Recv (commloc, _iproc, tag,
								nz*sizeof(dcmplx), (char *)diaggrcv, status);
	};

// Update

	nz = 0;

	for (ilist=0;ilist<_nlistloc;ilist++) {

		int iblk = _listloc[ilist];

		int nzloc = diasize[iblk];

		dcmplx *diaggloc;

		if (nfiles == 0) {
			diaggloc = diagg[iblk];
		} else {

			diaggloc = new dcmplx [nzloc];
			if (!diaggloc) MemoryFail (funcname);

			int ifileloc = blk2file[iblk*5+4];
			int ibsloc   = blk2bs  [iblk*5+4];

			FGet (files[ifileloc],nzloc,diaggloc,ibsloc);

		};

		if (_optype == '+') {
			for (int kii=0;kii<nzloc;kii++) diaggloc[kii] += diaggrcv[nz+kii];
		} else {
			for (int kii=0;kii<nzloc;kii++) diaggloc[kii] = diaggrcv[nz+kii];
		};

		if (nfiles > 0) {

			int ifileloc = blk2file[iblk*5+4];
			int ibsloc   = blk2bs  [iblk*5+4];

			FPut (files[ifileloc],nzloc,diaggloc,ibsloc);

			delete [] diaggloc;

		};

		nz += nzloc;

	};

	delete [] diaggrcv;

};

// Author: Kharchenko S.A.
// CMvmR: Receive vector data from prescribed processor
//========================================================================================
void CMvmR::ReceiveVectData (bool _is2index, const CMPIComm &_comm, char _recvtype, char _optype, // Receive vector data from prescribed processor
								int _myid, int _iproc, int _ishifttag, CSVector &_px,
								int _nlist, int *_list,
								int *_blksr, int *_sprndr, int *_bl2ndr) {

// Parameter _recvtype takes the following values:
//
// 'A' - receive all data of the list;
// 'L' - receive only those data of the list that belong to the current processor (id=_myid);
// 'l' - receive all data of the list except the one that belong to the current processor (id=_myid);
// 'R' - receive only those data of the list that belong to the remote processor (id=_iproc);
// 'r' - receive all data of the list except the one that belong to the remote processor (id=_iproc);

//	const char *funcname = "ReceiveVectData";

// Receive

	int nz = 0;

	int irecv = 0;

	double *pxrecv = vectsndrcv;

	int ilist;
	for (ilist=0;ilist<_nlist;ilist++) {

		int iblk = _list[ilist];

		if (_recvtype == 'A' || (_recvtype == 'L' && blk2cpu[iblk] == _myid) || (_recvtype == 'l' && blk2cpu[iblk] != _myid) ||
			(_recvtype == 'R' && blk2cpu[iblk] == _iproc) || (_recvtype == 'r' && blk2cpu[iblk] != _iproc)) {

			irecv++;
//			cout << " Recv vect from " << _iproc << " to " << _myid << " iblk = " << iblk << endl;

			int ni = _bl2ndr[iblk+1]-_bl2ndr[iblk];

			int nzloc = ni*nrhs;

			nz += nzloc;

		};

	};

	CMPIComm commloc = _comm;
	CMPIStatus status;

	if (nz > 0) {
		CMPIExchange::Recv (commloc, _iproc, _ishifttag,
								nz*sizeof(double), (char *)pxrecv, status);
	};

	int kii;

	nz = 0;

	for (ilist=0;ilist<_nlist;ilist++) {

		int iblk = _list[ilist];

		if (_recvtype == 'A' || (_recvtype == 'L' && blk2cpu[iblk] == _myid) || (_recvtype == 'l' && blk2cpu[iblk] != _myid) ||
			(_recvtype == 'R' && blk2cpu[iblk] == _iproc) || (_recvtype == 'r' && blk2cpu[iblk] != _iproc)) {

			int ni = _bl2ndr[iblk+1]-_bl2ndr[iblk];

			if (_optype == '=') {

				if (!_is2index) {

					for (int irhs=0;irhs<nrhs;irhs++) {
						int niloc=0;
						for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
							int blai=_sprndr[isup+1]-_sprndr[isup];
							int ibs = bspx[isup];
							for (int kii=0;kii<blai;kii++) {
								_px.vect[nlocext*irhs+ibs+kii] = pxrecv[nz+ni*irhs+niloc+kii];
							};
							niloc += blai;
						};
					};

				} else {

					for (int irhs=0;irhs<nrhs;irhs++) {
						int ibs = bspx[iblk];
						for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
							kii = isup-_blksr[iblk];
							_px.vect[nlocext*irhs+ibs+kii] = pxrecv[nz+ni*irhs+kii];
						};
					};
				};

			} else {

				if (!_is2index) {

					for (int irhs=0;irhs<nrhs;irhs++) {
						int niloc=0;
						for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
							int blai=_sprndr[isup+1]-_sprndr[isup];
							int ibs = bspx[isup];
							for (int kii=0;kii<blai;kii++) {
								_px.vect[nlocext*irhs+ibs+kii] += pxrecv[nz+ni*irhs+niloc+kii];
							};
							niloc += blai;
						};
					};

				} else {

					for (int irhs=0;irhs<nrhs;irhs++) {
						int ibs = bspx[iblk];
						for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
							kii = isup-_blksr[iblk];
							_px.vect[nlocext*irhs+ibs+kii] += pxrecv[nz+ni*irhs+kii];
						};
					};

				};

			};

			int nzloc = ni*nrhs;

			nz += nzloc;

		};

	};

};

// Author: Kharchenko S.A.
// CMvmR: Send vector data from prescribed processor
//========================================================================================
void CMvmR::SendVectData (bool _is2index, const CMPIComm &_comm, char _sendtype, // Send vector data to the prescribed processor
							int _myid, int _iproc, int _ishifttag, CSVector &_px,
							int _nlist, int *_list,
							int *_blksr, int *_sprndr, int *_bl2ndr) {

// Parameter _sendtype takes the following values:
//
// 'A' - send all data of the list;
// 'L' - send only those data of the list that belong to the current processor (id=_myid);
// 'l' - send all data of the list except the one that belong to the current processor (id=_myid);
// 'R' - send only those data of the list that belong to the remote processor (id=_iproc);
// 'r' - send all data of the list except the one that belong to the remote processor (id=_iproc);

//	const char *funcname = "SendVectData";

// Prepare and send data to the child nodes if necessary

	double *pxsend = vectsndrcv;

	int isend = 0;

	int nz = 0;

	int ilist, kii;

	for (ilist=0;ilist<_nlist;ilist++) {

		int iblk = _list[ilist];

		if (_sendtype == 'A' || (_sendtype == 'L' && blk2cpu[iblk] == _myid) || (_sendtype == 'l' && blk2cpu[iblk] != _myid) ||
			(_sendtype == 'R' && blk2cpu[iblk] == _iproc) || (_sendtype == 'r' && blk2cpu[iblk] != _iproc)) {

//			cout << " Send vect to " << _iproc << " from " << _myid << " iblk = " << iblk << endl;
			isend++;

			int ni = _bl2ndr[iblk+1]-_bl2ndr[iblk];

			if (!_is2index) {

				for (int irhs=0;irhs<nrhs;irhs++) {
					int niloc=0;
					for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
						int blai=_sprndr[isup+1]-_sprndr[isup];
						int ibs = bspx[isup];
						for (int kii=0;kii<blai;kii++) {
							pxsend[nz+ni*irhs+niloc+kii] = _px.vect[nlocext*irhs+ibs+kii];
						};
						niloc += blai;
					};
				};

			} else {

				for (int irhs=0;irhs<nrhs;irhs++) {
					int ibs = bspx[iblk];
					for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
						kii = isup-_blksr[iblk];
						pxsend[nz+ni*irhs+kii] = _px.vect[nlocext*irhs+ibs+kii];
					};
				};

			};

			int nzloc  = ni*nrhs;

			nz += nzloc;

		};

	};

	CMPIComm commloc = _comm;

	if (nz > 0) {
		CMPIExchange::Send (commloc, _iproc, _ishifttag,
									nz*sizeof(double), (char *)pxsend);
	};

};

// Author: Kharchenko S.A.
// CMvmC: Receive vector data from prescribed processor
//========================================================================================
void CMvmC::ReceiveVectData (const CMPIComm &_comm, char _recvtype, char _optype, // Receive vector data from prescribed processor
								int _myid, int _iproc, int _ishifttag, CSVectorC &_px,
								int _nlist, int *_list,
								int *_blksr, int *_sprndr, int *_bl2ndr) {

// Parameter _recvtype take the following values:
//
// 'A' - receive all data of the list;
// 'L' - receive only those data of the list that belong to the current processor (id=_myid);
// 'l' - receive all data of the list except the one that belong to the current processor (id=_myid);
// 'R' - receive only those data of the list that belong to the remote processor (id=_iproc);
// 'r' - receive all data of the list except the one that belong to the remote processor (id=_iproc);

//	const char *funcname = "ReceiveVectData";

// Receive

	int nz = 0;

	int irecv = 0;

	dcmplx *pxrecv = vectsndrcv;

	int ilist;
	for (ilist=0;ilist<_nlist;ilist++) {

		int iblk = _list[ilist];

		if (_recvtype == 'A' || (_recvtype == 'L' && blk2cpu[iblk] == _myid) || (_recvtype == 'l' && blk2cpu[iblk] != _myid) ||
			(_recvtype == 'R' && blk2cpu[iblk] == _iproc) || (_recvtype == 'r' && blk2cpu[iblk] != _iproc)) {

			irecv++;
//			cout << " Recv vect from " << _iproc << " to " << _myid << " iblk = " << iblk << endl;

			int ni = _bl2ndr[iblk+1]-_bl2ndr[iblk];

			int nzloc = ni*nrhs;

			nz += nzloc;

		};

	};

	CMPIComm commloc = _comm;
	CMPIStatus status;

	if (nz > 0) {
		CMPIExchange::Recv (commloc, _iproc, _ishifttag,
								nz*sizeof(dcmplx), (char *)pxrecv, status);
	};

	nz = 0;

	for (ilist=0;ilist<_nlist;ilist++) {

		int iblk = _list[ilist];

		if (_recvtype == 'A' || (_recvtype == 'L' && blk2cpu[iblk] == _myid) || (_recvtype == 'l' && blk2cpu[iblk] != _myid) ||
			(_recvtype == 'R' && blk2cpu[iblk] == _iproc) || (_recvtype == 'r' && blk2cpu[iblk] != _iproc)) {

			int ni = _bl2ndr[iblk+1]-_bl2ndr[iblk];

			if (_optype == '=') {

				for (int irhs=0;irhs<nrhs;irhs++) {
					int niloc=0;
					for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
						int blai=_sprndr[isup+1]-_sprndr[isup];
						int ibs = bspx[isup];
						for (int kii=0;kii<blai;kii++) {
							_px.vect[nlocext*irhs+ibs+kii] = pxrecv[nz+ni*irhs+niloc+kii];
						};
						niloc += blai;
					};
				};

			} else {

				for (int irhs=0;irhs<nrhs;irhs++) {
					int niloc=0;
					for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
						int blai=_sprndr[isup+1]-_sprndr[isup];
						int ibs = bspx[isup];
						for (int kii=0;kii<blai;kii++) {
							_px.vect[nlocext*irhs+ibs+kii] += pxrecv[nz+ni*irhs+niloc+kii];
						};
						niloc += blai;
					};
				};

			};

			int nzloc = ni*nrhs;

			nz += nzloc;

		};

	};

};

// Author: Kharchenko S.A.
// CMvmC: Send vector data from prescribed processor
//========================================================================================
void CMvmC::SendVectData (const CMPIComm &_comm, char _sendtype, // Send vector data to the prescribed processor
							int _myid, int _iproc, int _ishifttag, CSVectorC &_px,
							int _nlist, int *_list,
							int *_blksr, int *_sprndr, int *_bl2ndr) {

// Parameter _sendtype take the following values:
//
// 'A' - send all data of the list;
// 'L' - send only those data of the list that belong to the current processor (id=_myid);
// 'l' - send all data of the list except the one that belong to the current processor (id=_myid);
// 'R' - send only those data of the list that belong to the remote processor (id=_iproc);
// 'r' - send all data of the list except the one that belong to the remote processor (id=_iproc);

//	const char *funcname = "SendVectData";

// Prepare and send data to the child nodes if necessary

	dcmplx *pxsend = vectsndrcv;

	int isend = 0;

	int nz = 0;

	int ilist;

	for (ilist=0;ilist<_nlist;ilist++) {

		int iblk = _list[ilist];

		if (_sendtype == 'A' || (_sendtype == 'L' && blk2cpu[iblk] == _myid) || (_sendtype == 'l' && blk2cpu[iblk] != _myid) ||
			(_sendtype == 'R' && blk2cpu[iblk] == _iproc) || (_sendtype == 'r' && blk2cpu[iblk] != _iproc)) {

//			cout << " Send vect to " << _iproc << " from " << _myid << " iblk = " << iblk << endl;
			isend++;

			int ni = _bl2ndr[iblk+1]-_bl2ndr[iblk];

			for (int irhs=0;irhs<nrhs;irhs++) {
				int niloc=0;
				for (int isup=_blksr[iblk];isup<_blksr[iblk+1];isup++) {
					int blai=_sprndr[isup+1]-_sprndr[isup];
					int ibs = bspx[isup];
					for (int kii=0;kii<blai;kii++) {
						pxsend[nz+ni*irhs+niloc+kii] = _px.vect[nlocext*irhs+ibs+kii];
					};
					niloc += blai;
				};
			};

			int nzloc  = ni*nrhs;

			nz += nzloc;

		};

	};

	CMPIComm commloc = _comm;

	if (nz > 0) {
		CMPIExchange::Send (commloc, _iproc, _ishifttag,
									nz*sizeof(dcmplx), (char *)pxsend);
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Send matrix
//========================================================================================
void CGSMatrixR::SendMatrix2Index (const CMPIComm &_comm, int _iproc, int _msgtag, const CSMatrixR &_mtra) const { // Send matrix

// Pack the matrix

	int length;
	char *obj;

	CSMatrixR mtraloc = _mtra;

	mtraloc.PackMatrix (length, obj);

	CSMatrixR mtrdummy;

	mtraloc = mtrdummy;

// Send dimension

	int msgtag = 2*_msgtag;

	CMPIComm commloc = _comm;

	CMPIExchange::Send (commloc, _iproc, msgtag,
								sizeof(int), (char *)&length);

// Send matrix itself

	msgtag = 2*_msgtag+1;

	CMPIExchange::Send (commloc, _iproc, msgtag,
								length, obj);

// Free work data

	delete [] obj;

};

// Author: Kharchenko S.A.
// CGSMatrixR: Receive matrix
//========================================================================================
void CGSMatrixR::ReceiveMatrix2Index (const CMPIComm &_comm, int _iproc, int _msgtag, CSMatrixR &_mtra) const { // Send matrix

	const char *funcname = "ReceiveMatrix2Index";

// Receive the length of the array

	int length;

	int msgtag = 2*_msgtag;

	CMPIComm commloc = _comm;
	CMPIStatus status;

	CMPIExchange::Recv (commloc, _iproc, msgtag,
								sizeof(int), (char *)&length, status);

// Receive the data

	char *obj;

	obj = new char [length];
	if (!obj) MemoryFail (funcname);

	msgtag = 2*_msgtag+1;

	CMPIExchange::Recv (commloc, _iproc, msgtag,
								length, obj, status);

	_mtra.UnPackMatrix (length, obj);

	delete [] obj;

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Send matrix
//========================================================================================
void CGSMatrixRS::SendMatrix (const CMPIComm &_comm, int _iproc, int _msgtag, const CSMatrixRS &_mtra) const { // Send matrix

// Pack the matrix

	int length;
	char *obj;

	CSMatrixRS mtraloc = _mtra;

	mtraloc.PackMatrix (length, obj);

	CSMatrixRS mtrdummy;

	mtraloc = mtrdummy;

// Send dimension

	int msgtag = 2*_msgtag;

	CMPIComm commloc = _comm;

	CMPIExchange::Send (commloc, _iproc, msgtag,
								sizeof(int), (char *)&length);

// Send matrix itself

	msgtag = 2*_msgtag+1;

	CMPIExchange::Send (commloc, _iproc, msgtag,
								length, obj);

// Free work data

	delete [] obj;

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Receive matrix
//========================================================================================
void CGSMatrixRS::ReceiveMatrix (const CMPIComm &_comm, int _iproc, int _msgtag, CSMatrixRS &_mtra) const { // Send matrix

	const char *funcname = "ReceiveMatrix";

// Receive the length of the array

	int length;

	int msgtag = 2*_msgtag;

	CMPIComm commloc = _comm;
	CMPIStatus status;

	CMPIExchange::Recv (commloc, _iproc, msgtag,
								sizeof(int), (char *)&length, status);

// Receive the data

	char *obj;

	obj = new char [length];
	if (!obj) MemoryFail (funcname);

	msgtag = 2*_msgtag+1;

	CMPIExchange::Recv (commloc, _iproc, msgtag,
								length, obj, status);

	_mtra.UnPackMatrix (length, obj);

	delete [] obj;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Send matrix
//========================================================================================
void CGSMatrixCS::SendMatrix (const CMPIComm &_comm, int _iproc, int _msgtag, const CSMatrixCS &_mtra) const { // Send matrix

// Pack the matrix

	int length;
	char *obj;

	CSMatrixCS mtraloc = _mtra;

	mtraloc.PackMatrix (length, obj);

	CSMatrixCS mtrdummy;

	mtraloc = mtrdummy;

// Send dimension

	int msgtag = 2*_msgtag;

	CMPIComm commloc = _comm;

	CMPIExchange::Send (commloc, _iproc, msgtag,
								sizeof(int), (char *)&length);

// Send matrix itself

	msgtag = 2*_msgtag+1;

	CMPIExchange::Send (commloc, _iproc, msgtag,
								length, obj);

// Free work data

	delete [] obj;

};

// Author: Kharchenko S.A.
// CGSMatrixCS: Receive matrix
//========================================================================================
void CGSMatrixCS::ReceiveMatrix (const CMPIComm &_comm, int _iproc, int _msgtag, CSMatrixCS &_mtra) const { // Send matrix

	const char *funcname = "ReceiveMatrix";

// Receive the length of the array

	int length;

	int msgtag = 2*_msgtag;

	CMPIComm commloc = _comm;
	CMPIStatus status;

	CMPIExchange::Recv (commloc, _iproc, msgtag,
								sizeof(int), (char *)&length, status);

// Receive the data

	char *obj;

	obj = new char [length];
	if (!obj) MemoryFail (funcname);

	msgtag = 2*_msgtag+1;

	CMPIExchange::Recv (commloc, _iproc, msgtag,
								length, obj, status);

	_mtra.UnPackMatrix (length, obj);

	delete [] obj;

};

// Author: Kharchenko S.A.
// Perform requests allocation
//========================================================================================
void AllocateRequests (int _nsends, int _nrecvs, CMPIRequest *&_sndreqs, CMPIRequest *&_rcvreqs, // Perform requests allocation
						CMPIStatus *&_sndstat, CMPIStatus *&_rcvstat, const char *_funcname) {

//	const char *funcname = "AllocateRequests";

	_sndreqs = new CMPIRequest [_nsends+1];
	if (!_sndreqs) MemoryFail (_funcname);
	_rcvreqs = new CMPIRequest [_nrecvs+1];
	if (!_rcvreqs) MemoryFail (_funcname);
	_sndstat = new CMPIStatus  [_nsends+1];
	if (!_sndstat) MemoryFail (_funcname);
	_rcvstat = new CMPIStatus  [_nrecvs+1];
	if (!_rcvstat) MemoryFail (_funcname);

};

// Author: Kharchenko S.A.
// Wait for any receive from the list
//========================================================================================
void WaitReceive (int _nrecvs, CMPIRequest *_rcvreqs, int &_index, CMPIStatus *_rcvstat, // Wait for any receive from the list
					char *_funcname) {

	CMPIExchange::WaitAny (_nrecvs, _rcvreqs, _index, _rcvstat[0]);

};

// Author: Kharchenko S.A.
// Wait all sends for completion and free all requests
//========================================================================================
void WaitSendFreeRequests (int _nsends, int _nrecvs, CMPIRequest *_sndreqs, CMPIRequest *_rcvreqs, // Wait all sends for completion and free all requests
						CMPIStatus *_sndstat, CMPIStatus *_rcvstat, const char *_funcname) {

// Wait till all sends will be completed

	CMPIExchange::WaitAll (_nsends, _sndreqs, _sndstat);

// Free the memory

	delete [] _sndreqs;
	delete [] _rcvreqs;
	delete [] _sndstat;
	delete [] _rcvstat;

};

// Author: Kharchenko S.A.
// CSVectorC: Collect solution vector on root processor
//========================================================================================
void CSVectorC::CollectSolutionOnRoot (const CTree &_tree, const CGSMatrixCS &_gmtra, // Collect solution vector on root processor
										int *_bl2cpu, int _nrhs,
										dcmplx *_solloc, dcmplx *_soltot) {

	const char *funcname = "CollectSolutionOnRoot";

// Compute for each processor the local size of the block

	int nproc = _tree.nproc;

	int *len2cpu;

	len2cpu = new int [nproc+1];
	if (!len2cpu) MemoryFail (funcname);

	int i, iblk, iproc;

	for (i=0;i<nproc;i++) len2cpu[i] = 0;

	for (iblk=0;iblk<_gmtra.nblksc;iblk++) {
		iproc = _bl2cpu[iblk];
		len2cpu[iproc] += _gmtra.bl2ndc[iblk+1] - _gmtra.bl2ndc[iblk];
	};

// Allocate the memory for receives

	int *ibssol;
	dcmplx *solarr;

	ibssol = new int [nproc+1];
	if (!ibssol) MemoryFail (funcname);
	solarr = new dcmplx [_gmtra.n*_nrhs];
	if (!solarr) MemoryFail (funcname);

	ibssol[0] = 0;
	for (i=0;i<nproc;i++) ibssol[i+1] = ibssol[i] + len2cpu[i]*_nrhs;

// Send / receive data

	CMPIComm pcomm = _tree.GetComm ();

	if (_tree.myid == _tree.rootcpu) {

		int nsends = 0;
		int nreqvs = _tree.nproc-1;

		CMPIRequest *sndreqs=NULL, *rcvreqs=NULL;
		CMPIStatus  *sndstat=NULL, *rcvstat=NULL;

		AllocateRequests (nsends, nreqvs, sndreqs, rcvreqs, 
							sndstat, rcvstat, funcname);

		CMPIComm commloc = pcomm;

		int icount = 0;

		for (iproc=0;iproc<nproc;iproc++) {

			if (iproc != _tree.myid) {

				int tag = iproc;

				int ibs = ibssol[iproc];

//				MPI_Irecv (solarr+ibs, 2*len2cpu[iproc]*nrhs, MPI_DOUBLE, iproc, tag, 
//										*((MPI_Comm*)pcomm.comm), rcvreqs+icount);
				CMPIExchange::IRecv (commloc, iproc, tag,
											len2cpu[iproc]*nrhs*sizeof(dcmplx), (char *)(solarr+ibs), rcvreqs[icount]);


				icount++;

			} else {

				int ibs = ibssol[iproc];

				for (int kii=0;kii<len2cpu[iproc]*_nrhs;kii++) solarr[ibs+kii] = _solloc[kii];

			};

		};

		WaitSendFreeRequests (nreqvs, nsends, rcvreqs, sndreqs, rcvstat, sndstat, funcname);

	} else {

		int msgtag = _tree.myid;

		CMPIComm commloc = pcomm;

		CMPIExchange::Send (commloc, _tree.rootcpu, msgtag,
									len2cpu[_tree.myid]*_nrhs*sizeof(dcmplx), (char *)_solloc);

	};

// Store the data in the global array

	int ntotal = _gmtra.n;
	for (iproc=0;iproc<nproc;iproc++) {
		int mloc = 0;
		int ldloc = len2cpu[iproc];
		int ibs = ibssol[iproc];
		for (iblk=0;iblk<_gmtra.nblksc;iblk++) {
			if (_bl2cpu[iblk] == iproc) {
				int ni = _gmtra.bl2ndc[iblk+1]-_gmtra.bl2ndc[iblk];
				for (int irhs=0;irhs<_nrhs;irhs++) {
					for (int kii=0;kii<ni;kii++) {
						_soltot[irhs*ntotal+_gmtra.bl2ndc[iblk]+kii] = solarr[ibs+irhs*ldloc+mloc+kii];
					};
				};
				mloc += ni;
			};
		};
	};

// Free work memory

	delete [] len2cpu;
	delete [] ibssol;
	delete [] solarr;

};

// Author: Kharchenko S.A.
// Update QR decomposition of the complex block in parallel mode
//========================================================================================
void CQrd::UpdateQrd (const CTree &_tree, // Update QR decomposition of the block in parallel mode
						int _m, int _ibeg, int _iend, double *_a, int _lda, double *_tau) {

	const char *funcname = "UpdateQrd";

	double dzero = 0.0e0;

// Check input parameters

	if (_iend-_ibeg != blksiz) throw " Wrong ibeg and iend parameters in UpdateQrd function";
	if (_ibeg%blksiz != 0) throw " Wrong iend parameter in UpdateQrd function";

	int iblk=_ibeg / blksiz;

// Update first local part of the decomposition

	if (nfiles == 0) {

//		flowvision::EqnSolver::UpdateQrd (_m, _ibeg, _iend, _a, _lda, _tau);
		::UpdateQrd (_m, _ibeg, _iend, _a, _lda, _tau);

	} else {

		int ifile = bl2file[iblk];
		int ibs = bsblk[iblk];

		FPut (files[ifile], blksiz*_lda, _a, ibs);

//		flowvision::EqnSolver::UpdateQrd (_m, iblk, blksiz,
		::UpdateQrd (_m, iblk, blksiz,
						files, bl2file, bsblk, 
						_lda, _tau,
						_a);

	};

// Get R part of the data and put it into send array

	int ldqsend = (iblk+1)*blksiz;

	if (nfiles == 0) {
//		flowvision::EqnSolver::GetRPartQrd (_m, _ibeg, _iend, _a, _lda, qsendarr, ldqsend);
		::GetRPartQrd (_m, _ibeg, _iend, _a, _lda, qsendarr, ldqsend);
	} else {
//		flowvision::EqnSolver::GetRPartQrd (_m, iblk, iblk+1, blksiz, 
		::GetRPartQrd (_m, iblk, iblk+1, blksiz, 
						files, bl2file, bsblk, 
						_lda, qsendarr, ldqsend);
	};

// Allocate work array

	double *qblock;

	qblock = new double [(iblk+1)*blksiz*blksiz*nchildsmax];
	if (!qblock) MemoryFail (funcname);

// Search the tree

	CMPIComm pcomm = _tree.GetComm ();

	int nodeend = _tree.cpuidend[_tree.myid];

	int inode = nodeend;

	while (inode >= 0) {

		int ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

		int fathernode = _tree.nodes[inode].fatherid;
		int fathercpuid = _tree.nodes[inode].fathercpu;

// Receive data from the child nodes and init factorization data

		int ldq = (iblk+1)*blksiz*nchilds;

		for (int kii=0;kii<ldq*blksiz;kii++) qblock[kii] = dzero;

		for (int ichild=0;ichild<nchilds;ichild++) {

			int childnode = _tree.nodes[inode].childs[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			double *qrcv;

			if (childcpuid == _tree.myid) {
				qrcv = qsendarr;
			} else {

				CMPIStatus status;

				int msgtag = iblk+1;

				CMPIExchange::Recv (pcomm, childcpuid, msgtag,
											ldqsend*blksiz*sizeof(double), (char *)qrecvarr, status);

				qrcv = qrecvarr;

			};

			for (int jblk=0;jblk<=iblk;jblk++) {
				for (int irhs=0;irhs<blksiz;irhs++) {
					for (int kii=0;kii<blksiz;kii++) {
						qblock[irhs*ldq+jblk*blksiz*nchilds+blksiz*ichild+kii] = qrcv[irhs*ldqsend+jblk*blksiz+kii];
					};
				}
			};

		};

// Update QR decomposition

		UpdateQrdBlk (ilev, iblk, ldq, qblock);

// Send newly computed R data to the father cpu if necessary

		GetRPartQrd (_tree, ilev, _ibeg, _iend, qsendarr, ldqsend);

		if (fathercpuid != _tree.myid) {

			int msgtag = iblk+1;

			CMPIExchange::Send (pcomm, fathercpuid, msgtag,
										ldqsend*blksiz*sizeof(double), (char *)qsendarr);

		};

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			inode = fathernode;
			if (_tree.nodes[inode].nodecpu != _tree.myid) inode = -1;
		};

	};

// Free work arrays

	delete [] qblock;

};

// Author: Kharchenko S.A.
// Multiply Q factor by the current block in the parallel mode
//========================================================================================
void CQrd::MvmQBlk (const CTree &_tree, // Multiply Q factor by the current block in the parallel mode
						int _m, int _n, int _nrhs,
						double *_q, int _ldq, double *_tau, 
						double *_x, int _ldx, double *_qx, int _ldqx) { 

	const char *funcname = "MvmQBlk";

// Check input parameters

	if (_n%blksiz != 0) throw " Wrong n parameter in MvmQBlk function";

	if (_nrhs != blksiz) throw " Wrong nrhs parameter in MvmQBlk function";

	int iblk=_n / blksiz - 1;

// Allocate work array

	double *qblock;

	qblock = new double [(iblk+1)*blksiz*blksiz*nchildsmax];
	if (!qblock) MemoryFail (funcname);

// Assign initial qrecvarr

	int ldqrecv = (iblk+1)*blksiz;

	for (int irhs=0;irhs<blksiz;irhs++) {
		for (int kii=0;kii<_n;kii++) {
			qrecvarr[irhs*ldqrecv+kii] = _x[irhs*_ldx+kii];
		};
	};

// Search the tree

	int nodebeg = _tree.cpuidbeg[_tree.myid];

	CMPIComm pcomm = _tree.GetComm ();

	int inode = nodebeg;

	while (inode >= 0) {

		int ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

//		int fathernode = _tree.nodes[inode].fatherid;
		int fathercpuid = _tree.nodes[inode].fathercpu;

// Receive data from father cpu if necessary

		if (fathercpuid != _tree.myid) {

			CMPIStatus status;

			int msgtag = iblk+1;

			CMPIExchange::Recv (pcomm, fathercpuid, msgtag,
										ldqrecv*blksiz*sizeof(double), (char *)qrecvarr, status);

		};

// Multiply

		int ldq = (iblk+1)*blksiz*nchilds;

		MvmQBlk (ilev, iblk, _nrhs, qrecvarr, ldqrecv , qblock, ldq);

// Init multiplication data and send them to the child nodes

		for (int ichild=0;ichild<nchilds;ichild++) {

			int childnode = _tree.nodes[inode].childs[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			double *qsnd;

			if (childcpuid == _tree.myid) {
				qsnd = qrecvarr;
			} else {
				qsnd = qsendarr;
			};

			int ldqsend = (iblk+1)*blksiz;

			for (int jblk=0;jblk<=iblk;jblk++) {
				for (int irhs=0;irhs<blksiz;irhs++) {
					for (int kii=0;kii<blksiz;kii++) {
						qsnd[irhs*ldqsend+jblk*blksiz+kii] = qblock[irhs*ldq+jblk*blksiz*nchilds+blksiz*ichild+kii];
					};
				}
			};

			if (childcpuid != _tree.myid) {

				int msgtag = iblk+1;

				CMPIExchange::Send (pcomm, childcpuid, msgtag,
											ldqsend*blksiz*sizeof(double), (char *)qsendarr);

			};

		};

// Go down level

		if (nchilds == 1) {
			inode = -1;
		} else {
			int inodecurr = inode;
			inode = -1;
			for (int ichild=0;ichild<nchilds;ichild++) {
				int childnode = _tree.nodes[inodecurr].childs[ichild];
				if (_tree.nodes[childnode].nodecpu == _tree.myid) {
					inode = childnode;
					break;
				};
			};
		};

	};

// Perform multiplication of the local part of the decomposition

	if (nfiles == 0) {
		::MvmQBlk (_m, _n, _nrhs, _q, _ldq, _tau, 
					qrecvarr, ldqrecv, _qx, _ldqx);
	} else {

		int nblk = _n / _nrhs;

		::MvmQBlk (_m, nblk, blksiz, _nrhs,
					files, bl2file, bsblk, _ldq, _tau,
					qrecvarr, ldqrecv, _qx, _ldqx,
					_q);
	};

// Free work arrays

	delete [] qblock;

};

// Author: Kharchenko S.A.
// Update QR decomposition of the complex block in parallel mode
//========================================================================================
void CQrdC::UpdateQrd (const CTree &_tree, // Update QR decomposition of the block in parallel mode
						int _m, int _ibeg, int _iend, dcmplx *_a, int _lda, dcmplx *_tau) {

	const char *funcname = "UpdateQrd";

	dcmplx czero (0.0e0,0.0e0);

// Check input parameters

//	int myid = _tree.GetMyid ();

//	char strbuff[512];
//	sprintf(strbuff, "%s%d%s","ChkQrd_",myid,".dat");
//	ofstream fout (strbuff,ios::app);

	if (_iend-_ibeg != blksiz) throw " Wrong ibeg and iend parameters in UpdateQrd function";
	if (_ibeg%blksiz != 0) throw " Wrong iend parameter in UpdateQrd function";

	int iblk=_ibeg / blksiz;

// Update first local part of the decomposition

	if (nfiles == 0) {

//		flowvision::EqnSolver::UpdateQrd (_m, _ibeg, _iend, _a, _lda, _tau);
		::UpdateQrd (_m, _ibeg, _iend, _a, _lda, _tau);

	} else {

		int ifile = bl2file[iblk];
		int ibs = bsblk[iblk];

		FPut (files[ifile], blksiz*_lda, _a, ibs);

//		flowvision::EqnSolver::UpdateQrd (_m, iblk, blksiz,
		::UpdateQrd (_m, iblk, blksiz,
						files, bl2file, bsblk, 
						_lda, _tau,
						_a);
	};

// Get R part of the data and put it into send array

	int ldqsend = (iblk+1)*blksiz;

	if (nfiles == 0) {
//		flowvision::EqnSolver::GetRPartQrd (_m, _ibeg, _iend, _a, _lda, qsendarr, ldqsend);
		::GetRPartQrd (_m, _ibeg, _iend, _a, _lda, qsendarr, ldqsend);
	} else {
//		flowvision::EqnSolver::GetRPartQrd (_m, iblk, iblk+1, blksiz, 
		::GetRPartQrd (_m, iblk, iblk+1, blksiz, 
						files, bl2file, bsblk, 
						_lda, qsendarr, ldqsend);
	};

// Allocate work array

	dcmplx *qblock;

	qblock = new dcmplx [(iblk+1)*blksiz*blksiz*nchildsmax];
	if (!qblock) MemoryFail (funcname);

// Search the tree

	CMPIComm pcomm = _tree.GetComm ();

	int nodeend = _tree.cpuidend[_tree.myid];

	int inode = nodeend;

	while (inode >= 0) {

		int ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

		int fathernode = _tree.nodes[inode].fatherid;
		int fathercpuid = _tree.nodes[inode].fathercpu;

// Receive data from the child nodes and init factorization data

		int ldq = (iblk+1)*blksiz*nchilds;

		for (int kii=0;kii<ldq*blksiz;kii++) qblock[kii] = czero;

		for (int ichild=0;ichild<nchilds;ichild++) {

			int childnode = _tree.nodes[inode].childs[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			dcmplx *qrcv;

			if (childcpuid == _tree.myid) {
				qrcv = qsendarr;
			} else {

				CMPIStatus status;

				int msgtag = iblk+1;

				CMPIExchange::Recv (pcomm, childcpuid, msgtag,
											ldqsend*blksiz*sizeof(dcmplx), (char *)qrecvarr, status);

				qrcv = qrecvarr;

			};

			for (int jblk=0;jblk<=iblk;jblk++) {
				for (int irhs=0;irhs<blksiz;irhs++) {
					for (int kii=0;kii<blksiz;kii++) {
						qblock[irhs*ldq+jblk*blksiz*nchilds+blksiz*ichild+kii] = qrcv[irhs*ldqsend+jblk*blksiz+kii];
					};
				}
			};

		};

// Update QR decomposition

		UpdateQrdBlk (ilev, iblk, ldq, qblock);

// Send newly computed R data to the father cpu if necessary

		GetRPartQrd (_tree, ilev, _ibeg, _iend, qsendarr, ldqsend);

		if (fathercpuid != _tree.myid) {

			int msgtag = iblk+1;

			CMPIExchange::Send (pcomm, fathercpuid, msgtag,
										ldqsend*blksiz*sizeof(dcmplx), (char *)qsendarr);

		};

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			inode = fathernode;
			if (_tree.nodes[inode].nodecpu != _tree.myid) inode = -1;
		};

	};

// Free work arrays

	delete [] qblock;

};

// Author: Kharchenko S.A.
// Update QR decomposition of the complex block in parallel mode
//========================================================================================
void CQrdC::UpdateQrdRPart (const CTree &_tree, // Update QR decomposition of the block in parallel mode
						int _m, int _ibeg, int _iend, dcmplx *_a, int _lda, dcmplx *_tau) {

	const char *funcname = "UpdateQrdRPart";

	dcmplx czero (0.0e0,0.0e0);

// Check input parameters

	if (_iend-_ibeg != blksiz) throw " Wrong ibeg and iend parameters in UpdateQrd function";
	if (_ibeg%blksiz != 0) throw " Wrong iend parameter in UpdateQrd function";

	int iblk=_ibeg / blksiz;

// Get R part of the data and put it into send array

	int ldqsend = (iblk+1)*blksiz;

	if (nfiles == 0) {
//		flowvision::EqnSolver::GetRPartQrd (_m, _ibeg, _iend, _a, _lda, qsendarr, ldqsend);
		::GetRPartQrd (_m, _ibeg, _iend, _a, _lda, qsendarr, ldqsend);
	} else {
//		flowvision::EqnSolver::GetRPartQrd (_m, iblk, iblk+1, blksiz, 
		::GetRPartQrd (_m, iblk, iblk+1, blksiz, 
						files, bl2file, bsblk, 
						_lda, qsendarr, ldqsend);
	};

// Allocate work array

	dcmplx *qblock;

	qblock = new dcmplx [(iblk+1)*blksiz*blksiz*nchildsmax];
	if (!qblock) MemoryFail (funcname);

// Search the tree

	CMPIComm pcomm = _tree.GetComm ();

	int nodeend = _tree.cpuidend[_tree.myid];

	int inode = nodeend;

	while (inode >= 0) {

		int ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

		int fathernode = _tree.nodes[inode].fatherid;
		int fathercpuid = _tree.nodes[inode].fathercpu;

// Receive data from the child nodes and init factorization data

		int ldq = (iblk+1)*blksiz*nchilds;

		for (int kii=0;kii<ldq*blksiz;kii++) qblock[kii] = czero;

		for (int ichild=0;ichild<nchilds;ichild++) {

			int childnode = _tree.nodes[inode].childs[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			dcmplx *qrcv;

			if (childcpuid == _tree.myid) {
				qrcv = qsendarr;
			} else {

				CMPIStatus status;

				int msgtag = iblk+1;

				CMPIExchange::Recv (pcomm, childcpuid, msgtag,
											ldqsend*blksiz*sizeof(dcmplx), (char *)qrecvarr, status);

				qrcv = qrecvarr;

			};

			for (int jblk=0;jblk<=iblk;jblk++) {
				for (int irhs=0;irhs<blksiz;irhs++) {
					for (int kii=0;kii<blksiz;kii++) {
						qblock[irhs*ldq+jblk*blksiz*nchilds+blksiz*ichild+kii] = qrcv[irhs*ldqsend+jblk*blksiz+kii];
					};
				}
			};

		};

// Update QR decomposition

		UpdateQrdBlk (ilev, iblk, ldq, qblock);

// Send newly computed R data to the father cpu if necessary

		GetRPartQrd (_tree, ilev, _ibeg, _iend, qsendarr, ldqsend);

		if (fathercpuid != _tree.myid) {

			int msgtag = iblk+1;

			CMPIExchange::Send (pcomm, fathercpuid, msgtag,
										ldqsend*blksiz*sizeof(dcmplx), (char *)qsendarr);

		};

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			inode = fathernode;
			if (_tree.nodes[inode].nodecpu != _tree.myid) inode = -1;
		};

	};

// Free work arrays

	delete [] qblock;

};

// Author: Kharchenko S.A.
// Multiply Q factor by the current block in the parallel mode
//========================================================================================
void CQrdC::MvmQBlk (const CTree &_tree, // Multiply Q factor by the current block in the parallel mode
						int _m, int _n, int _nrhs,
						dcmplx *_q, int _ldq, dcmplx *_tau, 
						dcmplx *_x, int _ldx, dcmplx *_qx, int _ldqx) { 

	const char *funcname = "MvmQBlk";

	dcmplx czero (0.0e0,0.0e0);

// Check input parameters

	if (_n%blksiz != 0) throw " Wrong n parameter in MvmQBlk function";

	if (_nrhs != blksiz) throw " Wrong nrhs parameter in MvmQBlk function";

	int iblk=_n / blksiz - 1;

// Allocate work array

	dcmplx *qblock;

	qblock = new dcmplx [(iblk+1)*blksiz*blksiz*nchildsmax];
	if (!qblock) MemoryFail (funcname);

// Assign initial qrecvarr

	int ldqrecv = (iblk+1)*blksiz;

	for (int irhs=0;irhs<blksiz;irhs++) {
		for (int kii=0;kii<_n;kii++) {
			qrecvarr[irhs*ldqrecv+kii] = _x[irhs*_ldx+kii];
		};
	};

// Search the tree

	int nodebeg = _tree.cpuidbeg[_tree.myid];

	CMPIComm pcomm = _tree.GetComm ();

	int inode = nodebeg;

	while (inode >= 0) {

		int ilev = _tree.nodes[inode].nodelv;

		int nchilds = _tree.nodes[inode].nchilds;

//		int fathernode = _tree.nodes[inode].fatherid;
		int fathercpuid = _tree.nodes[inode].fathercpu;

// Receive data from father cpu if necessary

		if (fathercpuid != _tree.myid) {

			CMPIStatus status;

			int msgtag = iblk+1;

			CMPIExchange::Recv (pcomm, fathercpuid, msgtag,
										ldqrecv*blksiz*sizeof(dcmplx), (char *)qrecvarr, status);

		};

// Multiply

		int ldq = (iblk+1)*blksiz*nchilds;

		MvmQBlk (ilev, iblk, _nrhs, qrecvarr, ldqrecv , qblock, ldq);

// Init multiplication data and send them to the child nodes

		for (int ichild=0;ichild<nchilds;ichild++) {

			int childnode = _tree.nodes[inode].childs[ichild];
			int childcpuid = _tree.nodes[childnode].nodecpu;

			dcmplx *qsnd;

			if (childcpuid == _tree.myid) {
				qsnd = qrecvarr;
			} else {
				qsnd = qsendarr;
			};

			int ldqsend = (iblk+1)*blksiz;

			for (int jblk=0;jblk<=iblk;jblk++) {
				for (int irhs=0;irhs<blksiz;irhs++) {
					for (int kii=0;kii<blksiz;kii++) {
						qsnd[irhs*ldqsend+jblk*blksiz+kii] = qblock[irhs*ldq+jblk*blksiz*nchilds+blksiz*ichild+kii];
					};
				}
			};

			if (childcpuid != _tree.myid) {

				int msgtag = iblk+1;

				CMPIExchange::Send (pcomm, childcpuid, msgtag,
											ldqsend*blksiz*sizeof(dcmplx), (char *)qsendarr);

			};

		};

// Go down level

		if (nchilds == 1) {
			inode = -1;
		} else {
			int inodecurr = inode;
			inode = -1;
			for (int ichild=0;ichild<nchilds;ichild++) {
				int childnode = _tree.nodes[inodecurr].childs[ichild];
				if (_tree.nodes[childnode].nodecpu == _tree.myid) {
					inode = childnode;
					break;
				};
			};
		};

	};

// Perform multiplication of the local part of the decomposition

	if (nfiles == 0) {
		::MvmQBlk (_m, _n, _nrhs, _q, _ldq, _tau, 
					qrecvarr, ldqrecv, _qx, _ldqx);
	} else {

		int nblk = _n / _nrhs;

		::MvmQBlk (_m, nblk, blksiz, _nrhs,
					files, bl2file, bsblk, _ldq, _tau,
					qrecvarr, ldqrecv, _qx, _ldqx,
					_q);
	};

// Free work arrays

	delete [] qblock;

};

// Author: Kharchenko S.A.
// Add sparse integer arrays with no common entries
//========================================================================================
void AddSparseIArray (int _myid, int _nproc, int _nelem, int *_iarr) { // Add sparse integer arrays with no common entries

	const char *funcname = "AddSparseIArray";

// Count the local nomber of entries

	int i, nzloc=0;

	for (i=0;i<_nelem;i++) {
		if (_iarr[i] != 0) {
			nzloc++;
			if (_iarr[i] != 1) throw " Error in integer array elements: != 1";
		};
	};

// Exchange the numbers of elements 

	int *nzarr;
	int *nzsum;

	nzarr = new int [_nproc];
	if (!nzarr) MemoryFail (funcname);
	nzsum = new int [_nproc+1];
	if (!nzsum) MemoryFail (funcname);

	for (i=0;i<_nproc;i++) nzarr[i] = 0;
	nzarr[_myid] = nzloc;
	nzsum[0] = 0;

	CMPIComm commloc;

	CMPIExchange::ExchangeArrayMPI (commloc, INTEGERVALUE, ADD,
												_nproc, nzarr, nzsum+1);

// Prepare exchange arrays

	for (i=0;i<_nproc;i++) nzsum[i+1] += nzsum[i];

	int isize = nzsum[_nproc];
	int *iarrsnd;
	int *iarrrcv;

	iarrsnd = new int [isize];
	if (!iarrsnd) MemoryFail (funcname);
	iarrrcv = new int [isize];
	if (!iarrrcv) MemoryFail (funcname);

	for (i=0;i<isize;i++) iarrsnd[i] = 0;

	nzloc = nzsum[_myid];

	for (i=0;i<_nelem;i++) {
		if (_iarr[i] != 0) {
			iarrsnd[nzloc] = i;
			nzloc++;
		};
	};

// Exchange the data

	CMPIExchange::ExchangeArrayMPI (commloc, INTEGERVALUE, ADD,
												isize, iarrsnd, iarrrcv);

// Setup the final array

	for (i=0;i<_nelem;i++) _iarr[i] = 0;
	for (i=0;i<isize;i++) {
		if (_iarr[iarrrcv[i]] != 0) throw " Error in integer array elements: element assigned more than once";
		_iarr[iarrrcv[i]] = 1;
	};

// Free memory

	delete [] nzarr;
	delete [] nzsum;
	delete [] iarrsnd;
	delete [] iarrrcv;

};

// Author: Kharchenko S.A.
// Init timing tracing
//========================================================================================
void InitTracing (double *_times) { // Init timing tracing

	_times[0] = CMPIExchange::WallTimeMPI();
//	_times[0] = (double) clock () / (double) CLOCKS_PER_SEC;
//	_times[0] = (double) clock ();

};

// Author: Kharchenko S.A.
// Finalize timing tracing
//========================================================================================
void FinTracing (ofstream &_fout, char *_fname, double *_times) { // Finalize timing tracing

	double wtime1 = CMPIExchange::WallTimeMPI();
//	double wtime1 = (double) clock () / (double) CLOCKS_PER_SEC;
//	double wtime1 = (double) clock ();
	_fout << _fname << ": wall time = " << wtime1-_times[0] << std::endl;

};
