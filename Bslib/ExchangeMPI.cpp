//#include "stdafx.h"
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include "mpi.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "ExchangeMPI.h"
#include "globals.h"


//#define __NO_MPIFILE__

#ifdef _MPICH2
#ifdef SEEK_SET
   #undef SEEK_SET
#endif
#ifdef SEEK_CUR
   #undef SEEK_CUR
#endif
#ifdef SEEK_END
   #undef SEEK_END
#endif
#endif

using namespace std;

class CMPIComm MPIContext;

// using namespace EqnSolver;

// Author: Kharchenko S.A.
// Description: Zero length constructor
// CMPIComm::CMPIComm()
//========================================================================================
CMPIComm::CMPIComm () { // Zero length constructor

	comm = (void *) new MPI_Comm;
	*((MPI_Comm*)comm) = MPI_COMM_WORLD;

	ref = 0;

//	MPI_Comm_size (*((MPI_Comm*)comm),&procNumber);
//	MPI_Comm_rank (*((MPI_Comm*)comm),&procIndex);
	procNumber = 1;
	procIndex = 0;
	procGroup = 0;

};

// Author: Kharchenko S.A.
// Description: Destructor
// CMPIComm::~CMPIComm()
//========================================================================================
CMPIComm::~CMPIComm () { // Destructor
	MPI_Comm * pcomm = (MPI_Comm*) comm;

	if (ref != 0) {
		ref->DecreaseRef ();
		if (ref->GetCount () == 0) {
			MPI_Comm_free (pcomm);
			delete ref;
		};
	};
	delete pcomm;
};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CMPIComm::CMPIComm()
//========================================================================================
CMPIComm::CMPIComm (const CMPIComm &_comm) { // Copy constructor

	procIndex = _comm.procIndex;
	procGroup = _comm.procGroup;
	procNumber = _comm.procNumber;

	comm = new MPI_Comm;
	*((MPI_Comm*)comm) = *((MPI_Comm*)_comm.comm);

	if (_comm.ref == 0) {
		ref = 0;
	} else {
		ref = _comm.ref;
		ref->AddRef ();
	};

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CMPIComm::operator=()
//========================================================================================
CMPIComm &CMPIComm::operator= (const CMPIComm &_comm2) { // Equality operator

	MPI_Comm * pcomm = (MPI_Comm*) comm;

	if (ref != 0) {
		ref->DecreaseRef ();
		if (ref->GetCount () == 0) {
			MPI_Comm_free (pcomm);
			delete ref;
		};
	};

	delete pcomm;

	pcomm = new MPI_Comm;

	comm = (void *)pcomm;

	procIndex = _comm2.procIndex;
	procGroup = _comm2.procGroup;
	procNumber = _comm2.procNumber;

	*((MPI_Comm*)comm) = *((MPI_Comm*)_comm2.comm);

	if (_comm2.ref == 0) {
		ref = 0;
	} else {
		ref = _comm2.ref;
		ref->AddRef ();
	};

	return *this;

};

// Author: Kharchenko S.A.
// Description: Get the MPI communicator
// CMPIComm::GetCommMPI()
//========================================================================================
void * CMPIComm::GetCommMPI () { // Get the MPI communicator
	return comm;
};

// Author: Kharchenko S.A.
// Description: Set the MPI communicator
// CMPIComm::SetCommMPI_base()
//========================================================================================
void CMPIComm::SetCommMPI_base (void *_comm) { // Set the MPI communicator

	MPI_Comm * pcomm = (MPI_Comm*) comm;

	if (ref != 0) {
		ref->DecreaseRef ();
		if (ref->GetCount () == 0) {
			MPI_Comm_free (pcomm);
			delete ref;
		};
	};

	delete pcomm;
	comm = new MPI_Comm;
	*((MPI_Comm*)comm) = *((MPI_Comm*)_comm);

	ref = 0;

	procGroup = -1;
	procNumber = -1;
	procIndex = -1;

};

// Author: Kharchenko S.A.
// Description: Set the MPI communicator
// CMPIComm::SetCommMPI()
//========================================================================================
void CMPIComm::SetCommMPI (void *_comm) { // Set the MPI communicator

	SetCommMPI_base (_comm);

	MPI_Comm * pcomm = (MPI_Comm*) comm;

	MPI_Comm_size (*((MPI_Comm*)comm),&procNumber);
	MPI_Comm_rank (*((MPI_Comm*)comm),&procIndex);

	procGroup = 0;
//	cout << " In SetCommMpi: nproc = " << procNumber << " Myid = " << procIndex << endl;

};

// Author: Kharchenko S.A.
// Description: Set the group number
// CMPIComm::SetGroup()
//========================================================================================
void CMPIComm::SetGroup (int _igroup) { // Set group number

	procGroup = _igroup;

};

// Author: Kharchenko S.A.
// Description: Create communicator for a sublist of processors
// CMPIComm::CreateComm()
//========================================================================================
CMPIComm CMPIComm::CreateComm (int _nproc, int *_listcpu) { // Create communicator for a sublist of processors

	MPI_Comm * pcomm = (MPI_Comm*) comm;
	MPI_Group group;

	int info;

	info = MPI_Comm_group (*pcomm, &group);

	MPI_Group groupnew;

	info = MPI_Group_incl (group, _nproc, _listcpu, &groupnew);

	MPI_Comm commnew;

	info = MPI_Comm_create (*pcomm, groupnew, &commnew);

	info = MPI_Group_free (&group);
	info = MPI_Group_free (&groupnew);

	CMPIComm ccommnew;

	bool isCurrentCpu = false;

	int i;

	for (i=0;i<_nproc;i++) {
		if (_listcpu[i] == procIndex) isCurrentCpu = true;
	};

	if (isCurrentCpu) {
		ccommnew.SetCommMPI (&commnew);
//		ccommnew.ref = 0;
		ccommnew.ref = new CMPICommRef;
		ccommnew.ref->AddRef ();
	} else {
		ccommnew.SetCommMPI_base (&commnew);
	};

	return ccommnew;

};

// Author: Kharchenko S.A.
// Description: Create communicator for an initial sublist of processors
// CMPIComm::CreateCommInitial()
//========================================================================================
CMPIComm CMPIComm::CreateCommInitial (int _nproc) { // Create communicator for an initial sublist of processors

	int *listcpu;

	listcpu = new int [_nproc];

	int i;

	for (i=0;i<_nproc;i++) listcpu[i] = i;

	CMPIComm commnew = CreateComm (_nproc, listcpu);

	delete [] listcpu;

	return commnew;

};
/*
// Author: Kharchenko S.A.
// Description: Output communicator
// CMPIComm::>>()
//========================================================================================
ostream &operator<< (ostream &_stream, const CMPIComm &_comm) { // Output communicator

	_stream << " CMPIComm:" << endl;

	_stream << " procNumber = " << _comm.procNumber << endl;
	_stream << " procIndex = " << _comm.procIndex << endl;
	_stream << " procGroup = " << _comm.procGroup << endl;
	_stream << " MpiComm = " << *((MPI_Comm*)_comm.comm) << endl;

	return _stream;

};
*/
// Author: Kharchenko S.A.
// Description: Default constructor
// CMPIRequest::CMPIRequest()
//========================================================================================
CMPIRequest::CMPIRequest () { // Default constructor

	req = new MPI_Request;

};

// Author: Kharchenko S.A.
// Description: Destructor
// CMPIRequest::~CMPIRequest()
//========================================================================================
CMPIRequest::~CMPIRequest () { // Destructor

	delete (MPI_Request*) req;

};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CMPIRequest::CMPIRequest()
//========================================================================================
CMPIRequest::CMPIRequest (const CMPIRequest &_req) { // Copy constructor

	req = new MPI_Request;
	*((MPI_Request*)req) = *((MPI_Request*)_req.req);

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CMPIRequest::operator=()
//========================================================================================
CMPIRequest &CMPIRequest::operator= (const CMPIRequest &_req2) { // Equality operator

	delete (MPI_Request*) req;

	req = new MPI_Request;
	*((MPI_Request*)req) = *((MPI_Request*)_req2.req);

	return *this;

};

// Author: Kharchenko S.A.
// Description: Get the request handler pointer
// CMPIRequest::GetRequest()
//========================================================================================
void * CMPIRequest::GetRequest () { // Get the request handler pointer

	return req;

};

// Author: Kharchenko S.A.
// Description: Set the request handler pointer
// CMPIRequest::SetRequest()
//========================================================================================
void CMPIRequest::SetRequest (void *_req) { // Set the request handler pointer

	MPI_Request * preq = (MPI_Request*) req;
	delete preq;
	req = new MPI_Request;
	*((MPI_Request*)req) = *((MPI_Request*)_req);

};
/*
// Author: Kharchenko S.A.
// Description: Output request handler
// CMPIRequest::<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CMPIRequest &_req) { // Output request handler

	_stream << " CMPIRequest:" << endl;

	_stream << " req = " << *((MPI_Request*)_req.req) << endl;

	return _stream;

};
*/
// Author: Kharchenko S.A.
// Description: Default constructor
// CMPIRequest::CMPIRequest()
//========================================================================================
CMPIStatus::CMPIStatus () { // Default constructor

	stat = new MPI_Status;

};

// Author: Kharchenko S.A.
// Description: Destructor
// CMPIStatus::~CMPIStatus()
//========================================================================================
CMPIStatus::~CMPIStatus () { // Destructor

	delete (MPI_Status*) stat;

};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CMPIStatus::CMPIStatus()
//========================================================================================
CMPIStatus::CMPIStatus (const CMPIStatus &_stat) { // Copy constructor

	stat = new MPI_Status;
	*((MPI_Status*)stat) = *((MPI_Status*)_stat.stat);

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CMPIStatus::operator=()
//========================================================================================
CMPIStatus &CMPIStatus::operator= (const CMPIStatus &_stat2) { // Equality operator

	delete (MPI_Status*) stat;

	stat = new MPI_Status;
	*((MPI_Status*)stat) = *((MPI_Status*)_stat2.stat);

	return *this;

};

// Author: Kharchenko S.A.
// Description: Get the status pointer
// CMPIStatus::GetStatus()
//========================================================================================
void * CMPIStatus::GetStatus () { // Get the status handler pointer

	return stat;

};

// Author: Kharchenko S.A.
// Description: Set the status handler pointer
// CMPIStatus::SetStatus()
//========================================================================================
void CMPIStatus::SetStatus (void *_stat) { // Set the status handler pointer

	MPI_Status * pstat = (MPI_Status*) stat;
	delete pstat;
	stat = new MPI_Status;
	*((MPI_Status*)stat) = *((MPI_Status*)_stat);

};
/*
// Author: Kharchenko S.A.
// Description: Output status handler
// CMPIStatus::<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CMPIStatus &_stat) { // Output status handler

	_stream << " CMPIStatus:" << endl;

//	_stream << " stat = " << *((MPI_Status*)_stat.stat) << endl;

	return _stream;

};
*/
// Author: Kharchenko S.A.
// Description: Default constructor
// CMPIFile::CMPIFile()
//========================================================================================
CMPIFile::CMPIFile () { // Default constructor

	file = 0;

};

// Author: Kharchenko S.A.
// Description: Destructor
// CMPIFile::~CMPIFile()
//========================================================================================
CMPIFile::~CMPIFile () { // Destructor
};

// Author: Kharchenko S.A.
// Description: Copy constructor
// CMPIFile::CMPIFile()
//========================================================================================
CMPIFile::CMPIFile (const CMPIFile &_file) { // Copy constructor

	file = _file.file;

};

// Author: Kharchenko S.A.
// Description: Equality operator
// CMPIFile::operator=()
//========================================================================================
CMPIFile &CMPIFile::operator= (const CMPIFile &_file2) { // Equality operator

	file = _file2.file;

	return *this;

};

// Author: Kharchenko S.A.
// Description: Get the file handler pointer
// CMPIFile::GetFile()
//========================================================================================
void * CMPIFile::GetFile () { // Get the file handler pointer
	return file;
};

// Author: Kharchenko S.A.
// Description: Set the file handler pointer
// CMPIFile::SetFile()
//========================================================================================
void CMPIFile::SetFile (void *_file) { // Set the file handler pointer

	file = _file;

};

// Author: Kharchenko S.A.
// Description: Output file handler
// CMPIFile::<<()
//========================================================================================
ostream &operator<< (ostream &_stream, const CMPIFile &_file) { // Output file handler

	_stream << " CMPIFile:" << endl;

#ifndef __NO_MPIFILE__
	_stream << " MpiFile = " << *((MPI_File*)_file.file) << endl;
#endif

	return _stream;

};

// Author: Kharchenko S.A.
// Description: Initialize exchanges on all CPUs
// CMPIExchange::InitializeMPI()
//========================================================================================
int CMPIExchange::InitializeMPI (int* argc, char **argv[]) { // Initialize exchanges on all CPUs
	MPI_Init (argc,argv);
	return 0;
};

// Author: Kharchenko S.A.
// Description: Finalize exchanges on all CPUs
// CMPIExchange::FinalizeMPI()
//========================================================================================
int CMPIExchange::FinalizeMPI () { // Finalize exchanges on all CPUs
	MPI_Finalize();
	return 0;
};

// Author: Kharchenko S.A.
// Description: Get the standard MPI communicator
// CMPIExchange::GetCommMPI()
//========================================================================================
CMPIComm CMPIExchange::GetCommMPI () { // Get the standard MPI communicator
	CMPIComm commloc;
	MPI_Comm * pcomm = (MPI_Comm*) commloc.comm;
	delete pcomm;
	commloc.comm = new MPI_Comm;
	*((MPI_Comm*)commloc.comm) = MPI_COMM_WORLD;
	MPI_Comm_size (*((MPI_Comm*)commloc.comm),&commloc.procNumber);
	MPI_Comm_rank (*((MPI_Comm*)commloc.comm),&commloc.procIndex);
	commloc.procGroup = 0;

	return commloc;
};

// Author: Kharchenko S.A.
// Description: Get the group based MPI communicator
// CMPIExchange::GetCommMPI()
//========================================================================================
CMPIComm CMPIExchange::GetCommMPI (int _ngroups, int &_igroup) { // Get the group based MPI communicator

	CMPIComm commloc;
	MPI_Comm * pcomm = (MPI_Comm*) commloc.comm;
	delete pcomm;
	commloc.comm = new MPI_Comm;

	int ngroupsloc = _ngroups;
	if (ngroupsloc <= 1) ngroupsloc = 1;

	int nprocloc;
	int myidloc;

	MPI_Comm_size (MPI_COMM_WORLD,&nprocloc);
	MPI_Comm_rank (MPI_COMM_WORLD,&myidloc);

	_igroup = myidloc % ngroupsloc;
	int ikey = myidloc;

	MPI_Comm_split (MPI_COMM_WORLD, _igroup, ikey, (MPI_Comm*)commloc.comm);

	MPI_Comm_size (*((MPI_Comm*)commloc.comm),&commloc.procNumber);
	MPI_Comm_rank (*((MPI_Comm*)commloc.comm),&commloc.procIndex);

	commloc.procGroup = 0;

	return commloc;
};

// Author: Kharchenko S.A.
// Description: Free MPI communicator
// CMPIExchange::FreeCommMPI()
//========================================================================================
void CMPIExchange::FreeCommMPI (CMPIComm &_comm) { // Free MPI communicator
	MPI_Comm_free ((MPI_Comm*)_comm.comm);
};

// Author: Kharchenko S.A.
// Description: Get the number of CPUs
// CMPIExchange::GetNprocMPI()
//========================================================================================
int CMPIExchange::GetNprocMPI (CMPIComm &_comm) { // Get the number of CPUs
	int nproc;
	MPI_Comm_size (*((MPI_Comm*)_comm.comm),&nproc);
	return nproc;
};

// Author: Kharchenko S.A.
// Description: Get the ID number of the CPU
// CMPIExchange::GetMyidMPI()
//========================================================================================
int CMPIExchange::GetMyidMPI (CMPIComm &_comm) { // Get the ID number of the CPU
	int myid = 0;
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm),&myid);
	return myid;
};

// Author: Kharchenko S.A.
// Description: Get the global number of CPUs
// CMPIExchange::GetNprocMPIGlobal()
//========================================================================================
int CMPIExchange::GetNprocMPIGlobal () { // Get the global number of CPUs
	int nproc = 0;
	MPI_Comm_size (MPI_COMM_WORLD,&nproc);
	return nproc;
};

// Author: Kharchenko S.A.
// Description: Get the global ID number of the CPU
// CMPIExchange::GetMyidMPIGlobal()
//========================================================================================
int CMPIExchange::GetMyidMPIGlobal () { // Get the global ID number of the CPU
	int myid = 0;
	MPI_Comm_rank (MPI_COMM_WORLD,&myid);
	return myid;
};

// Author: Kharchenko S.A.
// Description: Get the processor name
// CMPIExchange::GetCPUNameMPI()
//========================================================================================
char * CMPIExchange::GetCPUNameMPI () { // Get the processor name
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name (processor_name,&namelen);
	char *name;
	name = new char [namelen+1];
	strcpy (name, processor_name);
	return name;
};

// Author: Kharchenko S.A.
// Description: Barrier synchronization of the CPUs
// CMPIExchange::SynchronizeMPI()
//========================================================================================
int CMPIExchange::SynchronizeMPI (const CMPIComm &_comm) { // Barrier synchronization of the CPUs
	MPI_Barrier (*((MPI_Comm*)_comm.comm));
	return 0;
};

// Author: Kharchenko S.A.
// Description: Open file
// CMPIExchange::MPIFileOpen()
//========================================================================================
int CMPIExchange::MPIFileOpen (CMPIComm &_comm, char* _filename, int _mode, CMPIFile &_File) { // Open file

#ifndef __NO_MPIFILE__

	MPI_Comm *pcomm = (MPI_Comm *) _comm.GetCommMPI ();
	MPI_File *pfile = new MPI_File;

// Assign open mode

	int amode = 0;

	if ((_mode & MPIFILE_READ) != 0)
	{
		if ((_mode & MPIFILE_WRITE) != 0)
			amode = MPI_MODE_RDWR | MPI_MODE_CREATE;
		else
			amode = MPI_MODE_RDONLY;
	}
	else
	{
		if ((_mode & MPIFILE_WRITE) != 0) amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
	}

	if ((_mode & MPIFILE_APPEND) != 0)
	{
		if ((_mode & MPIFILE_WRITE) != 0) amode |= MPI_MODE_APPEND;
	}

// Open the file

	int info;

	info = MPI_File_open (*pcomm, _filename, amode, MPI_INFO_NULL, pfile);
	if (info != 0) return info;

//	info = MPI_File_set_view (*pfile, 0, MPI_CHAR, MPI_CHAR, "Bytes", MPI_INFO_NULL);
//	if (info != 0) return info;

	_File.SetFile ((void *)pfile);

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description: Close file
// CMPIExchange::MPIFileClose()
//========================================================================================
int CMPIExchange::MPIFileClose (CMPIFile &_file) { // Close file

#ifndef __NO_MPIFILE__

	MPI_File *pfile = (MPI_File *)_file.GetFile ();

	int info;

	info = MPI_File_close (pfile);
	if (info != 0) return info;

	delete pfile;

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description: Write into the file
// CMPIExchange::MPIFileWrite()
//========================================================================================
int CMPIExchange::MPIFileWrite  (CMPIFile &_file, char* _buffer, size_t _size) { // Write into the file

#ifndef __NO_MPIFILE__

	MPI_File *pfile = (MPI_File *)_file.GetFile ();

	int info;
	MPI_Status status;

	info = MPI_File_write_shared (*pfile, _buffer, (int)_size, MPI_CHAR, &status);
	if (info != 0) return info;

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description: Read from the file
// CMPIExchange::MPIFileRead()
//========================================================================================
int CMPIExchange::MPIFileRead  (CMPIFile &_file, char* _buffer, size_t _size) { // Read from the file

#ifndef __NO_MPIFILE__

	MPI_File *pfile = (MPI_File *)_file.GetFile ();

	int info;
	MPI_Status status;

	info = MPI_File_read_shared (*pfile, _buffer, (int)_size, MPI_CHAR, &status);
	if (info != 0) return info;

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description: Position cursor of the file
// CMPIExchange::MPIFileSeek()
//========================================================================================
int CMPIExchange::MPIFileSeek  (CMPIFile &_file, int _offset, EMPISeekOrigin _fromWhere) { // Position cursor of the file

#ifndef __NO_MPIFILE__

	MPI_File *pfile = (MPI_File *)_file.GetFile ();

	int info;

	int whence;

	switch (_fromWhere)
	{
		case MPISEEK_BEGIN:
			whence = MPI_SEEK_SET;
		case MPISEEK_CURRENT:
			whence = MPI_SEEK_CUR;
		case MPISEEK_END:
			whence = MPI_SEEK_END;
	};

	info = MPI_File_seek_shared (*pfile, _offset, whence);
	if (info != 0) return info;

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description: Position cursor of the file
// CMPIExchange::MPIFileSeek()
//========================================================================================
int CMPIExchange::MPIFileSetSize  (CMPIFile &_file, off_t _offset) { // Trancate the file up to _offset

#ifndef __NO_MPIFILE__

	MPI_File *pfile = (MPI_File *)_file.GetFile ();

	int info;

	info = MPI_File_set_size (*pfile, _offset);
	if (info != 0) return info;

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description: Position cursor of the file
// CMPIExchange::MPIFileSeek()
//========================================================================================
int CMPIExchange::MPIFileDelete (char* _filename) { // Delete file

#ifndef __NO_MPIFILE__

// Delete the file

	int info;

	info = MPI_File_delete (_filename, MPI_INFO_NULL);
	if (info != 0) return info;

#endif

	return 0;

};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::CharAdd()
//========================================================================================
void CharAdd(void *in, void *inout, int *len, MPI_Datatype *dptr )
{
	char *inL = (char*)in;
	char *inoutL = (char*)inout;
	int i; 
	char c; 

	for (i=0; i< *len; ++i) { 
		c = *inoutL+*inL;
		*inoutL = c; 
		inL++; inoutL++; 
	};
};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::CharMaximum()
//========================================================================================
void CharMaximum(void *in, void *inout, int *len, MPI_Datatype *dptr )
{
	char *inL = (char*)in;
	char *inoutL = (char*)inout;
	int i; 
	char c; 

	for (i=0; i< *len; ++i) { 
		c = *inoutL>*inL ? *inoutL : *inL;
		*inoutL = c; 
		inL++; inoutL++; 
	};
};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::CharMinimum()
//========================================================================================
void CharMinimum(void *in, void *inout, int *len, MPI_Datatype *dptr )
{
	char *inL = (char*)in;
	char *inoutL = (char*)inout;
	int i; 
	char c; 

	for (i=0; i< *len; ++i) { 
		c = *inoutL<*inL ? *inoutL : *inL;
		*inoutL = c; 
		inL++; inoutL++; 
	};
};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::UCharAdd()
//========================================================================================
void UCharAdd(void *in, void *inout, int *len, MPI_Datatype *dptr )
{
	unsigned char *inL = (unsigned char*)in;
	unsigned char *inoutL = (unsigned char*)inout;
	int i; 
	unsigned char c; 

	for (i=0; i< *len; ++i) { 
		c = *inoutL+*inL;
		*inoutL = c; 
		inL++; inoutL++; 
	};
};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::UCharMaximum()
//========================================================================================
void UCharMaximum(void *in, void *inout, int *len, MPI_Datatype *dptr )
{
	unsigned char *inL = (unsigned char*)in;
	unsigned char *inoutL = (unsigned char*)inout;
	int i; 
	unsigned char c; 

	for (i=0; i< *len; ++i) { 
		c = *inoutL>*inL ? *inoutL : *inL;
		*inoutL = c; 
		inL++; inoutL++; 
	};
};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::UCharMinimum()
//========================================================================================
void UCharMinimum(void *in, void *inout, int *len, MPI_Datatype *dptr )
{
	unsigned char *inL = (unsigned char*)in;
	unsigned char *inoutL = (unsigned char*)inout;
	int i; 
	unsigned char c; 

	for (i=0; i< *len; ++i) { 
		c = *inoutL<*inL ? *inoutL : *inL;
		*inoutL = c; 
		inL++; inoutL++; 
	};
};

// Author: Kharchenko S.A.
// Description: Exchange array of small length among all processor with prescribed operation on all elements
// CMPIExchange::ExchangeArrayMPI()
//========================================================================================
int CMPIExchange::ExchangeArrayMPI(const CMPIComm &_comm, 
									MPIExchangeDataType _Datatype, MPIExchangeOperation _Operation,
									int _Length, void *_ArrayIni, void *_ArrayFin) 
{
	MPI_Datatype datatype;
	MPI_Op op;
	int nprocloc = _comm.GetNproc ();
	int size = 0;
	switch (_Datatype)
	{
	case CHARVALUE:
		datatype = MPI_CHAR;
		size = sizeof(char);
		switch (_Operation)
		{
		case UNDEFIN:
			return 0;
		case ADD:
//       op = MPI_SUM;
			MPI_Op_create (CharAdd, true, &op);
			break;
		case MAXIMUM:
//       op = MPI_MAX;
			MPI_Op_create (CharMaximum, true, &op);
			break;
		case MINIMUM:
//         op = MPI_MIN;
			MPI_Op_create (CharMinimum, true, &op);
			break;
		default:
			throw " CExchange::ExchangeArray: unknown operation";
		}
		break;
	case UNSIGNEDCHARVALUE:
		datatype = MPI_UNSIGNED_CHAR;
		size = sizeof(unsigned char);
		switch (_Operation)
		{
		case UNDEFIN:
			return 0;
		case ADD:
//       op = MPI_SUM;
			MPI_Op_create (UCharAdd, true, &op);
			break;
		case MAXIMUM:
//       op = MPI_MAX;
			MPI_Op_create (UCharMaximum, true, &op);
			break;
		case MINIMUM:
//       op = MPI_MIN;
			MPI_Op_create (UCharMinimum, true, &op);
			break;
		default:
			throw " CExchange::ExchangeArray: unknown operation";
		}
		break;
	case INTEGERVALUE:
		datatype = MPI_INT;
		size = sizeof(int);
		break;
	case UNSIGNEDINTEGERVALUE:
		datatype = MPI_UNSIGNED;
		size = sizeof(unsigned int);
		break;
	case FLOATVALUE:
		datatype = MPI_FLOAT;
		size = sizeof(float);
		break;
	case DOUBLEVALUE:
		datatype = MPI_DOUBLE;
		size = sizeof(double);
		break;
	case SHORTVALUE:
		datatype = MPI_SHORT;
		size = sizeof(short);
		break;
	case UNSIGNEDSHORTVALUE:
		datatype = MPI_UNSIGNED_SHORT;
		size = sizeof(unsigned short);
		break;
	default:
		throw " CExchange::ExchangeArray: unknown data type";
	}
	if (_Datatype != CHARVALUE && _Datatype != UNSIGNEDCHARVALUE) {
		switch (_Operation)
		{
		case UNDEFIN:
			return 0;
		case ADD:
			op = MPI_SUM;
			break;
		case MAXIMUM:
			op = MPI_MAX;
			break;
		case MINIMUM:
			op = MPI_MIN;
			break;
		default:
			throw " CExchange::ExchangeArray: unknown operation";
		}
	}
	if (_ArrayIni != _ArrayFin) 
	{
		if (nprocloc == 1) {
			memcpy(_ArrayFin, _ArrayIni, _Length*size);
		} else {
			MPI_Allreduce(_ArrayIni, _ArrayFin, _Length, datatype, op, *((MPI_Comm*)_comm.comm));
		};
	}
	else
	{
		if (nprocloc != 1) {
			void* buffer = new char[_Length*size];
			MPI_Allreduce(_ArrayIni, buffer, _Length, datatype, op, *((MPI_Comm*)_comm.comm));
			memcpy(_ArrayIni, buffer, _Length*size);
			delete [] (char *) buffer;
		};
	}
	return 0;
};

// Author: Kharchenko S.A.
// Description: Exchange array of small length among all processor with prescribed operation on all elements
// CMPIExchange::ExchangeArrayMPIListCpu()
//========================================================================================
int CMPIExchange::ExchangeArrayMPIListCpu(const CMPIComm &_comm, int _NListCpu, int *_ListCpu,
												MPIExchangeDataType _Datatype, MPIExchangeOperation _Operation,
												int _Length, void *_ArrayIni, void *_ArrayFin) 
{

	const char *funcname = "ExchangeArrayMPIListCpu";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	MPI_Comm_size (*((MPI_Comm*)_comm.comm), &nproc);
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm), &myid);

// Check parameters on entry

	if (_Datatype != INTEGERVALUE || _Operation != ADD) {
		throw " CMPIExchange::ExchangeArrayMPIListCpu: wrong params on entry ";
		return 1;
	};

// Allocate the destination array

	int isize = _NListCpu * _Length;
	int *iarr;

	iarr = new int [isize];
	if (!iarr) MemoryFail (funcname);

// Allocate and init sends and receives

// Count the total number of sends/receives

	int nsends = _NListCpu-1;
	int nrecvs = _NListCpu-1;

	MPI_Request *sndreqs, *rcvreqs;
	MPI_Status *sndstat, *rcvstat;

	sndreqs = new MPI_Request [nsends+1];
	if (!sndreqs) MemoryFail (funcname);
	rcvreqs = new MPI_Request [nrecvs+1];
	if (!rcvreqs) MemoryFail (funcname);
	sndstat = new MPI_Status  [nsends+1];
	if (!sndstat) MemoryFail (funcname);
	rcvstat = new MPI_Status  [nrecvs+1];
	if (!rcvstat) MemoryFail (funcname);

// Init sends

	int ilist, iproc, i, j;

	int msgtag, ierr;

	nsends = 0;

	for (ilist=0;ilist<_NListCpu;ilist++) {
		iproc = _ListCpu[ilist];

		if (iproc != myid) {

			msgtag = myid;

			ierr = MPI_Isend (_ArrayIni, _Length, MPI_INT, iproc, msgtag, *((MPI_Comm*)_comm.comm), sndreqs+nsends);

			nsends++;

		};
	};

// Init receives

	nrecvs = 0;

	int *piarr;

	for (ilist=0;ilist<_NListCpu;ilist++) {
		iproc = _ListCpu[ilist];

		if (iproc != myid) {

			msgtag = iproc;

			ierr = MPI_Irecv (iarr+_Length*ilist, _Length, MPI_INT, iproc, msgtag, *((MPI_Comm*)_comm.comm), rcvreqs+nrecvs);

			nrecvs++;
		} else {
			piarr = (int *) _ArrayIni;
			for (i=0;i<_Length;i++) iarr[_Length*ilist+i] = piarr[i];
		};

	};

// Wait for completion of sends/receives

	ierr = MPI_Waitall (nrecvs, rcvreqs, rcvstat);

	ierr = MPI_Waitall (nsends, sndreqs, sndstat);

// Compute the result

	int *resarr;

	resarr = new int [_Length];
	if (!resarr) MemoryFail (funcname);

	for (i=0;i<_Length;i++) resarr[i] = 0;

	for (i=0;i<_NListCpu;i++) {
		for (j=0;j<_Length;j++) {
			resarr[j] += iarr[i*_Length+j];
		};
	};

// Store result

	piarr = (int *) _ArrayFin;
	for (i=0;i<_Length;i++) piarr[i] = resarr[i];

// Free work arrays

	delete [] iarr;
	delete [] sndreqs;
	delete [] rcvreqs;
	delete [] sndstat;
	delete [] rcvstat;
	delete [] resarr;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Exchange data requests
// CMPIExchange::DataExchangeRequestMPI()
//========================================================================================
int CMPIExchange::DataExchangeRequestMPI (CMPIComm &_comm, // Exchange data requests
										int _NObjRecv, int *_ObjTypeRecv, int *_ObjIDRecv, int *_CpuIDRecv,
										int &_NObjSend, int *&_ObjTypeSend, int *&_ObjIDSend, int *&_CpuIDSend) {

	const char *funcname = "DataExchangeRequest";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	MPI_Comm_size (*((MPI_Comm*)_comm.comm), &nproc);
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm), &myid);

// Collect the numbers of objects from all cpu's

	int *cpu2obj;

	cpu2obj = new int [nproc+1];
	if (!cpu2obj) MemoryFail (funcname);

	int i, iproc;

	for (i=0;i<=nproc;i++) cpu2obj[i] = 0;

	cpu2obj[myid+1] = _NObjRecv;

	if (nproc > 1) {

		int *cpu2objloc;

		cpu2objloc = new int [nproc+1];
		if (!cpu2objloc) MemoryFail (funcname);

		MPI_Allreduce (cpu2obj, cpu2objloc, nproc+1,
						MPI_INT, MPI_SUM, *((MPI_Comm*)_comm.comm));

		for (i=0;i<=nproc;i++) cpu2obj[i] = cpu2objloc[i];

		delete [] cpu2objloc;

	};

	for (i=0;i<nproc;i++) cpu2obj[i+1] = cpu2obj[i] + cpu2obj[i+1];

// Allocate and combine all the data into one collected set of arrays

	int nobjtot = cpu2obj[nproc];

	int *objtypetot;
	int *objidtot;
	int *cpuidtot;

	objtypetot = new int [nobjtot];
	if (!objtypetot) MemoryFail (funcname);
	objidtot = new int [nobjtot];
	if (!objidtot) MemoryFail (funcname);
	cpuidtot = new int [nobjtot];
	if (!cpuidtot) MemoryFail (funcname);

	for (i=0;i<nobjtot;i++) {

		objtypetot[i] = 0;
		objidtot  [i] = 0;
		cpuidtot  [i] = 0;

	};

	int ibeg = cpu2obj[myid];

	for (i=cpu2obj[myid];i<cpu2obj[myid+1];i++) {

		objtypetot[i] = _ObjTypeRecv[i-ibeg];
		objidtot  [i] = _ObjIDRecv  [i-ibeg];
		cpuidtot  [i] = _CpuIDRecv  [i-ibeg];

	};

	if (nproc > 1) {

		int *listini, *listfin;

		listini = new int [3*nobjtot];
		if (!listini) MemoryFail (funcname);
		listfin = new int [3*nobjtot];
		if (!listfin) MemoryFail (funcname);

		for (i=0;i<nobjtot;i++) {

			listini[          i] = objtypetot [i];
			listini[  nobjtot+i] = objidtot   [i];
			listini[2*nobjtot+i] = cpuidtot   [i];

		};

		MPI_Allreduce (listini, listfin, 3*nobjtot,
						MPI_INT, MPI_SUM, *((MPI_Comm*)_comm.comm));

		for (i=0;i<nobjtot;i++) {

			objtypetot[i] = listfin[          i];
			objidtot  [i] = listfin[  nobjtot+i];
			cpuidtot  [i] = listfin[2*nobjtot+i];

		};

		delete [] listini;
		delete [] listfin;

	};

// Prepare the final data

	_NObjSend = 0;

	for (i=0;i<nobjtot;i++) {
		if (cpuidtot[i] == myid) _NObjSend++;
	};

	_ObjTypeSend = new int [_NObjSend];
	if (!_ObjTypeSend) MemoryFail (funcname);
	_ObjIDSend = new int [_NObjSend];
	if (!_ObjIDSend) MemoryFail (funcname);
	_CpuIDSend = new int [_NObjSend];
	if (!_CpuIDSend) MemoryFail (funcname);

	_NObjSend = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		for (i=cpu2obj[iproc];i<cpu2obj[iproc+1];i++) {
			if (cpuidtot[i] == myid) {
				_ObjTypeSend[_NObjSend] = objtypetot[i];
				_ObjIDSend[_NObjSend] = objidtot[i];
				_CpuIDSend[_NObjSend] = iproc;
				_NObjSend++;
			};
		};
	};

// Free work arrays

	delete [] cpu2obj;
	delete [] objtypetot;
	delete [] objidtot;
	delete [] cpuidtot;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Exchange sizes data requests
// CMPIExchange::SizesDataExchangeRequest()
//========================================================================================
int CMPIExchange::SizesDataExchangeRequest (CMPIComm &_comm, // Exchange sizes data requests
											int _NObjRecv, int *_ObjTypeRecv, int *_ObjIDRecv, int *_ObjSizeRecv, int *_CpuIDRecv,
											int &_NObjSend, int *&_ObjTypeSend, int *&_ObjIDSend, int *&_ObjSizeSend, int *&_CpuIDSend) {

	const char *funcname = "SizesDataExchangeRequest";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	MPI_Comm_size (*((MPI_Comm*)_comm.comm), &nproc);
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm), &myid);

// Collect the numbers of objects from all cpu's

	int *cpu2obj;

	cpu2obj = new int [nproc+1];
	if (!cpu2obj) MemoryFail (funcname);

	int i, iproc;

	for (i=0;i<=nproc;i++) cpu2obj[i] = 0;

	cpu2obj[myid+1] = _NObjRecv;

	if (nproc > 1) {

		int *cpu2objloc;

		cpu2objloc = new int [nproc+1];
		if (!cpu2objloc) MemoryFail (funcname);

		MPI_Allreduce (cpu2obj, cpu2objloc, nproc+1,
						MPI_INT, MPI_SUM, *((MPI_Comm*)_comm.comm));

		for (i=0;i<=nproc;i++) cpu2obj[i] = cpu2objloc[i];

		delete [] cpu2objloc;

	};

	for (i=0;i<nproc;i++) cpu2obj[i+1] = cpu2obj[i] + cpu2obj[i+1];

// Allocate and combine all the data into one collected set of arrays

	int nobjtot = cpu2obj[nproc];

	int *objtypetot;
	int *objidtot;
	int *objsizetot;
	int *cpuidtot;

	objtypetot = new int [nobjtot];
	if (!objtypetot) MemoryFail (funcname);
	objidtot = new int [nobjtot];
	if (!objidtot) MemoryFail (funcname);
	objsizetot = new int [nobjtot];
	if (!objsizetot) MemoryFail (funcname);
	cpuidtot = new int [nobjtot];
	if (!cpuidtot) MemoryFail (funcname);

	for (i=0;i<nobjtot;i++) {

		objtypetot[i] = 0;
		objidtot  [i] = 0;
		cpuidtot  [i] = 0;
		objsizetot[i] = 0;

	};

	int ibeg = cpu2obj[myid];

	for (i=cpu2obj[myid];i<cpu2obj[myid+1];i++) {

		objtypetot[i] = _ObjTypeRecv[i-ibeg];
		objidtot  [i] = _ObjIDRecv  [i-ibeg];
		cpuidtot  [i] = _CpuIDRecv  [i-ibeg];
		objsizetot[i] = _ObjSizeRecv[i-ibeg];

	};

	if (nproc > 1) {

		int *listini, *listfin;

		listini = new int [4*nobjtot];
		if (!listini) MemoryFail (funcname);
		listfin = new int [4*nobjtot];
		if (!listfin) MemoryFail (funcname);

		for (i=0;i<nobjtot;i++) {

			listini[          i] = objtypetot [i];
			listini[  nobjtot+i] = objidtot   [i];
			listini[2*nobjtot+i] = cpuidtot   [i];
			listini[3*nobjtot+i] = objsizetot [i];

		};

		MPI_Allreduce (listini, listfin, 4*nobjtot,
						MPI_INT, MPI_SUM, *((MPI_Comm*)_comm.comm));

		for (i=0;i<nobjtot;i++) {

			objtypetot[i] = listfin[          i];
			objidtot  [i] = listfin[  nobjtot+i];
			cpuidtot  [i] = listfin[2*nobjtot+i];
			objsizetot[i] = listfin[3*nobjtot+i];

		};

		delete [] listini;
		delete [] listfin;

	};

// Prepare the final data

	_NObjSend = 0;

	for (i=0;i<nobjtot;i++) {
		if (cpuidtot[i] == myid) _NObjSend++;
	};

	_ObjTypeSend = new int [_NObjSend];
	if (!_ObjTypeSend) MemoryFail (funcname);
	_ObjIDSend = new int [_NObjSend];
	if (!_ObjIDSend) MemoryFail (funcname);
	_ObjSizeSend = new int [_NObjSend];
	if (!_ObjSizeSend) MemoryFail (funcname);
	_CpuIDSend = new int [_NObjSend];
	if (!_CpuIDSend) MemoryFail (funcname);

	_NObjSend = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		for (i=cpu2obj[iproc];i<cpu2obj[iproc+1];i++) {
			if (cpuidtot[i] == myid) {
				_ObjTypeSend[_NObjSend] = objtypetot[i];
				_ObjIDSend[_NObjSend] = objidtot[i];
				_ObjSizeSend[_NObjSend] = objsizetot[i];
				_CpuIDSend[_NObjSend] = iproc;
				_NObjSend++;
			};
		};
	};

// Free work arrays

	delete [] cpu2obj;
	delete [] objtypetot;
	delete [] objidtot;
	delete [] objsizetot;
	delete [] cpuidtot;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Exchange data
// CMPIExchange::DataExchangeMPI()
//========================================================================================
int CMPIExchange::DataExchangeMPI (CMPIComm &_comm, // Exchange data
								int _NObjSend, int *_ObjTypeSend, int *_ObjIDSend, int *_CpuIDSend,
								int *_ObjSizeSend, char **_ObjSend,
								int &_NObjRecv, int *&_ObjTypeRecv, int *&_ObjIDRecv, int *&_CpuIDRecv,
								int *&_ObjSizeRecv, char **&_ObjRecv) {

	const char *funcname = "DataExchangeMPI";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	MPI_Comm_size (*((MPI_Comm*)_comm.comm), &nproc);
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm), &myid);

// Exchange send data and transform it into receive data

	int ierr;

	ierr = CMPIExchange::SizesDataExchangeRequest (_comm,
													_NObjSend, _ObjTypeSend, _ObjIDSend, _ObjSizeSend, _CpuIDSend,
													_NObjRecv, _ObjTypeRecv, _ObjIDRecv, _ObjSizeRecv, _CpuIDRecv);

// Allocate receive data

	int i, iproc;

	_ObjRecv = new char * [_NObjRecv];
	if (!_ObjRecv) MemoryFail (funcname);

	for (i=0;i<_NObjRecv;i++) {
		_ObjRecv[i] = new char [_ObjSizeRecv[i]];
		if (!_ObjRecv[i]) MemoryFail (funcname);
	};

// Create the set of send/receive MPI datatypes

	int *nobjsend, *nobjrecv;
	int *sizesendtot, *sizerecvtot;
	int *objsizes;
	char **objs;
	MPI_Datatype **mpisend, **mpirecv;

	int nloc = _NObjSend;
	if (_NObjRecv > nloc) nloc = _NObjRecv;

	objsizes = new int [nloc];
	if (!objsizes) MemoryFail (funcname);
	objs = new char * [nloc];
	if (!objs) MemoryFail (funcname);
	nobjsend = new int [nproc];
	if (!nobjsend) MemoryFail (funcname);
	nobjrecv = new int [nproc];
	if (!nobjrecv) MemoryFail (funcname);
	sizesendtot = new int [nproc];
	if (!sizesendtot) MemoryFail (funcname);
	sizerecvtot = new int [nproc];
	if (!sizerecvtot) MemoryFail (funcname);
	mpisend = new MPI_Datatype * [nproc];
	if (!mpisend) MemoryFail (funcname);
	mpirecv = new MPI_Datatype * [nproc];
	if (!mpirecv) MemoryFail (funcname);

	for (iproc=0;iproc<nproc;iproc++) {
		nloc = 0;
		sizesendtot[iproc] = 0;
		for (i=0;i<_NObjSend;i++) {
			if (_CpuIDSend[i] == iproc) {
				objsizes[nloc] = _ObjSizeSend[i];
				sizesendtot[iproc] += _ObjSizeSend[i];
				objs[nloc] = _ObjSend[i];
				nloc++;
			};
		};
		nobjsend[iproc] = nloc;
		mpisend[iproc] = (MPI_Datatype *) CMPIExchange::Char2Mpi (nloc, objsizes, objs);
		nloc = 0;
		sizerecvtot[iproc] = 0;
		for (i=0;i<_NObjRecv;i++) {
			if (_CpuIDRecv[i] == iproc) {
				objsizes[nloc] = _ObjSizeRecv[i];
				sizerecvtot[iproc] += _ObjSizeRecv[i];
				objs[nloc] = _ObjRecv[i];
				nloc++;
			};
		};
		nobjrecv[iproc] = nloc;
		mpirecv[iproc] = (MPI_Datatype *) CMPIExchange::Char2Mpi (nloc, objsizes, objs);
	};

// Perform local copy

	int ipsend, iprecv;

	if (nobjsend[myid] != 0) {
		if (nobjsend[myid] != nobjrecv[myid]) throw " DataExchange: different numbers of sends and receives to itself";
		ipsend = 0;
		iprecv = 0;
		nloc = 0;
		while (nloc != nobjsend[myid]) {
			while (_CpuIDSend[ipsend] != myid) ipsend++;
			while (_CpuIDRecv[iprecv] != myid) iprecv++;
			memcpy (_ObjRecv[iprecv], _ObjSend[ipsend], _ObjSizeSend[ipsend] * sizeof(char));
			ipsend++;
			iprecv++;
			nloc++;
		};
	};

// Count the total number of sends/receives

	int nsends = 0;
	int nrecvs = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (nobjsend[iproc] != 0) nsends++;
		if (nobjrecv[iproc] != 0) nrecvs++;
	};

	MPI_Request *sndreqs, *rcvreqs;
	MPI_Status *sndstat, *rcvstat;

	sndreqs = new MPI_Request [nsends+1];
	if (!sndreqs) MemoryFail (funcname);
	rcvreqs = new MPI_Request [nrecvs+1];
	if (!rcvreqs) MemoryFail (funcname);
	sndstat = new MPI_Status  [nsends+1];
	if (!sndstat) MemoryFail (funcname);
	rcvstat = new MPI_Status  [nrecvs+1];
	if (!rcvstat) MemoryFail (funcname);

	int msgtag;

// Init sends

	nsends = 0;

	for (iproc=0;iproc<nproc;iproc++) {

		if (iproc != myid && nobjsend[iproc] != 0 && sizesendtot[iproc] > 0) {

			msgtag = myid;

			ierr = MPI_Isend (MPI_BOTTOM, 1, *((MPI_Datatype*)mpisend[iproc]), iproc, msgtag, *((MPI_Comm*)_comm.comm), sndreqs+nsends);

			nsends++;

		};
	};

// Init receives

	nrecvs = 0;

	for (iproc=0;iproc<nproc;iproc++) {

		if (iproc != myid && nobjrecv[iproc] != 0 && sizerecvtot[iproc] > 0) {

			msgtag = iproc;

			ierr = MPI_Irecv (MPI_BOTTOM, 1, *((MPI_Datatype*)mpirecv[iproc]), iproc, msgtag, *((MPI_Comm*)_comm.comm), rcvreqs+nrecvs);

			nrecvs++;
		};

	};

// Wait for completion of sends/receives

	ierr = MPI_Waitall (nrecvs, rcvreqs, rcvstat);

	ierr = MPI_Waitall (nsends, sndreqs, sndstat);

// Free work arrays

	delete [] nobjsend;
	delete [] nobjrecv;
	delete [] sizesendtot;
	delete [] sizerecvtot;
	delete [] objsizes;
	delete [] objs;

	for (iproc=0;iproc<nproc;iproc++) {
		ierr = MPI_Type_free (mpisend[iproc]);
		if (ierr != MPI_SUCCESS) {
			cout << " MPI error occured in MPI_Type_free " << funcname << ": " << ierr << endl;
			throw " MPI error occured in MPI_Type_free ";
		};
		ierr = MPI_Type_free (mpirecv[iproc]);
		if (ierr != MPI_SUCCESS) {
			cout << " MPI error occured in MPI_Type_free " << funcname << ": " << ierr << endl;
			throw " MPI error occured in MPI_Type_free ";
		};
		delete mpisend[iproc];
		delete mpirecv[iproc];
	};

	delete [] mpisend;
	delete [] mpirecv;

	delete [] sndreqs;
	delete [] rcvreqs;
	delete [] sndstat;
	delete [] rcvstat;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Exchange sizes data requests
// CMPIExchange::SizesDataExchangeRequestListCpu()
//========================================================================================
int CMPIExchange::SizesDataExchangeRequestListCpu (CMPIComm &_comm, int _NListCpu, int *_ListCpu, // Exchange sizes data requests
																	int _NObjRecv, int *_ObjTypeRecv, int *_ObjIDRecv, int *_ObjSizeRecv, int *_CpuIDRecv,
																	int &_NObjSend, int *&_ObjTypeSend, int *&_ObjIDSend, int *&_ObjSizeSend, int *&_CpuIDSend) {

	const char *funcname = "SizesDataExchangeRequestListCpu";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	MPI_Comm_size (*((MPI_Comm*)_comm.comm), &nproc);
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm), &myid);

// Collect the numbers of objects from all cpu's

	int *cpu2obj;

	cpu2obj = new int [nproc+1];
	if (!cpu2obj) MemoryFail (funcname);

	int i, iproc;

	for (i=0;i<=nproc;i++) cpu2obj[i] = 0;

	cpu2obj[myid+1] = _NObjRecv;

	if (nproc > 1) {

		ExchangeArrayMPIListCpu (_comm, _NListCpu, _ListCpu,
											INTEGERVALUE, ADD,
											nproc+1, cpu2obj, cpu2obj);

	};

	for (i=0;i<nproc;i++) cpu2obj[i+1] = cpu2obj[i] + cpu2obj[i+1];

// Allocate and combine all the data into one collected set of arrays

	int nobjtot = cpu2obj[nproc];

	int *objtypetot;
	int *objidtot;
	int *objsizetot;
	int *cpuidtot;

	objtypetot = new int [nobjtot];
	if (!objtypetot) MemoryFail (funcname);
	objidtot = new int [nobjtot];
	if (!objidtot) MemoryFail (funcname);
	objsizetot = new int [nobjtot];
	if (!objsizetot) MemoryFail (funcname);
	cpuidtot = new int [nobjtot];
	if (!cpuidtot) MemoryFail (funcname);

	for (i=0;i<nobjtot;i++) {

		objtypetot[i] = 0;
		objidtot  [i] = 0;
		cpuidtot  [i] = 0;
		objsizetot[i] = 0;

	};

	int ibeg = cpu2obj[myid];

	for (i=cpu2obj[myid];i<cpu2obj[myid+1];i++) {

		objtypetot[i] = _ObjTypeRecv[i-ibeg];
		objidtot  [i] = _ObjIDRecv  [i-ibeg];
		cpuidtot  [i] = _CpuIDRecv  [i-ibeg];
		objsizetot[i] = _ObjSizeRecv[i-ibeg];

	};

	if (nproc > 1) {

		int *listini;

		listini = new int [4*nobjtot];
		if (!listini) MemoryFail (funcname);

		for (i=0;i<nobjtot;i++) {

			listini[          i] = objtypetot [i];
			listini[  nobjtot+i] = objidtot   [i];
			listini[2*nobjtot+i] = cpuidtot   [i];
			listini[3*nobjtot+i] = objsizetot [i];

		};

		ExchangeArrayMPIListCpu (_comm, _NListCpu, _ListCpu,
											INTEGERVALUE, ADD,
											4*nobjtot, listini, listini);

		for (i=0;i<nobjtot;i++) {

			objtypetot[i] = listini[          i];
			objidtot  [i] = listini[  nobjtot+i];
			cpuidtot  [i] = listini[2*nobjtot+i];
			objsizetot[i] = listini[3*nobjtot+i];

		};

		delete [] listini;

	};

// Prepare the final data

	_NObjSend = 0;

	for (i=0;i<nobjtot;i++) {
		if (cpuidtot[i] == myid) _NObjSend++;
	};

	_ObjTypeSend = new int [_NObjSend];
	if (!_ObjTypeSend) MemoryFail (funcname);
	_ObjIDSend = new int [_NObjSend];
	if (!_ObjIDSend) MemoryFail (funcname);
	_ObjSizeSend = new int [_NObjSend];
	if (!_ObjSizeSend) MemoryFail (funcname);
	_CpuIDSend = new int [_NObjSend];
	if (!_CpuIDSend) MemoryFail (funcname);

	_NObjSend = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		for (i=cpu2obj[iproc];i<cpu2obj[iproc+1];i++) {
			if (cpuidtot[i] == myid) {
				_ObjTypeSend[_NObjSend] = objtypetot[i];
				_ObjIDSend[_NObjSend] = objidtot[i];
				_ObjSizeSend[_NObjSend] = objsizetot[i];
				_CpuIDSend[_NObjSend] = iproc;
				_NObjSend++;
			};
		};
	};

// Free work arrays

	delete [] cpu2obj;
	delete [] objtypetot;
	delete [] objidtot;
	delete [] objsizetot;
	delete [] cpuidtot;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Exchange data
// CMPIExchange::DataExchangeMPIListCpu()
//========================================================================================
int CMPIExchange::DataExchangeMPIListCpu (CMPIComm &_comm, int _NListCpu, int *_ListCpu, // Exchange data
												int _NObjSend, int *_ObjTypeSend, int *_ObjIDSend, int *_CpuIDSend,
												int *_ObjSizeSend, char **_ObjSend,
												int &_NObjRecv, int *&_ObjTypeRecv, int *&_ObjIDRecv, int *&_CpuIDRecv,
												int *&_ObjSizeRecv, char **&_ObjRecv) {

	const char *funcname = "DataExchangeMPIListCpu";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	MPI_Comm_size (*((MPI_Comm*)_comm.comm), &nproc);
	MPI_Comm_rank (*((MPI_Comm*)_comm.comm), &myid);

// Check that all send data match the cpu list

	int i, icpu, ilist, iproc;
	bool is_cpu_in_list;

	for (i=0;i<_NObjSend;i++) {
		icpu = _CpuIDSend[i];
		is_cpu_in_list = false;
		for (ilist=0;ilist<_NListCpu;ilist++) {
			iproc = _ListCpu[ilist];
			if (icpu == iproc) is_cpu_in_list = true;
		};
		if (!is_cpu_in_list) throw " CMPIExchange::DataExchangeMPIListCpu: Error in the cpu list ";
	};

// Exchange send data and transform it into receive data

	int ierr;

	ierr = CMPIExchange::SizesDataExchangeRequestListCpu (_comm, _NListCpu, _ListCpu,
																			_NObjSend, _ObjTypeSend, _ObjIDSend, _ObjSizeSend, _CpuIDSend,
																			_NObjRecv, _ObjTypeRecv, _ObjIDRecv, _ObjSizeRecv, _CpuIDRecv);

// Allocate receive data

	_ObjRecv = new char * [_NObjRecv];
	if (!_ObjRecv) MemoryFail (funcname);

	for (i=0;i<_NObjRecv;i++) {
		_ObjRecv[i] = new char [_ObjSizeRecv[i]];
		if (!_ObjRecv[i]) MemoryFail (funcname);
	};

// Create the set of send/receive MPI datatypes

	int *nobjsend, *nobjrecv;
	int *sizesendtot, *sizerecvtot;
	int *objsizes;
	char **objs;
	MPI_Datatype **mpisend, **mpirecv;

	int nloc = _NObjSend;
	if (_NObjRecv > nloc) nloc = _NObjRecv;

	objsizes = new int [nloc];
	if (!objsizes) MemoryFail (funcname);
	objs = new char * [nloc];
	if (!objs) MemoryFail (funcname);
	nobjsend = new int [nproc];
	if (!nobjsend) MemoryFail (funcname);
	nobjrecv = new int [nproc];
	if (!nobjrecv) MemoryFail (funcname);
	sizesendtot = new int [nproc];
	if (!sizesendtot) MemoryFail (funcname);
	sizerecvtot = new int [nproc];
	if (!sizerecvtot) MemoryFail (funcname);
	mpisend = new MPI_Datatype * [nproc];
	if (!mpisend) MemoryFail (funcname);
	mpirecv = new MPI_Datatype * [nproc];
	if (!mpirecv) MemoryFail (funcname);

	for (iproc=0;iproc<nproc;iproc++) {
		nloc = 0;
		sizesendtot[iproc] = 0;
		for (i=0;i<_NObjSend;i++) {
			if (_CpuIDSend[i] == iproc) {
				objsizes[nloc] = _ObjSizeSend[i];
				sizesendtot[iproc] += _ObjSizeSend[i];
				objs[nloc] = _ObjSend[i];
				nloc++;
			};
		};
		nobjsend[iproc] = nloc;
		mpisend[iproc] = (MPI_Datatype *) CMPIExchange::Char2Mpi (nloc, objsizes, objs);
		nloc = 0;
		sizerecvtot[iproc] = 0;
		for (i=0;i<_NObjRecv;i++) {
			if (_CpuIDRecv[i] == iproc) {
				objsizes[nloc] = _ObjSizeRecv[i];
				sizerecvtot[iproc] += _ObjSizeRecv[i];
				objs[nloc] = _ObjRecv[i];
				nloc++;
			};
		};
		nobjrecv[iproc] = nloc;
		mpirecv[iproc] = (MPI_Datatype *) CMPIExchange::Char2Mpi (nloc, objsizes, objs);
	};

// Perform local copy

	int ipsend, iprecv;

	if (nobjsend[myid] != 0) {
		if (nobjsend[myid] != nobjrecv[myid]) throw " DataExchange: different numbers of sends and receives to itself";
		ipsend = 0;
		iprecv = 0;
		nloc = 0;
		while (nloc != nobjsend[myid]) {
			while (_CpuIDSend[ipsend] != myid) ipsend++;
			while (_CpuIDRecv[iprecv] != myid) iprecv++;
			memcpy (_ObjRecv[iprecv], _ObjSend[ipsend], _ObjSizeSend[ipsend] * sizeof(char));
			ipsend++;
			iprecv++;
			nloc++;
		};
	};

// Count the total number of sends/receives

	int nsends = 0;
	int nrecvs = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (nobjsend[iproc] != 0) nsends++;
		if (nobjrecv[iproc] != 0) nrecvs++;
	};

	MPI_Request *sndreqs, *rcvreqs;
	MPI_Status *sndstat, *rcvstat;

	sndreqs = new MPI_Request [nsends+1];
	if (!sndreqs) MemoryFail (funcname);
	rcvreqs = new MPI_Request [nrecvs+1];
	if (!rcvreqs) MemoryFail (funcname);
	sndstat = new MPI_Status  [nsends+1];
	if (!sndstat) MemoryFail (funcname);
	rcvstat = new MPI_Status  [nrecvs+1];
	if (!rcvstat) MemoryFail (funcname);

	int msgtag;

// Init sends

	nsends = 0;

	for (ilist=0;ilist<_NListCpu;ilist++) {
		iproc = _ListCpu[ilist];

		if (iproc != myid && nobjsend[iproc] != 0 && sizesendtot[iproc] > 0) {

			msgtag = myid;

			ierr = MPI_Isend (MPI_BOTTOM, 1, *((MPI_Datatype*)mpisend[iproc]), iproc, msgtag, *((MPI_Comm*)_comm.comm), sndreqs+nsends);

			nsends++;

		};
	};

// Init receives

	nrecvs = 0;

	for (ilist=0;ilist<_NListCpu;ilist++) {
		iproc = _ListCpu[ilist];

		if (iproc != myid && nobjrecv[iproc] != 0 && sizerecvtot[iproc] > 0) {

			msgtag = iproc;

			ierr = MPI_Irecv (MPI_BOTTOM, 1, *((MPI_Datatype*)mpirecv[iproc]), iproc, msgtag, *((MPI_Comm*)_comm.comm), rcvreqs+nrecvs);

			nrecvs++;
		};

	};

// Wait for completion of sends/receives

	ierr = MPI_Waitall (nrecvs, rcvreqs, rcvstat);

	ierr = MPI_Waitall (nsends, sndreqs, sndstat);

// Free work arrays

	delete [] nobjsend;
	delete [] nobjrecv;
	delete [] sizesendtot;
	delete [] sizerecvtot;
	delete [] objsizes;
	delete [] objs;

	for (iproc=0;iproc<nproc;iproc++) {
		ierr = MPI_Type_free (mpisend[iproc]);
		if (ierr != MPI_SUCCESS) {
			cout << " MPI error occured in MPI_Type_free " << funcname << ": " << ierr << endl;
			throw " MPI error occured in MPI_Type_free ";
		};
		ierr = MPI_Type_free (mpirecv[iproc]);
		if (ierr != MPI_SUCCESS) {
			cout << " MPI error occured in MPI_Type_free " << funcname << ": " << ierr << endl;
			throw " MPI error occured in MPI_Type_free ";
		};
		delete mpisend[iproc];
		delete mpirecv[iproc];
	};

	delete [] mpisend;
	delete [] mpirecv;

	delete [] sndreqs;
	delete [] rcvreqs;
	delete [] sndstat;
	delete [] rcvstat;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Send data to the other CPU
// CMPIExchange::Send()
//========================================================================================
int CMPIExchange::Send (CMPIComm &_comm, int _rank, int _msgtag, // Send data to the other CPU
								int _length, char *_buffer) {

//	const char *funcname = "Send";

	int ierr = MPI_Send (_buffer, _length, MPI_CHAR, _rank, _msgtag, *((MPI_Comm*)_comm.comm));

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Send:" << ierr << endl;
		throw " MPI error occured in MPI_Send ";
	};

	return 0;

};

// Author: Kharchenko S.A.
// Description: Asynchronous send data to the other CPU
// CMPIExchange::ISend()
//========================================================================================
int CMPIExchange::ISend (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous send data to the other CPU
									int _length, char *_buffer, CMPIRequest &_req) {

//	const char *funcname = "ISend";

	int ierr = MPI_Isend (_buffer, _length, MPI_CHAR, _rank, _msgtag, *((MPI_Comm*)_comm.comm), ((MPI_Request*)_req.req));

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_ISend:" << ierr << endl;
		throw " MPI error occured in MPI_ISend ";
	};

	return 0;

};

// Author: Kharchenko S.A.
// Description: Recv data from the other CPU
// CMPIExchange::Recv()
//========================================================================================
int CMPIExchange::Recv (CMPIComm &_comm, int _rank, int _msgtag, // Recv data from the other CPU
								int _length, char *_buffer, CMPIStatus &_stat) {

//	const char *funcname = "Recv";

	int ierr = MPI_Recv (_buffer, _length, MPI_CHAR, _rank, _msgtag, *((MPI_Comm*)_comm.comm), ((MPI_Status*)_stat.stat));

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Recv:" << ierr << endl;
		throw " MPI error occured in MPI_Recv ";
	};

	return 0;

};

// Author: Kharchenko S.A.
// Description: Asynchronous recv data from the other CPU
// CMPIExchange::IRecv()
//========================================================================================
int CMPIExchange::IRecv (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous recv data from the other CPU
								int _length, char *_buffer, CMPIRequest &_req) {

//	const char *funcname = "IRecv";

	int ierr = MPI_Irecv (_buffer, _length, MPI_CHAR, _rank, _msgtag, *((MPI_Comm*)_comm.comm), ((MPI_Request*)_req.req));

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_IRecv:" << ierr << endl;
		throw " MPI error occured in MPI_IRecv ";
	};

	return 0;

};

// Author: Kharchenko S.A.
// Description: Wait for completion of exchanges
// CMPIExchange::WaitAll()
//========================================================================================
int CMPIExchange::WaitAll  (int _count, CMPIRequest *_recvarr, CMPIStatus *_statarr) { // Wait for completion of exchanges

	const char *funcname = "WaitAll";

// Create local arrays

	MPI_Request *reqvarrloc;
	MPI_Status *statarrloc;

	reqvarrloc = new MPI_Request [_count];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new MPI_Status [_count];
	if (!statarrloc) MemoryFail (funcname);

// Copy initial data

	int i;

	for (i=0;i<_count;i++) {
		reqvarrloc[i] = *((MPI_Request *)_recvarr[i].req);
	};

// Wait

	int ierr = MPI_Waitall (_count, reqvarrloc, statarrloc);

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Waitall:" << ierr << endl;
		throw " MPI error occured in MPI_Waitall ";
	};

// Copy resulting data

	for (i=0;i<_count;i++) {
		*((MPI_Status *)_statarr[i].stat) = statarrloc[i];
	};

// Free work arrays

	delete [] reqvarrloc;
	delete [] statarrloc;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Wait for completion of any one exchange
// CMPIExchange::WaitAny()
//========================================================================================
int CMPIExchange::WaitAny  (int _count, CMPIRequest *_recvarr, int &_index, CMPIStatus &_stat) { // Wait for completion of any one exchange

	const char *funcname = "WaitAny";

// Create local arrays

	MPI_Request *reqvarrloc;
	MPI_Status statloc;

	reqvarrloc = new MPI_Request [_count];
	if (!reqvarrloc) MemoryFail (funcname);

// Copy initial data

	int i;

	for (i=0;i<_count;i++) {
		reqvarrloc[i] = *((MPI_Request *)_recvarr[i].req);
	};

// Wait

	int ierr = MPI_Waitany (_count, reqvarrloc, &_index, &statloc);

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Waitany:" << ierr << endl;
		throw " MPI error occured in MPI_Waitany ";
	};

// Copy resulting data

	*((MPI_Status *)_stat.stat) = statloc;

// Free work arrays

	delete [] reqvarrloc;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Waits for completion of some set of exchanges
// CMPIExchange::WaitSome()
//========================================================================================
int CMPIExchange::WaitSome  (int _count, CMPIRequest *_recvarr, int &_ncompl, int *_indarr, CMPIStatus *_statarr) { // Waits for completion of some set of exchanges

	const char *funcname = "WaitSome";

// Create local arrays

	MPI_Request *reqvarrloc;
	MPI_Status *statarrloc;

	reqvarrloc = new MPI_Request [_count];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new MPI_Status [_count];
	if (!statarrloc) MemoryFail (funcname);

// Copy initial data

	int i;

	for (i=0;i<_count;i++) {
		reqvarrloc[i] = *((MPI_Request *)_recvarr[i].req);
	};

// Wait

	int ierr = MPI_Waitsome (_count, reqvarrloc, &_ncompl, _indarr, statarrloc);

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Waitsome:" << ierr << endl;
		throw " MPI error occured in MPI_Waitsome ";
	};

// Copy resulting data

	for (i=0;i<_ncompl;i++) {
		*((MPI_Status *)_statarr[i].stat) = statarrloc[i];
	};

// Free work arrays

	delete [] reqvarrloc;
	delete [] statarrloc;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Test for completion of any exchange
// CMPIExchange::Test()
//========================================================================================
int CMPIExchange::Test  (CMPIRequest &_req, bool &_flag, CMPIStatus &_stat) { // Test for completion of any exchange

	int iflag;

	MPI_Test (((MPI_Request*)_req.req), &iflag, ((MPI_Status*)_stat.stat));

	if (iflag != 0) {
		_flag = true;
	} else {
		_flag = false;
	};

	return 0;

};

// Author: Kharchenko S.A.
// Description: Test for completion of any one exchange
// CMPIExchange::TestAny()
//========================================================================================
int CMPIExchange::TestAny  (int _count, CMPIRequest *_recvarr, bool &_flag, int &_index, CMPIStatus &_stat) { // Wait for completion of any one exchange

	const char *funcname = "TestAny";

// Create local arrays

	MPI_Request *reqvarrloc;
	MPI_Status statloc;

	reqvarrloc = new MPI_Request [_count];
	if (!reqvarrloc) MemoryFail (funcname);

// Copy initial data

	int i;

	for (i=0;i<_count;i++) {
		reqvarrloc[i] = *((MPI_Request *)_recvarr[i].req);
	};

// Test

	int iflag;

	int ierr = MPI_Testany (_count, reqvarrloc, &_index, &iflag, &statloc);

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Testany:" << ierr << endl;
		throw " MPI error occured in MPI_Testany ";
	};

// Copy resulting data

	_flag = false;

	if (iflag != 0) {
		_flag = true;
		*((MPI_Status *)_stat.stat) = statloc;
	};

// Free work arrays

	delete [] reqvarrloc;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Tests for completion of some set of exchanges
// CMPIExchange::TestSome()
//========================================================================================
int CMPIExchange::TestSome  (int _count, CMPIRequest *_recvarr, int &_ncompl, int *_indarr, CMPIStatus *_statarr) { // Test for completion of some set of exchanges

	const char *funcname = "TestSome";

// Create local arrays

	MPI_Request *reqvarrloc;
	MPI_Status *statarrloc;

	reqvarrloc = new MPI_Request [_count];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new MPI_Status [_count];
	if (!statarrloc) MemoryFail (funcname);

// Copy initial data

	int i;

	for (i=0;i<_count;i++) {
		reqvarrloc[i] = *((MPI_Request *)_recvarr[i].req);
	};

// Wait

	int ierr = MPI_Testsome (_count, reqvarrloc, &_ncompl, _indarr, statarrloc);

	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Testsome:" << ierr << endl;
		throw " MPI error occured in MPI_Testsome ";
	};

// Copy resulting data

	for (i=0;i<_ncompl;i++) {
		*((MPI_Status *)_statarr[i].stat) = statarrloc[i];
	};

// Free work arrays

	delete [] reqvarrloc;
	delete [] statarrloc;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Synchronous check if the expected message is received
// CMPIExchange::Probe()
//========================================================================================
int CMPIExchange::Probe (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous check if the expected message is received
									CMPIStatus &_stat) {

//	const char *funcname = "Probe";

	int rankloc = _rank;
	int msgtagloc = _msgtag;

	if (rankloc < 0) rankloc = MPI_ANY_SOURCE;
	if (msgtagloc < 0) msgtagloc = MPI_ANY_TAG;

	MPI_Probe (rankloc, msgtagloc, *((MPI_Comm*)_comm.comm), ((MPI_Status*)_stat.stat));

	return 0;

};

// Author: Kharchenko S.A.
// Description: Asynchronous check if the expected message is received
// CMPIExchange::IProbe()
//========================================================================================
int CMPIExchange::IProbe (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous check if the expected message is received
									bool &_flag, CMPIStatus &_stat) {

//	const char *funcname = "IProbe";

	int iflag;

	int rankloc = _rank;
	int msgtagloc = _msgtag;

	if (rankloc < 0) rankloc = MPI_ANY_SOURCE;
	if (msgtagloc < 0) msgtagloc = MPI_ANY_TAG;

	MPI_Iprobe (rankloc, msgtagloc, *((MPI_Comm*)_comm.comm), &iflag, ((MPI_Status*)_stat.stat));

	if (iflag != 0) {
		_flag = true;
	} else {
		_flag = false;
	};

	return 0;

};

// Author: Kharchenko S.A.
// Description: Returns the source of the message
// CMPIExchange::GetSource()
//========================================================================================
int CMPIExchange::GetSource (CMPIStatus &_stat, int &_source) { // Returns the source of the message

	MPI_Status *pstat = (MPI_Status*)_stat.stat;

	_source = pstat->MPI_SOURCE;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Returns the tag of the message
// CMPIExchange::GetTag()
//========================================================================================
int CMPIExchange::GetTag (CMPIStatus &_stat, int &_tag) { // Returns the tag of the message

	MPI_Status *pstat = (MPI_Status*)_stat.stat;

	_tag = pstat->MPI_TAG;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Returns the number of received elements in char's
// CMPIExchange::GetCount()
//========================================================================================
int CMPIExchange::GetCount (CMPIStatus &_stat, int &_count) { // Returns the number of received elements in char's

	MPI_Get_count (((MPI_Status*)_stat.stat), MPI_CHAR, &_count);

	return 0;

};

// Author: Kharchenko S.A.
// Description: Create MPI datatype for the set of character arrays
// CMPIExchange::Char2Mpi()
//========================================================================================
void * CMPIExchange::Char2Mpi (int _NObjs, int *_ObjSizes, char **_Objs) { // Create MPI datatype for the set of character arrays

	const char *funcname = "Char2Mpi";

// Init MPI arrays

	int *mpiblen;
	MPI_Aint *mpidisp;
	MPI_Datatype *mpitype;

	mpiblen = new int [_NObjs];
	if (!mpiblen) MemoryFail (funcname);
	mpidisp = new MPI_Aint [_NObjs];
	if (!mpidisp) MemoryFail (funcname);
	mpitype = new MPI_Datatype [_NObjs];
	if (!mpitype) MemoryFail (funcname);

	int i, ierr;

	for (i=0;i<_NObjs;i++) {
		mpitype[i] = MPI_CHAR;
		mpiblen[i] = _ObjSizes[i];

		ierr = MPI_Address (_Objs[i], &mpidisp[i]);
		if (ierr != MPI_SUCCESS) {
			cout << " MPI error occured in MPI_Address " << funcname << ": " << ierr << endl;
			throw " MPI error occured in MPI_Address ";
		};
	};

	MPI_Datatype *temp;

	temp = new MPI_Datatype;

	ierr = MPI_Type_struct (_NObjs, mpiblen, mpidisp, mpitype, temp);
	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Type_struct " << funcname << ": " << ierr << endl;
		throw " MPI error occured in MPI_Type_struct ";
	};

	ierr = MPI_Type_commit (temp);
	if (ierr != MPI_SUCCESS) {
		cout << " MPI error occured in MPI_Type_commit " << funcname << ": " << ierr << endl;
		throw " MPI error occured in MPI_Type_commit ";
	};

	delete [] mpiblen;
	delete [] mpidisp;
	delete [] mpitype;

	return temp;

};

// Author: Kharchenko S.A.
// Description: Create communicator for a group of processors
// CMPIExchange::CreateCommForGroup()
//========================================================================================
CMPIComm CMPIExchange::CreateCommForGroup (CMPIComm &_comm, int _NGroups) { // Create communicator for a group of processors

	int ngroupsloc = _NGroups;
	if (ngroupsloc <= 1) ngroupsloc = 1;

	int nprocloc;
	int myidloc;

	MPI_Comm_size (MPI_COMM_WORLD,&nprocloc);
	MPI_Comm_rank (MPI_COMM_WORLD,&myidloc);

	int ncpugroup = (nprocloc+ngroupsloc-1) / ngroupsloc;

	int imyid = myidloc % ncpugroup;
	int igroup = (myidloc-imyid) / ncpugroup;

	int ibeg = igroup * ncpugroup;
	int iend = (igroup+1) * ncpugroup-1;
	if (iend > nprocloc-1) iend = nprocloc-1;

	int ranges[][3] = {ibeg,iend,1};

	MPI_Group group, groupnew;

	MPI_Comm_group (*((MPI_Comm*)_comm.comm), &group);

	MPI_Group_range_incl (group, 1, ranges, &groupnew);

	CMPIComm commloc;
	MPI_Comm * pcomm = (MPI_Comm*) commloc.comm;
	delete pcomm;
	commloc.comm = new MPI_Comm;

	MPI_Comm_create (*((MPI_Comm*)_comm.comm), groupnew, (MPI_Comm*)commloc.comm);

	MPI_Comm_size (*((MPI_Comm*)commloc.comm),&commloc.procNumber);
	MPI_Comm_rank (*((MPI_Comm*)commloc.comm),&commloc.procIndex);

	commloc.procGroup = igroup;

	return commloc;

};

// Author: Kharchenko S.A.
// Description:
// CMPIExchange::WallTimeMPI()
//========================================================================================
double CMPIExchange::WallTimeMPI()
{
	return MPI_Wtime();
}

extern "C" {

void ParMETIS_V3_PartKway (int *vtxdist, int *xadj, int *adjncy, 
                           int *vwgt, int *adjwgt, int *wgtflag, 
                           int *numflag, int *ncon, int *nparts, 
                           float *tpwgts, float *ubvec,
                           int *options, int *edgecut, int *part, MPI_Comm *comm); 

};

// Author: Kharchenko S.A.
// Description: Interface for ParMETIS
// CParMetis::PartKway()
//========================================================================================
void CParMETIS::PartKway (CMPIComm &_comm, int *_vtxdist, int *_xadj, int *_adjncy, 
                          int *_vwgt, int *_adjwgt, int *_wgtflag, 
                          int *_numflag, int *_ncon, int *_nparts, 
                          float *_tpwgts, float *_ubvec,
                          int *_options, int *_edgecut, int *_part)
{
   void *pcomm = _comm.GetCommMPI ();
//   ParMETIS_V3_PartKway (_vtxdist, _xadj, _adjncy, _vwgt, _adjwgt, _wgtflag, 
//                         _numflag, _ncon, _nparts, _tpwgts, _ubvec,
//                         _options, _edgecut, _part, ((MPI_Comm *) pcomm)); 
}

#ifdef __FV_2004__
} // namespace flowvision
#endif
