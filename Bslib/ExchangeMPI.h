#include <iostream>

//#ifndef __FV_2004__
//#define __FV_2004__
//#endif

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

// mpiexchange.h: Description of MPI exchange module
//
//////////////////////////////////////////////////////////////////////////////

//#ifndef _MPICH2
//#define _MPICH2
//#endif

#ifndef __MPIExchange
#define __MPIExchange

// Preliminary declarations

enum MPIExchangeDataType {CHARVALUE=0, UNSIGNEDCHARVALUE, INTEGERVALUE, UNSIGNEDINTEGERVALUE, 
                          FLOATVALUE, DOUBLEVALUE, SHORTVALUE, UNSIGNEDSHORTVALUE};

// MPIExchangeDataType structure specifies the data type in the exchange operation:
//
// CHARVALUE - char
// UNSIGNEDCHAR - unsigned char
// INTEGERVALUE - int
// UNSIGNEDINTEGERVALUE - unsigned int
// FLOATVALUE   - float
// DOUBLEVALUE  - double
// SHORTVALUE   - short
// UNSIGNEDSHORTVALUE - unsigned short

enum MPIExchangeOperation {UNDEFIN=-1, ADD=0, MAXIMUM, MINIMUM};

// MPIExchangeOperation structure specifies the exchange operation to be performed:
//
// UNDEFIN  - undefined uperation, do nothing
// ADD      - addition of each value
// MAXIMUM  - maximum  of each value
// MINIMUM  - minimum  of each value

enum EMPIFileOpenModeFlags
{
   MPIFILE_READ = 0x01,
   MPIFILE_WRITE = 0x02,
   MPIFILE_APPEND = 0x04
};

// EMPIFileOpenModeFlags structure specifies the file open flags
//
// MPIFILE_READ   - file for reading
// MPIFILE_WRITE  - file for writing
// MPIFILE_APPEND - file is appended on open

enum EMPISeekOrigin
{
   MPISEEK_BEGIN,
   MPISEEK_CURRENT,
   MPISEEK_END
};

// EMPISeekOrigin structure specifies the seek operations
//
// MPISEEK_BEGIN   - go into the beginning of the file
// MPISEEK_CURRENT - go into the current position in the file
// MPISEEK_END     - go into the ending    of the file

class CMPICommRef
{
public:
	int count; // Reference counter
// Constructors and destructor
	CMPICommRef () {count = 0;}; // Default constructor
	~CMPICommRef () {}; // Destructor
// Interface functions
	void AddRef () {count++;};
	void DecreaseRef () {count--;};
	int GetCount () {return count;};
};

class CMPIComm
{
public:
	void *comm;
	CMPICommRef *ref;
	int procIndex;
	int procGroup;
	int procNumber;
// Functions
// Constructors and destructor
	CMPIComm (); // Default constructor
	CMPIComm (const CMPIComm &_aa); // Copy constructor
	~CMPIComm (); // Destructor
// Get/set functions
	void * GetCommMPI (); // Get the MPI communicator
	int GetGroup () const {return procGroup;}; // Get group
	int GetNproc () const {return procNumber;}; // Get nproc
	int GetMyid () const {return procIndex;}; // Get myid
	void SetCommMPI (void *_comm); // Set the MPI communicator
	void SetCommMPI_base (void *_comm); // Set the MPI communicator
	void SetGroup (int _igroup); // Set the group number
// Operator functions
	CMPIComm &operator= (const CMPIComm &_comm); // Equality operator
// CMPIComm functions
	CMPIComm CreateComm (int _nproc, int *_listcpu); // Create communicator for a sublist of processors
	CMPIComm CreateCommInitial (int _nproc); // Create communicator for an initial sublist of processors
// Friend classes
	friend class CMPIExchange;
//	friend std::ostream &operator<< (std::ostream &_stream, const CMPIComm &_comm); // Output communicator
};

class CMPIRequest
{
public:
	void *req;
// Functions
// Constructors and destructor
	CMPIRequest (); // Default constructor
	CMPIRequest (const CMPIRequest &_req); // Copy constructor
	~CMPIRequest (); // Destructor
// Get/set functions
	void * GetRequest (); // Get the request handler pointer
	void SetRequest (void *_req); // Set the request handler pointer
// Operator functions
	CMPIRequest &operator= (const CMPIRequest &_req); // Equality operator
// Friend classes
	friend class CMPIExchange;
//	friend std::ostream &operator<< (std::ostream &_stream, const CMPIRequest &_req); // Output request
};

class CMPIStatus
{
public:
	void *stat;
// Functions
// Constructors and destructor
	CMPIStatus (); // Default constructor
	CMPIStatus (const CMPIStatus &_stat); // Copy constructor
	~CMPIStatus (); // Destructor
// Get/set functions
	void * GetStatus (); // Get the request handler pointer
	void SetStatus (void *_stat); // Set the request handler pointer
// Operator functions
	CMPIStatus &operator= (const CMPIStatus &_stat); // Equality operator
// Friend classes
	friend class CMPIExchange;
//	friend std::ostream &operator<< (std::ostream &_stream, const CMPIStatus &_stat); // Output status
};

class CMPIFile
{
public:
	void *file;
// Functions
// Constructors and destructor
	CMPIFile (); // Default constructor
	CMPIFile (const CMPIFile &_aa); // Copy constructor
	~CMPIFile (); // Destructor
// Get/set functions
	void * GetFile (); // Get the file handler
	void SetFile (void *_file); // Set the FILE handler
// Operator functions
	CMPIFile &operator= (const CMPIFile &_file); // Equality operator
// Friend classes
	friend class CMPIExchange;
	friend std::ostream &operator<< (std::ostream &_stream, const CMPIFile &_file); // Output file
};

extern class CMPIComm MPIContext;

class CMPIExchange
{
// Data
	public:
// User functions
// Control functions
	static int InitializeMPI (int* argc, char **argv[]); // Initialize exchanges on all CPUs
	static int FinalizeMPI (); // Finalize exchanges on all CPUs
	static int SynchronizeMPI (const CMPIComm &_comm); // Barrier synchronization of the CPUs
	static CMPIComm CreateCommForGroup (CMPIComm &_comm, int _NGroups); // Create communicator for a group of processors
// Set & get functions
	static CMPIComm GetCommMPI (); // Get the standard MPI communicator
	static CMPIComm GetCommMPI (int _ngroups, int &_igroup); // Get the group based MPI communicator
	static void FreeCommMPI (CMPIComm &_comm); // Free MPI communicator
	static int GetNprocMPI (CMPIComm &_comm); // Get the number of CPUs
	static int GetMyidMPI (CMPIComm &_comm); // Get the ID number of the CPU
	static int GetNprocMPIGlobal (); // Get the global number of CPUs
	static int GetMyidMPIGlobal (); // Get the ID number of the CPU
	static char * GetCPUNameMPI (); // Get the processor name
// Exchange functions
	static int ExchangeArrayMPI (const CMPIComm &_comm, // Exchange array of small length among all processor with prescribed operation on all elements
								MPIExchangeDataType _Datatype, MPIExchangeOperation _Operation,
								int _Length, void *_ArrayIni, void *_ArrayFin);
	static int ExchangeArrayMPIListCpu (const CMPIComm &_comm, // Exchange array of small length among all processor with prescribed operation on all elements
													int _NListCpu, int *_ListCpu,
													MPIExchangeDataType _Datatype, MPIExchangeOperation _Operation,
													int _Length, void *_ArrayIni, void *_ArrayFin);
	static int DataExchangeRequestMPI (CMPIComm &_comm, // Exchange data requests
									int _NObjRecv, int *_ObjTypeRecv, int *_ObjIDRecv, int *_CpuIDRecv,
									int &_NObjSend, int *&_ObjTypeSend, int *&_ObjIDSend, int *&_CpuIDSend);
	static int SizesDataExchangeRequestListCpu (CMPIComm &_comm, int _NListCpu, int *_ListCpu, // Exchange sizes data requests
																int _NObjRecv, int *_ObjTypeRecv, int *_ObjIDRecv, int *_ObjSizeRecv, int *_CpuIDRecv,
																int &_NObjSend, int *&_ObjTypeSend, int *&_ObjIDSend, int *&_ObjSizeSend, int *&_CpuIDSend);
	static int DataExchangeMPI (CMPIComm &_comm, // Exchange data
											int _NObjSend, int *_ObjTypeSend, int *_ObjIDSend, int *_CpuIDSend,
											int *_ObjSizeSend, char **_ObjSend,
											int &_NObjRecv, int *&_ObjTypeRecv, int *&_ObjIDRecv, int *&_CpuIDRecv,
											int *&_ObjSizeRecv, char **&_ObjRecv);
	static int DataExchangeMPIListCpu (CMPIComm &_comm, int _NListCpu, int *_ListCpu, // Exchange data
													int _NObjSend, int *_ObjTypeSend, int *_ObjIDSend, int *_CpuIDSend,
													int *_ObjSizeSend, char **_ObjSend,
													int &_NObjRecv, int *&_ObjTypeRecv, int *&_ObjIDRecv, int *&_CpuIDRecv,
													int *&_ObjSizeRecv, char **&_ObjRecv);
// Send/recv/test/wait/probe functions
	static int Send (CMPIComm &_comm, int _rank, int _msgtag, // Send data to the other CPU
								int _length, char *_buffer);
	static int ISend (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous send data to the other CPU
								int _length, char *_buffer, CMPIRequest &_recv);
	static int Recv (CMPIComm &_comm, int _rank, int _msgtag, // Receive data from the other CPU
								int _length, char *_buffer, CMPIStatus &_stat);
	static int IRecv (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous receive data from the other CPU
								int _length, char *_buffer, CMPIRequest &_recv);
	static int WaitAll  (int _count, CMPIRequest *_recvarr, CMPIStatus *_statarr); // Wait for completion of exchanges
	static int WaitAny  (int _count, CMPIRequest *_recvarr, int &_index, CMPIStatus &_stat); // Wait for completion of any one exchange
	static int WaitSome (int _count, CMPIRequest *_recvarr, int &_ncompl, int *_indarr, CMPIStatus *_statarr); // Waits for completion of some set of exchanges
	static int Test     (CMPIRequest &_req, bool &_flag, CMPIStatus &_stat); // Test for completion of any exchange
	static int TestAny  (int _count, CMPIRequest *_recvarr, bool &_flag, int &_index, CMPIStatus &_stat); // Test for completion of any one exchange
	static int TestSome (int _count, CMPIRequest *_recvarr, int &_ncompl, int *_indarr, CMPIStatus *_statarr); // Tests for completion of some set of exchanges
	static int Probe    (CMPIComm &_comm, int _rank, int _msgtag, // Synchronous check if the expected message is received
								CMPIStatus &_stat);
	static int IProbe   (CMPIComm &_comm, int _rank, int _msgtag, // Asynchronous check if the expected message is received
								bool &_flag, CMPIStatus &_stat);
	static int GetSource (CMPIStatus &_stat, int &_source); // Returns the source of the message
	static int GetTag    (CMPIStatus &_stat, int &_tag);    // Returns the tag of the message
	static int GetCount  (CMPIStatus &_stat, int &_count);  // Returns the number of received elements in char's
// IO functions
	static int MPIFileOpen   (CMPIComm &_comm, char* _filename, int _mode, CMPIFile &_File); // Open file
	static int MPIFileClose  (CMPIFile &_handle); // Close file
	static int MPIFileWrite  (CMPIFile &_file, char* _buffer, size_t _size); // Write into the file
	static int MPIFileRead   (CMPIFile &_file, char* _buffer, size_t _size); // Read from the file
	static int MPIFileSeek   (CMPIFile &_file, int _offset, EMPISeekOrigin _fromWhere); // Position cursor of the file
	static int MPIFileSetSize(CMPIFile &_file, /*off_t*/long _offset); // Trancate the file up to _offset
	static int MPIFileDelete (char* _filename); // Delete file
// Miscelleneous functions
	static double WallTimeMPI();
private:
// Internal functions
	static int SizesDataExchangeRequest (CMPIComm &_comm, // Exchange sizes data requests
									int _NObjRecv, int *_ObjTypeRecv, int *_ObjIDRecv, int *_ObjSizeRecv, int *_CpuIDRecv,
									int &_NObjSend, int *&_ObjTypeSend, int *&_ObjIDSend, int *&_ObjSizeSend, int *&_CpuIDSend);
	static void * Char2Mpi (int _NObjs, int *_ObjSizes, char **_Objs); // Create MPI datatype for the set of character arrays
};

class CParMETIS
{
   public:
   static void PartKway (CMPIComm &_comm, int *_vtxdist, int *_xadj, int *_adjncy, 
                         int *_vwgt, int *_adjwgt, int *_wgtflag, 
                         int *_numflag, int *_ncon, int *_nparts, 
                         float *_tpwgts, float *_ubvec,
                         int *_options, int *_edgecut, int *_part);
};

#endif
