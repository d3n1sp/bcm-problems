#include "mpi.h"

int MPI_Init(int *argc, char ***argv) {return 0;};

int MPI_Get_processor_name (char *name, int *resultlen) {
		char *nameloc = "Test\0";
		*name = *nameloc;
		*resultlen = 5;
		return 0;
	};

int MPI_Comm_size (MPI_Comm comm, int *size) {*size = 1; return 0;};

int MPI_Comm_rank (MPI_Comm comm, int *rank) {*rank = 0; return 0;};

int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm) {return 0;};

int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) {return 0;};

int MPI_Irecv (void *buf, int count, MPI_Datatype datatype, int source, 
				int tag, MPI_Comm comm, MPI_Request *request) {return 0;};

int MPI_Isend (void *buf, int count, MPI_Datatype datatype, int dest, int tag,
				MPI_Comm comm, MPI_Request *request) {return 0;};

int MPI_Waitany (int count, MPI_Request array_of_requests[], int *index, MPI_Status *status) {return 0;};

int MPI_Request_free (MPI_Request *request) {return 0;};

int MPI_Waitall (int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) {return 0;};

int MPI_Allreduce (void *snd , void *rcv, int count, 
					MPI_Datatype datatype, MPI_Op oper, MPI_Comm comm) {return 0;};

int MPI_Barrier (MPI_Comm comm) {return 0;};

int MPI_Type_commit(MPI_Datatype *) {return 0;};

int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *) {return 0;};

int MPI_Address(void*, MPI_Aint *) {return 0;};

int MPI_Finalize () {return 0;};

