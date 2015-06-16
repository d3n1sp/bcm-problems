//------------------------------------------------------------------------------------------------
// File: Decomp.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include <iostream>

#include "globals.h"
#include "ExchangeMPI.h"

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

// decomp.h: Description of decomposition module
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Decomp
#define __Decomp

class CMPIComm;
class CSolverPar;
class CCpuTree;

// Preliminary declarations

class CSMatrix;

class CDecomp
{
// Data
public:
// User functions
	static int DistributeCellsOverCPUs (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
									int _Nproc, double *_ProcWeight,
									int _Nx, int _Ny, int _Nz, int _Collap,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU);
	static int DistributeCellsOverCPUsForSortedTree (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
									int _Nproc, double *_ProcWeight,
									int _Nx, int _Ny, int _Nz, int _Collap,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU);
	static int DistributeCellsOverCPUs (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU);
	static int DistributeCellsOverCPUsForSortedTree (CMPIComm &_comm, // Compute distributions of the cells over CPU's
											int _Nproc, double *_ProcWeight,
											int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
											CSolverPar &_Param,
											CCpuTree &_CpuTree, int *&_CellsOrder, 
											int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU);
	static int DistributeCellsOverCPUs (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU, int *&_HyperCellWeight);
	static int DistributeCellsOverCPUsForSortedTree (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU, int *&_HyperCellWeight);
	static int DistributeCellsOverCPUsNestedDessectionWeights (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU, int *&_HyperCellWeight);
	static int DistributeCellsOverCPUsNestedDessection (CMPIComm &_comm, // Compute distributions of the cells over CPU's
										int _Nproc, double *_ProcWeight,
										int _NCells, int *_IANeiCells, int *_JANeiCells,
										CSolverPar &_Param,
										CCpuTree &_CpuTree, int *&_CellsOrder, 
										int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU);
	static int DistributeCellsOverCPUsNestedDessectionParMetis (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																					int _N1LevelCells, int *_Index1LevelCells,
																					int *_IA1Level2Real, int *_Real2Index,
																					int *_IARealCells, int *_JAReal2Index,
																					CSolverPar &_Param, 
																					CCpuTree &_CpuTree,
																					int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																					int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
																					int *&_Index1LevelCellsNew, int *&_IA1Level2RealNew);
	static int DistributeCellsOverCPUsNestedDessectionParMetisWeights (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																					int _N1LevelCells, int *_Index1LevelCells,
																					int *_IA1Level2Real, int *_Real2Index,
																					int *_IARealCells, int *_JAReal2Index,
																					CSolverPar &_Param, 
																					CCpuTree &_CpuTree,
																					int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																					int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
																					int *&_Index1LevelCellsNew);
	static int DistributeCellsNDParMetisBoundary1Level (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																			int _N1LevelCells, int *_Index1LevelCells,
																			int *_IA1Level2Real, int *_Real2Index,
																			int *_IARealCells, int *_JAReal2Index,
																			CSolverPar &_Param, 
																			CCpuTree &_CpuTree,
																			int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																			int &_NHyperCellsReal, int *&_HyperCell2CPUReal, int *&_HyperCellsInRealCells,
																			int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells, int *&_IA1Level2RealNew,
																			int *&_Index1LevelCellsNew, int *&_IndexRealCellsNew, int *&_IndexRealCellsIni);
	static int DistributeCellsNDParMetisBoundary1LevelWeights (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																			int _N1LevelCells, int *_Index1LevelCells,
																			int *_IA1Level2Real, int *_Real2Index,
																			int *_IARealCells, int *_JAReal2Index,
																			CSolverPar &_Param, 
																			CCpuTree &_CpuTree,
																			int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																			int &_NHyperCellsReal, int *&_HyperCell2CPUReal, int *&_HyperCellsInRealCells,
																			int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
																			int *&_Index1LevelCellsNew, int *&_IndexRealCellsNew, int *&_IndexRealCellsIni);
	static int DistributeRegular3DGridWeights (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
													int _Nx, int _Ny, int _Nz,
													CSolverPar &_Param, 
													CCpuTree &_CpuTree,
													int &_NHyperCells, int *&_HyperCell2CPU,
													int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
													int *&_Index1LevelCellsNew);
	static int DistributeRegular3DGrid (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
													int _Nx, int _Ny, int _Nz,
													CSolverPar &_Param, 
													CCpuTree &_CpuTree,
													int &_NHyperCells, int *&_HyperCell2CPU,
													int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
													int *&_Index1LevelCellsNew);
	static int DistributeRegular3DGridDirect (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
													int _Nx, int _Ny, int _Nz,
													CSolverPar &_Param, 
													CCpuTree &_CpuTree,
													int &_NHyperCells, int *&_HyperCell2CPU,
													int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
													int *&_Index1LevelCellsNew);
private:
	static void PartitionIndices (CMPIComm &_comm, // Compute partitioning of the indices for balanced load
									int _nprocpart, int _ncells, int *_listindex, int *_sprnds,
									int *_blksindex);
	static CSMatrix TransformSparsity (CMPIComm &_comm, // Create CSMatrix 2index structure that is based on pair {1levelcell,realcell} instead of {hypercell,realcell1}
													int _N1LevelCells, int *_sprnds, int *_Index1LevelCells,
													int *_Real2Index,
													int *_IARealCells, int *_JAReal2Index);
	static void SpPartitioning (CMPIComm &_comm, // Compute supernodes ordering and blocks partitioning
											int _nlistsp, int *_listsp, int *_sprndssp,
											int *_ordernd, int _nblks, int *_blksnd, int *_blks2cpu,
											int *_ordersp, int *_blkssp);

private:
// Internal functions
};

class CDecompHCells
{
// Data
	int nhcells;    // nhcells is the number of hypercells
	int *hcells;    // hcells[nhcells+1] is IA type array that describes the number of cells in the hypercells
	char ***cellschars; // cellschars array contains the char descriptions of all cells
	int sizeofbox;      // sizeofbox is the sizeof () of the box structure
	char **cellsboxes;  // cellsboxes array contains the bounding boxes of all cells
public:
// Functions
// Constructors and destructor
	CDecompHCells (); // Zero data constructor
	CDecompHCells (int _nhcells, int *_hcells); // Memory allocation constructor with zero data
	~CDecompHCells (); // Destructor
// Operator functions
	CDecompHCells &operator= (const CDecompHCells &_objhcells); // Equality operator
// Get/set functions
	int GetNhcells () const {return nhcells;};  // Get nhcells
	int GetHcellSize (int _ihcell) const {return hcells[_ihcell+1]-hcells[_ihcell];};  // Get hcell size
	char ** GetCellschars (int _ihcell) const {return cellschars[_ihcell];}; // Get cellschars
	void * GetCellsboxes (int _ihcell) const {return (void *)cellsboxes[_ihcell];}; // Get cellsboxes
	void SetSizeofbox (int _sizeofbox) {sizeofbox = _sizeofbox;};  // Set sizeofbox
	void SetCellsboxes (int _ihcell, char *_ptr) {cellsboxes[_ihcell] = _ptr;}; // Set cellsboxes
// Init function
	void FillBoxesData (CMPIComm &_comm, int _nblks, int *_blks, // Fill the data that describe the boxes
								int _nlist, int *_list, void *_boxes, int _sizeofbox);
// Input and output
	friend std::ostream &operator<< (std::ostream &_stream, const CDecompHCells &_obj); // Output cells
// Friend classes
};

extern CDecompHCells HCellsShow;

#endif
