//------------------------------------------------------------------------------------------------
// File: Decomp.cpp
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

#include "globals.h"
#include "ExchangeMPI.h"
#include "Decomp.h"
#include "CpuTree.h"
#include "tree.h"
#include "smatrix.h"
#include "SolverPar.h"

using namespace std;

// Author: Kharchenko S.A.
// Description: Compute distributions of the cells over CPU's for regular 3D mesh
// CDecomp::DistributeCellsOverCPUsForSortedTree()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsForSortedTree (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
									int _Nproc, double *_ProcWeight,
									int _Nx, int _Ny, int _Nz, int _Collap,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU) {

	const char *funcname = "DistributeCellsOverCPUsForSortedTree";

// Create collapsed XYZ partitioning

	int nxloc = (_Nx+_Collap-1) / _Collap;
	int nyloc = (_Ny+_Collap-1) / _Collap;
	int nzloc = (_Nz+_Collap-1) / _Collap;

	int *blksx;
	int *blksy;
	int *blksz;

	blksx = new int [nxloc+1];
	if (!blksx) MemoryFail (funcname);
	blksy = new int [nyloc+1];
	if (!blksy) MemoryFail (funcname);
	blksz = new int [nzloc+1];
	if (!blksz) MemoryFail (funcname);

	int i;

	for (i=0;i<nxloc;i++) blksx[i] = i*_Collap;
	blksx[nxloc] = _Nx;
	for (i=0;i<nyloc;i++) blksy[i] = i*_Collap;
	blksy[nyloc] = _Ny;
	for (i=0;i<nzloc;i++) blksz[i] = i*_Collap;
	blksz[nzloc] = _Nz;

// Compute sparsity of the collapsed matrix

	CSMatrix axyz = CSMatrix::Sparsity3D (nxloc, nyloc, nzloc);

// Compute cpu's partitioning

	int ntot = axyz.GetN ();
	int nztot = axyz.GetNzja ();
	int *ialoc = axyz.GetIa ();
	int *jaloc = axyz.GetJa ();

	int *CellsWeight;
	int *JaNeiWeight;

	CellsWeight = new int [ntot];
	if (!CellsWeight) MemoryFail (funcname);
	JaNeiWeight = new int [nztot];
	if (!JaNeiWeight) MemoryFail (funcname);

	int collap3 = _Collap * _Collap * _Collap;

	for (i=0;i<ntot;i++) CellsWeight[i] = collap3;
	for (i=0;i<nztot;i++) JaNeiWeight[i] = collap3;

	int info;

	int *CellsOrder;
	int *HyperCells;

	info = DistributeCellsOverCPUsForSortedTree (_comm,
										_Nproc, _ProcWeight, ntot, CellsWeight, 
										ialoc, jaloc, JaNeiWeight, _Param,
										_CpuTree, CellsOrder, 
										_NHyperCells, HyperCells, _HyperCell2CPU);
	if (info != 0) return info;

// Compute inverse order

	int nxyz = _Nx * _Ny * _Nz;

	int *iorder;

	iorder = new int [nxyz];
	if (!iorder) MemoryFail (funcname);

	int j;

	for (i=0;i<ntot;i++) iorder[CellsOrder[i]] = i;

// Create global arrays CellsOrder and HyperCells

	int nxy = _Nx * _Ny;

	_CellsOrder = new int [nxyz];
	if (!_CellsOrder) MemoryFail (funcname);
	_HyperCells = new int [_NHyperCells+1];
	if (!_HyperCells) MemoryFail (funcname);

	int nz = 0;
	_HyperCells[0] = 0;

	int iblk, ind, ixc, iyc, izc, itemp, ind1;
	int ix, iy, iz;

	for (iblk=0;iblk<_NHyperCells;iblk++) {
		for (j=HyperCells[iblk];j<HyperCells[iblk+1];j++) {
			ind = iorder[j];
			itemp = ind / nxloc;
			ixc = ind - itemp * nxloc;
			izc = itemp / nyloc;
			iyc = itemp - izc * nyloc;
			for (ix=blksx[ixc];ix<blksx[ixc+1];ix++) {
				for (iy=blksy[iyc];iy<blksy[iyc+1];iy++) {
					for (iz=blksz[izc];iz<blksz[izc+1];iz++) {
						ind1 = iz*nxy+iy*_Nx+ix;
						_CellsOrder[nz] = ind1;
						nz++;
					};
				};
			};
		};
		_HyperCells[iblk+1] = nz;
	};

	for (i=0;i<nxyz;i++) iorder[_CellsOrder[i]] = i;
	for (i=0;i<nxyz;i++) _CellsOrder[i] = iorder[i];

// Free work arrays

	delete [] blksx;
	delete [] blksy;
	delete [] blksz;
	delete [] CellsWeight;
	delete [] JaNeiWeight;
	delete [] CellsOrder;
	delete [] HyperCells;
	delete [] iorder;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's for regular 3D mesh
// CDecomp::DistributeCellsOverCPUs()
//========================================================================================
int CDecomp::DistributeCellsOverCPUs (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
									int _Nproc, double *_ProcWeight,
									int _Nx, int _Ny, int _Nz, int _Collap,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU) {

	const char *funcname = "DistributeCellsOverCPUs";

// Create collapsed XYZ partitioning

	int nxloc = (_Nx+_Collap-1) / _Collap;
	int nyloc = (_Ny+_Collap-1) / _Collap;
	int nzloc = (_Nz+_Collap-1) / _Collap;

	int *blksx;
	int *blksy;
	int *blksz;

	blksx = new int [nxloc+1];
	if (!blksx) MemoryFail (funcname);
	blksy = new int [nyloc+1];
	if (!blksy) MemoryFail (funcname);
	blksz = new int [nzloc+1];
	if (!blksz) MemoryFail (funcname);

	int i;

	for (i=0;i<nxloc;i++) blksx[i] = i*_Collap;
	blksx[nxloc] = _Nx;
	for (i=0;i<nyloc;i++) blksy[i] = i*_Collap;
	blksy[nyloc] = _Ny;
	for (i=0;i<nzloc;i++) blksz[i] = i*_Collap;
	blksz[nzloc] = _Nz;

// Compute sparsity of the collapsed matrix

	CSMatrix axyz = CSMatrix::Sparsity3D (nxloc, nyloc, nzloc);

// Compute cpu's partitioning

	int ntot = axyz.GetN ();
	int nztot = axyz.GetNzja ();
	int *ialoc = axyz.GetIa ();
	int *jaloc = axyz.GetJa ();

	int *CellsWeight;
	int *JaNeiWeight;

	CellsWeight = new int [ntot];
	if (!CellsWeight) MemoryFail (funcname);
	JaNeiWeight = new int [nztot];
	if (!JaNeiWeight) MemoryFail (funcname);

	int collap3 = _Collap * _Collap * _Collap;

	for (i=0;i<ntot;i++) CellsWeight[i] = collap3;
	for (i=0;i<nztot;i++) JaNeiWeight[i] = collap3;

	int info;

	int *CellsOrder;
	int *HyperCells;

	info = DistributeCellsOverCPUs (_comm,
										_Nproc, _ProcWeight, ntot, CellsWeight, 
										ialoc, jaloc, JaNeiWeight, _Param,
										_CpuTree, CellsOrder, 
										_NHyperCells, HyperCells, _HyperCell2CPU);
	if (info != 0) return info;

// Compute inverse order

	int nxyz = _Nx * _Ny * _Nz;

	int *iorder;

	iorder = new int [nxyz];
	if (!iorder) MemoryFail (funcname);

	int j;

	for (i=0;i<ntot;i++) iorder[CellsOrder[i]] = i;

// Create global arrays CellsOrder and HyperCells

	int nxy = _Nx * _Ny;

	_CellsOrder = new int [nxyz];
	if (!_CellsOrder) MemoryFail (funcname);
	_HyperCells = new int [_NHyperCells+1];
	if (!_HyperCells) MemoryFail (funcname);

	int nz = 0;
	_HyperCells[0] = 0;

	int iblk, ind, ixc, iyc, izc, itemp, ind1;
	int ix, iy, iz;

	for (iblk=0;iblk<_NHyperCells;iblk++) {
		for (j=HyperCells[iblk];j<HyperCells[iblk+1];j++) {
			ind = iorder[j];
			itemp = ind / nxloc;
			ixc = ind - itemp * nxloc;
			izc = itemp / nyloc;
			iyc = itemp - izc * nyloc;
			for (ix=blksx[ixc];ix<blksx[ixc+1];ix++) {
				for (iy=blksy[iyc];iy<blksy[iyc+1];iy++) {
					for (iz=blksz[izc];iz<blksz[izc+1];iz++) {
						ind1 = iz*nxy+iy*_Nx+ix;
						_CellsOrder[nz] = ind1;
						nz++;
					};
				};
			};
		};
		_HyperCells[iblk+1] = nz;
	};

	for (i=0;i<nxyz;i++) iorder[_CellsOrder[i]] = i;
	for (i=0;i<nxyz;i++) _CellsOrder[i] = iorder[i];

// Free work arrays

	delete [] blksx;
	delete [] blksy;
	delete [] blksz;
	delete [] CellsWeight;
	delete [] JaNeiWeight;
	delete [] CellsOrder;
	delete [] HyperCells;
	delete [] iorder;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's
// CDecomp::DistributeCellsOverCPUs()
//========================================================================================
int CDecomp::DistributeCellsOverCPUs (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU, int *&_HyperCellWeight) {

	const char *funcname = "DistributeCellsOverCPUs";

// Compute cpu's partitioning

	int info;

	info = DistributeCellsOverCPUs (_comm,
										_Nproc, _ProcWeight, _NCells, _CellsWeight, 
										_IANeiCells, _JANeiCells, _JANeiWeight, _Param,
										_CpuTree, _CellsOrder, 
										_NHyperCells, _HyperCells, _HyperCell2CPU);
	if (info != 0) return info;

// Compute additional statistics array

	_HyperCellWeight = new int [_NHyperCells];
	if (!_HyperCellWeight) MemoryFail (funcname);

	int iHCell;

	for (iHCell=0;iHCell<_NHyperCells;iHCell++) _HyperCellWeight[iHCell] = 0;

	int *iorder;

	iorder = new int [_NCells];
	if (!iorder) MemoryFail (funcname);

	int i, j;

	for (i=0;i<_NCells;i++) iorder[_CellsOrder[i]] = i;

	for (iHCell=0;iHCell<_NHyperCells;iHCell++) {
		for (i=_HyperCells[iHCell];i<_HyperCells[iHCell+1];i++) {
			j = iorder[i];
			_HyperCellWeight[iHCell] += _CellsWeight[j];
		};
	};

	delete [] iorder;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's
// CDecomp::DistributeCellsOverCPUsForSortedTree()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsForSortedTree (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU, int *&_HyperCellWeight) {

	const char *funcname = "DistributeCellsOverCPUsForSortedTree";

// Compute cpu's partitioning

	int info;

	info = DistributeCellsOverCPUsForSortedTree (_comm,
										_Nproc, _ProcWeight, _NCells, _CellsWeight, 
										_IANeiCells, _JANeiCells, _JANeiWeight, _Param,
										_CpuTree, _CellsOrder, 
										_NHyperCells, _HyperCells, _HyperCell2CPU);
	if (info != 0) return info;

// Compute additional statistics array

	_HyperCellWeight = new int [_NHyperCells];
	if (!_HyperCellWeight) MemoryFail (funcname);

	int iHCell;

	for (iHCell=0;iHCell<_NHyperCells;iHCell++) _HyperCellWeight[iHCell] = 0;

	int *iorder;

	iorder = new int [_NCells];
	if (!iorder) MemoryFail (funcname);

	int i, j;

	for (i=0;i<_NCells;i++) iorder[_CellsOrder[i]] = i;

	for (iHCell=0;iHCell<_NHyperCells;iHCell++) {
		for (i=_HyperCells[iHCell];i<_HyperCells[iHCell+1];i++) {
			j = iorder[i];
			_HyperCellWeight[iHCell] += _CellsWeight[j];
		};
	};

	delete [] iorder;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's
// CDecomp::DistributeCellsOverCPUs()
//========================================================================================
int CDecomp::DistributeCellsOverCPUs (CMPIComm &_comm, // Compute distributions of the cells over CPU's
										int _Nproc, double *_ProcWeight,
										int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
										CSolverPar &_Param,
										CCpuTree &_CpuTree, int *&_CellsOrder, 
										int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU) {

	const char *funcname = "DistributeCellsOverCPUs";

// External try/catch block

	try {

// Take the number of cpus and cpu id

		int nproc;
		int myid;

		nproc = _comm.GetNproc ();
		myid = _comm.GetMyid ();

//	if (nproc != _Nproc) throw " DistributeCellsOverCPUs: wrong number of processors on entry";

// Create matrix for computing decomposition data

		int nz = _IANeiCells[_NCells];

		CSMatrix aloc (_NCells,nz);

		int i;

		for (i=0;i< _NCells;i++) aloc.list[i] = i;
		for (i=0;i<=_NCells;i++) aloc.ia[i] = _IANeiCells[i];
		for (i=0;i<nz;i++) aloc.ja[i] = _JANeiCells[i];

// Sort column indices and weights

		CIntInt *iiarr;

		iiarr = new CIntInt [_NCells];
		if (!iiarr) MemoryFail (funcname);

		int j, jloc, ibeg;

		for (i=0;i<_NCells;i++) {
			ibeg = aloc.ia[i];
			for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
				jloc = j-ibeg;
				iiarr[jloc].intvalue = aloc.ja[j];
				iiarr[jloc].int2value = _JANeiWeight[j];
			};
	//		qsort (aloc.ja+aloc.ia[i], aloc.ia[i+1]-aloc.ia[i], sizeof(int), compint);
			qsort (iiarr, aloc.ia[i+1]-aloc.ia[i], sizeof(CIntInt), CIntInt::CompareIntInt);
			for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
				jloc = j-ibeg;
				aloc.ja[j] = iiarr[jloc].intvalue;
				_JANeiWeight[j] = iiarr[jloc].int2value;
			};
		};

		delete [] iiarr;

// Create local tree of processors

		int nchilds = _Param.slvparam.nchilds;

		double *memory2cpu;

		memory2cpu = new double [_Nproc+1];
		if (!memory2cpu) MemoryFail (funcname);

		for (i=0;i<_Nproc;i++) memory2cpu[i] = 1.0e0;

		CTree treeloc (_Nproc, nchilds, _ProcWeight, memory2cpu);

		treeloc.SetMyid (myid);

		delete [] memory2cpu;

// Allocate CellsOrder array

		_CellsOrder = new int [_NCells];
		if (!_CellsOrder) MemoryFail (funcname);

// Perform initial tree partitioning and ordering

		CTree treelocext;
		int ncycle = _Param.slvparam.ncycle;

//		for(i=0;i<nz;i++) _JANeiWeight[i] = 1;

		if (ncycle < 0) {
			aloc.PartOrdMtr (treeloc, _CellsWeight, _JANeiWeight, _CellsOrder);
			treelocext = treeloc;
		} else {
			aloc.PartOrdMtrSchur (treeloc, treelocext, _CellsWeight, _JANeiWeight, _CellsOrder, _Param.slvparam);
		};

//		treeloc.SetNproc (_Nproc);
//		treelocext.SetNproc (_Nproc);

// Decompose the data

		int weightcellvolume   = _Param.slvparam.weightcellvolume;
		int weightcellboundary = _Param.slvparam.weightcellboundary;

		if (ncycle < 0) {
			aloc.PartitionOrdMtr (treelocext, _CellsWeight, _JANeiWeight, _CellsOrder, 
									weightcellvolume, weightcellboundary, _NHyperCells, _HyperCells);
		} else {
			aloc.PartitionOrdMtrSchur (treelocext, _CellsWeight, _JANeiWeight, _CellsOrder, 
									weightcellvolume, weightcellboundary, _NHyperCells, _HyperCells);
		};

// Compute and store array HyperCell2CPU

		_HyperCell2CPU = new int [_NHyperCells];
		if (!_HyperCell2CPU) MemoryFail (funcname);

		for (i=0;i<_NHyperCells;i++) _HyperCell2CPU[i] = -1;

		int inode;

		for (inode=0;inode<treelocext.nnodes;inode++) 
		{
			int itype = treelocext.nodes[inode].nodetype;
			int isupbeg = treelocext.nodes[inode].indbeg;
			int isupend = treelocext.nodes[inode].indend;
			if (isupbeg <= isupend && itype != exchangenode) {
				int ihcellbeg = -1;
				int ihcellend = -1;
				for (i=0;i<_NHyperCells;i++) 
				{
					if (_HyperCells[i] == isupbeg) ihcellbeg = i;
					if (_HyperCells[i+1]-1 == isupend) ihcellend = i;
				};
				if (ihcellbeg == -1 || ihcellend == -1) throw " DistributeCellsOverCPUs: hypercell was not found";
				int icpu = treelocext.nodes[inode].nodecpu;
				for (i=ihcellbeg;i<=ihcellend;i++) _HyperCell2CPU[i] = icpu;
			};
		};

// Transform supernodes to blocks

		treelocext.TreeInBlocks (_NHyperCells,_HyperCells);

// Collapse extended tree

		treeloc = treelocext.CollapseTree ();

// Store current tree in the cputree structure

		_CpuTree.tree    = treeloc;
		_CpuTree.treeext = treelocext;

		_CpuTree.SetComm (_comm);

	}

	catch (char *strerr) {
		std::cout << " Error happenned: " << strerr << std::endl;
		return 1;
	}
	catch (...) {
		std::cout << " Unknown error happenned: " << std::endl;
		return 1;
	};
	_CpuTree.SetComm (_comm);

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's
// CDecomp::DistributeCellsOverCPUsForSortedTree()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsForSortedTree (CMPIComm &_comm, // Compute distributions of the cells over CPU's
										int _Nproc, double *_ProcWeight,
										int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells, int *_JANeiWeight,
										CSolverPar &_Param,
										CCpuTree &_CpuTree, int *&_CellsOrder, 
										int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU) {

	const char *funcname = "DistributeCellsOverCPUsForSortedTree";

// External try/catch block

	try {

// Take the number of cpus and cpu id

		int nproc;
		int myid;

		nproc = _comm.GetNproc ();
		myid = _comm.GetMyid ();

//	if (nproc != _Nproc) throw " DistributeCellsOverCPUs: wrong number of processors on entry";

// Create matrix for computing decomposition data

		int nz = _IANeiCells[_NCells];

		CSMatrix aloc (_NCells,nz);

		int i;

		for (i=0;i< _NCells;i++) aloc.list[i] = i;
		for (i=0;i<=_NCells;i++) aloc.ia[i] = _IANeiCells[i];
		for (i=0;i<nz;i++) aloc.ja[i] = _JANeiCells[i];

// Sort column indices and weights

		CIntInt *iiarr;

		iiarr = new CIntInt [_NCells];
		if (!iiarr) MemoryFail (funcname);

		int j, jloc, ibeg;

		for (i=0;i<_NCells;i++) {
			ibeg = aloc.ia[i];
			for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
				jloc = j-ibeg;
				iiarr[jloc].intvalue = aloc.ja[j];
				iiarr[jloc].int2value = _JANeiWeight[j];
			};
	//		qsort (aloc.ja+aloc.ia[i], aloc.ia[i+1]-aloc.ia[i], sizeof(int), compint);
			qsort (iiarr, aloc.ia[i+1]-aloc.ia[i], sizeof(CIntInt), CIntInt::CompareIntInt);
			for (j=aloc.ia[i];j<aloc.ia[i+1];j++) {
				jloc = j-ibeg;
				aloc.ja[j] = iiarr[jloc].intvalue;
				_JANeiWeight[j] = iiarr[jloc].int2value;
			};
		};

		delete [] iiarr;

// Create local tree of processors

		int nchilds = _Param.slvparam.nchilds;

		double *memory2cpu;

		memory2cpu = new double [_Nproc+1];
		if (!memory2cpu) MemoryFail (funcname);

		for (i=0;i<_Nproc;i++) memory2cpu[i] = 1.0e0;

		CTree treeloc (_Nproc, nchilds, _ProcWeight, memory2cpu);

		treeloc.SetMyid (myid);

		delete [] memory2cpu;

// Allocate CellsOrder array

		_CellsOrder = new int [_NCells];
		if (!_CellsOrder) MemoryFail (funcname);

// Perform initial tree partitioning and ordering

		CTree treelocext;
		int ncycle = _Param.slvparam.ncycle;

		if (ncycle < 0) {
			aloc.PartOrdMtr (treeloc, _CellsWeight, _JANeiWeight, _CellsOrder);
			treelocext = treeloc;
		} else {
			aloc.PartOrdMtrSchur (treeloc, treelocext, _CellsWeight, _JANeiWeight, _CellsOrder, _Param.slvparam);
		};

		treeloc.SetNproc (_Nproc);
		treelocext.SetNproc (_Nproc);

// Sort the tree

		treelocext.SortNodes ();

// Compute collapsed tree

		treeloc = treelocext.CollapseGeneralTree ();

// Decompose the data

		int weightcellvolume   = _Param.slvparam.weightcellvolume;
		int weightcellboundary = _Param.slvparam.weightcellboundary;

		if (ncycle < 0) {
			aloc.PartitionOrdMtrForSortedTree (treelocext, _CellsWeight, _JANeiWeight, _CellsOrder, 
									weightcellvolume, weightcellboundary, _NHyperCells, _HyperCells);
		} else {
			aloc.PartitionOrdMtrSchurForSortedTree (treelocext, _CellsWeight, _JANeiWeight, _CellsOrder, 
									weightcellvolume, weightcellboundary, _NHyperCells, _HyperCells);
		};

// Compute and store array HyperCell2CPU

		_HyperCell2CPU = new int [_NHyperCells];
		if (!_HyperCell2CPU) MemoryFail (funcname);

		for (i=0;i<_NHyperCells;i++) _HyperCell2CPU[i] = -1;

		int inode;

		for (inode=0;inode<treelocext.nnodes;inode++) 
		{
			int itype = treelocext.nodes[inode].nodetype;
			int ihcellbeg = treelocext.nodes[inode].indbeg;
			int ihcellend = treelocext.nodes[inode].indend;
			if (itype != exchangenode) {
				int icpu = treelocext.nodes[inode].nodecpu;
				for (i=ihcellbeg;i<=ihcellend;i++) _HyperCell2CPU[i] = icpu;
			};
		};

// Collapse extended tree

		treeloc = treelocext.CollapseGeneralTree ();

// Store current tree in the cputree structure

		_CpuTree.tree    = treeloc;
		_CpuTree.treeext = treelocext;

// Catch error

	}

	catch (char *strerr) {
		std::cout << " Error happenned: " << strerr << std::endl;
		return 1;
	}
	catch (...) {
		std::cout << " Unknown error happenned: " << std::endl;
		return 1;
	};

	_CpuTree.SetComm (_comm);

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's
// CDecomp::DistributeCellsOverCPUsNestedDessectionWeights()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsNestedDessectionWeights (CMPIComm &_comm, // Compute distributions of the cells over CPU's
									int _Nproc, double *_ProcWeight,
									int _NCells, int *_CellsWeight, int *_IANeiCells, int *_JANeiCells,
									CSolverPar &_Param,
									CCpuTree &_CpuTree, int *&_CellsOrder, 
									int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU, int *&_HyperCellWeight) {

	const char *funcname = "DistributeCellsOverCPUsNestedDessectionWeights";

// Compute the total cells weight

	int weight_total = 0;

	int i;

	for (i=0;i<_NCells;i++) weight_total += _CellsWeight[i];

// Compute the power of 2 for ND distribution parametrization

	CSlvParam *pparam = _Param.GetParam ();

	int weightcellvolumeloc = pparam->weightcellvolume;

	int npower2 = (weight_total+weightcellvolumeloc-1) / weightcellvolumeloc;

	if (npower2 <= 0) npower2 = 1;

	int ip = 1;

	while (ip < npower2) ip *= 2;

	npower2 = ip;

	pparam->ncpuext = npower2;

// Compute initial cpu's partitioning without weights usage

	int nhc, *hc, *hc2cpu;

	int info;

	info = DistributeCellsOverCPUsNestedDessection (_comm,
										_Nproc, _ProcWeight, _NCells,
										_IANeiCells, _JANeiCells, _Param,
										_CpuTree, _CellsOrder, 
										nhc, hc, hc2cpu);

	if (info != 0) return info;

// Compute inverse order for weights computations

	int *iorder;

	iorder = new int [_NCells];
	if (!iorder) MemoryFail (funcname);

	int j;

	for (i=0;i<_NCells;i++) iorder[_CellsOrder[i]] = i;

// Allocate working arrays that support additional hypercells cutting and cut

	int *hcnew, *hcwork, *hcwork1, *imaskhc, *hco2n;

	hcnew = new int [_NCells+1];
	if (!hcnew) MemoryFail (funcname);
	hcwork = new int [_NCells+1];
	if (!hcwork) MemoryFail (funcname);
	hcwork1 = new int [_NCells+1];
	if (!hcwork1) MemoryFail (funcname);
	imaskhc = new int [_NCells+1];
	if (!imaskhc) MemoryFail (funcname);
	hco2n = new int [nhc+1];
	if (!hco2n) MemoryFail (funcname);

	int nhcnew = 0;

	hcnew[0] = 0;
	hco2n[0] = 0;

	int ihc, ibeg, iend, nhcloc, nsplit, ibeg1, iend1, weightloc, jj, nhcloc1;
	int niloc, nparts, niloc1, iproc, j1;
	bool is_completed;

	for (ihc=0;ihc<nhc;ihc++) {
		ibeg = hc[ihc];
		iend = hc[ihc+1]-1;
		hcwork[0] = 0;
		hcwork[1] = iend-ibeg+1;
		nhcloc = 1;
		is_completed = false;
		while (!is_completed) {

// Find hcells still to be cut and register them

			nsplit = 0;

			for (i=0;i<nhcloc;i++) {
				ibeg1 = hcwork[i];
				iend1 = hcwork[i+1]-1;
				weightloc = 0;
				for (j=ibeg1;j<=iend1;j++) {
					j1 = ibeg+j;
					jj = iorder[j1];
					weightloc += _CellsWeight[jj];
				};
				if (weightloc > weightcellvolumeloc && iend1 > ibeg1) {
					imaskhc[i] = weightloc;
					nsplit++;
				} else {
					imaskhc[i] = -1;
				};
			};
			if (nsplit == 0) is_completed = true;
			if (!is_completed) {
				nhcloc1 = 0;
				hcwork1[0] = 0;
				for (i=0;i<nhcloc;i++) {
					ibeg1 = hcwork[i];
					iend1 = hcwork[i+1]-1;
					niloc = iend1-ibeg1+1;
					if (imaskhc[i] > 0) {
						nparts = (imaskhc[i]+weightcellvolumeloc-1) / weightcellvolumeloc;
						if (nparts < 2) nparts = 2;
						niloc1 = niloc / nparts;
						if (niloc1 < 1) niloc = 1;
						for (j=0;j<nparts-1;j++) {
							hcwork1[nhcloc1+1] = hcwork1[nhcloc1] + niloc1;
							nhcloc1++;
						};
						hcwork1[nhcloc1+1] = hcwork1[nhcloc1] + (niloc-(nparts-1)*niloc1);
						nhcloc1++;
					} else {
						hcwork1[nhcloc1+1] = hcwork1[nhcloc1] + niloc;
						nhcloc1++;
					};
				};
				nhcloc = nhcloc1;
				for (i=0;i<=nhcloc;i++) hcwork[i] = hcwork1[i];
			};
		};
		for (i=0;i<nhcloc;i++) {
			ibeg1 = hcwork[i];
			iend1 = hcwork[i+1]-1;
			hcnew[nhcnew+1] = hcnew[nhcnew] + (iend1-ibeg1+1);
			nhcnew++;
		};
		hco2n[ihc+1] = nhcnew;
	};

// Recompute hypercells

	_NHyperCells = nhcnew;

	_HyperCells = new int [_NHyperCells+1];
	if (!_HyperCells) MemoryFail (funcname);
	_HyperCell2CPU = new int [_NHyperCells];
	if (!_HyperCell2CPU) MemoryFail (funcname);

	for (i=0;i<=_NHyperCells;i++) _HyperCells[i] = hcnew[i];

	for (ihc=0;ihc<nhc;ihc++) {
		iproc = hc2cpu[ihc];
		ibeg = hco2n[ihc];
		iend = hco2n[ihc+1]-1;
		for (i=ibeg;i<=iend;i++) _HyperCell2CPU[i] = iproc;
	};

// Register new hypercells

	CTree *ptree = _CpuTree.GetTree ();
	CTree *ptreeext = _CpuTree.GetTreeext ();

	ptree->ReassignBlockNumbers (hco2n);
	ptreeext->ReassignBlockNumbers (hco2n);

// Free work arrays

	delete [] hc;
	delete [] hc2cpu;
	delete [] hcnew;
	delete [] hcwork;
	delete [] hcwork1;
	delete [] imaskhc;
	delete [] hco2n;

// Compute additional statistics array

	_HyperCellWeight = new int [_NHyperCells];
	if (!_HyperCellWeight) MemoryFail (funcname);

	int iHCell;

	for (iHCell=0;iHCell<_NHyperCells;iHCell++) _HyperCellWeight[iHCell] = 0;

	for (iHCell=0;iHCell<_NHyperCells;iHCell++) {
		for (i=_HyperCells[iHCell];i<_HyperCells[iHCell+1];i++) {
			j = iorder[i];
			_HyperCellWeight[iHCell] += _CellsWeight[j];
		};
	};

	delete [] iorder;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's
// CDecomp::DistributeCellsOverCPUsNestedDessection()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsNestedDessection (CMPIComm &_comm, // Compute distributions of the cells over CPU's
										int _Nproc, double *_ProcWeight,
										int _NCells, int *_IANeiCells, int *_JANeiCells,
										CSolverPar &_Param,
										CCpuTree &_CpuTree, int *&_CellsOrder, 
										int &_NHyperCells, int *&_HyperCells, int *&_HyperCell2CPU) {

	const char *funcname = "DistributeCellsOverCPUsNestedDessection";

// Take the number of cpus and cpu id

	int nproc;
	int myid;

	nproc = _comm.GetNproc ();
	myid = _comm.GetMyid ();

// Modify the number of processors used for distribution

	int nproc_power2 = 1;
	while (nproc_power2*2 <= _Nproc) nproc_power2 *= 2;

// Create matrix for computing decomposition data

	int nz = _IANeiCells[_NCells];

	CSMatrix aloc (_NCells,nz);

	int i;

	for (i=0;i< _NCells;i++) aloc.list[i] = i;
	for (i=0;i<=_NCells;i++) aloc.ia[i] = _IANeiCells[i];
	for (i=0;i<nz;i++) aloc.ja[i] = _JANeiCells[i];

// Sort column indices and weights

	for (i=0;i<_NCells;i++) {
		qsort (aloc.ja+aloc.ia[i], aloc.ia[i+1]-aloc.ia[i], sizeof(int), compint);
	};

// Create local and extended tree of processors

	int nchilds = 2;

	int nprocext = nproc_power2;

	if (nprocext < _Param.GetNcpuext ()) nprocext = _Param.GetNcpuext ();

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nprocext+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocext+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nprocext;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocext;i++) procweight[i] = 1.0e0;

	CTree treeloc (nproc_power2, nchilds, _ProcWeight, memory2cpu);
	CTree treelocexthc (nprocext, nchilds, procweight, memory2cpu);

	treeloc.SetMyid (myid);
	treelocexthc.SetMyid (myid);

	delete [] memory2cpu;
	delete [] procweight;

// Sort the tree

	treeloc.SortNodes ();
	treelocexthc.SortNodes ();

// Allocate CellsOrder array

	_CellsOrder = new int [_NCells];
	if (!_CellsOrder) MemoryFail (funcname);

// Compute ND ordering and reorder

	CSlvParam *pparam = _Param.GetParam ();

	bool schur = false;

	if (pparam->schurtype == 1) {
		schur = true;
	};

	CSMatrix asord;

	if (!schur) {

		aloc.OrderPrfMtr (-1, _CellsOrder);
		asord = aloc.OrdMtr (_CellsOrder);
		asord.PartBinaryTreeND (treelocexthc);

	} else {

		CSlvParam *pparam = _Param.GetParam ();

		aloc.PartBinaryTreeNDSchur (treelocexthc, _CellsOrder, *pparam);
		asord = aloc.OrdMtr (_CellsOrder);

	};

// Return HC partitioning and CPU distribution

	treelocexthc.PartitionTreeSchur (_NHyperCells, _HyperCells, _HyperCell2CPU);

// Compute current computations tree

	treelocexthc.IncludeBinaryTreeSchur (treeloc, _NHyperCells, _HyperCells, asord);

// Add blocks partitioning into the tree

	int *pblkstree = treeloc.GetBlkstree ();

	delete [] pblkstree;

	pblkstree = new int [_NHyperCells+1];
	if (!pblkstree) MemoryFail (funcname);

	for (i=0;i<=_NHyperCells;i++) pblkstree[i] = _HyperCells[i];

	treeloc.SetNblkstree (_NHyperCells);
	treeloc.SetBlkstree (pblkstree);

// Init additional trees

	CTree treelocext;

	treelocext = treeloc;

	treeloc.SetNproc (nproc_power2);
	treelocext.SetNproc (nproc_power2);
	treelocexthc.SetNproc (nprocext);

// Add list of cpus for all nodes

	treeloc.CpuListsForNodesSchur ();
	treelocext.CpuListsForNodesSchur ();

// Compute current distribution

	treeloc.Block2CpuSchur (_HyperCell2CPU);

//	ofstream fout ("ChkDecomp.dat");

//	OutArr (fout, " HyperCells ",_NHyperCells+1,_HyperCells);
//	OutArr (fout, " HyperCells2CPU ",_NHyperCells,_HyperCell2CPU);
//	OutArr (fout, " Order ",_NCells,_CellsOrder);

// Store current tree in the cputree structure

	_CpuTree.tree    = treeloc;
	_CpuTree.treeext = treelocext;
	_CpuTree.treehc  = treelocexthc;

//	_CpuTree.SetCommSchur (_comm);

	return 0;

};

// Author: Kharchenko S.A.
// Description: Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
// CDecomp::DistributeCellsOverCPUsNestedDessectionParMetis()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsNestedDessectionParMetis (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																					int _N1LevelCells, int *_Index1LevelCells,
																					int *_IA1Level2Real, int *_Real2Index,
																					int *_IARealCells, int *_JAReal2Index,
																					CSolverPar &_Param, 
																					CCpuTree &_CpuTree,
																					int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																					int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
																					int *&_Index1LevelCellsNew, int *&_IA1Level2RealNew) {

	const char *funcname = "DistributeCellsOverCPUsNestedDessectionParMetis";

	int myid = _comm.GetMyid ();
	int nproc = _comm.GetNproc ();

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkDecomp_",myid,".dat");
//	ofstream fout (strbuff);

// Sort the indices of 1 level cells

	CIntInt *iiarr;

	iiarr = new CIntInt [_N1LevelCells];
	if (!iiarr) MemoryFail (funcname);

	int i;

	for (i=0;i<_N1LevelCells;i++) {
		iiarr[i].intvalue = _Index1LevelCells[i*2];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, _N1LevelCells, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *listindex, *attrindex, *ordindexinv;

	listindex = new int [_N1LevelCells];
	if (!listindex) MemoryFail (funcname);
	attrindex = new int [_N1LevelCells];
	if (!attrindex) MemoryFail (funcname);
	ordindexinv = new int [_N1LevelCells];
	if (!ordindexinv) MemoryFail (funcname);

	for (i=0;i<_N1LevelCells;i++) {
		listindex[i] = iiarr[i].intvalue;
		ordindexinv[i] = iiarr[i].int2value;
	};

	delete [] iiarr;

// Check that array Index1LevelCells is correct

	int indexmax = 0;
	if (_N1LevelCells > 0) indexmax = listindex[_N1LevelCells-1];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, MAXIMUM, 1, &indexmax, &indexmax);

	int ncellstot = _N1LevelCells;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &ncellstot, &ncellstot);

	if (indexmax != ncellstot-1) throw " CDecomp::DistributeCellsOverCPUsNestedDessectionParMetis: array Index1LevelCells on entry is incorrect ";

// Create ordered sprnds

	int *sprnds;
	int *sprndsord;

	sprnds = new int [_N1LevelCells+1];
	if (!sprnds) MemoryFail (funcname);
	sprndsord = new int [_N1LevelCells+1];
	if (!sprndsord) MemoryFail (funcname);

	sprnds[0] = 0;
	sprndsord[0] = 0;

	int ind;

	for (i=0;i<_N1LevelCells;i++) {
		ind = ordindexinv[i];
		sprnds[i+1] = sprnds[i] + (_IA1Level2Real[i+1]-_IA1Level2Real[i]);
		sprndsord[i+1] = sprndsord[i] + (_IA1Level2Real[ind+1]-_IA1Level2Real[ind]);
	};

// Partition the indices into the parts for balanced load

	CSlvParam *pparam = _Param.GetParam ();

	int nprocpart = pparam->ncpupart;

	if (nprocpart <= 0) {
		nprocpart = 1;
		while (nprocpart*2 <= nproc) nprocpart *= 2;
	};

	if (nprocpart > nproc) {
		nprocpart = 1;
		while (nprocpart*2 <= nproc) nprocpart *= 2;
	};

	int *blksindex;

	blksindex = new int [nprocpart+1];
	if (!blksindex) MemoryFail (funcname);

	PartitionIndices (_comm, 
							nprocpart, _N1LevelCells, listindex, sprndsord,
							blksindex);

// Create CSMatrix 2index structure that is based on pair {1levelcell,realcell} instead of {hypercell,realcell1}

	CSMatrix mtr;

	for (i=0;i<_N1LevelCells;i++) {
		attrindex[i] = _Index1LevelCells[2*i];
	};

	mtr = TransformSparsity (_comm, 
										_N1LevelCells, sprnds, attrindex,
										_Real2Index,
										_IARealCells, _JAReal2Index);

// Exchange matrix data

	int nint = 0;
	int ibeg = 0;
	int iend = -1;

	if (myid < nprocpart) {
		nint = 1;
		ibeg = blksindex[myid];
		iend = blksindex[myid+1]-1;
	};

	int nsuploc = iend-ibeg+1;

	int *listsp;
	int *ordersp;
	int *attrinterval;

	listsp = new int [nsuploc];
	if (!listsp) MemoryFail (funcname);
	ordersp = new int [nsuploc];
	if (!ordersp) MemoryFail (funcname);
	attrinterval = new int [nsuploc];
	if (!attrinterval) MemoryFail (funcname);

	for (i=0;i<nsuploc;i++) listsp[i] = ibeg+i;

	CSMatrix mtrnew;

	mtrnew = mtr.GetSubmatrixViaIntervals2Index (_comm,
																nint, &ibeg, &iend);

	CSMatrix mtrdummy;

	mtr = mtrdummy;

// Compute global attributes array for the condenced interval [ibeg,iend]

	int *index2cpuini;

	index2cpuini = new int [_N1LevelCells];
	if (!index2cpuini) MemoryFail (funcname);

	int *iacpu;
	int *iptrcpu;
	int **plists;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);
	plists = new int * [nproc];
	if (!plists) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iacpu[i] = 0;

	int iproc, ni, icpubeg, icpuend, icpu;

	for (i=0;i<_N1LevelCells;i++) {
		ind = _Index1LevelCells[i*2];
		icpubeg = 0;
		icpuend = nprocpart-1;
		while (icpubeg != icpuend) {
			icpu = (icpubeg+icpuend)/2;
			if (ind >= blksindex[icpu] && ind < blksindex[icpu+1]) {
				icpubeg = icpu;
				icpuend = icpu;
			} else if (ind < blksindex[icpu]) {
				icpuend = icpu-1;
			} else if (ind >= blksindex[icpu+1]) {
				icpubeg = icpu+1;
			};
		};
		index2cpuini[i] = icpubeg;
		iacpu[icpubeg+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	for (i=0;i<nproc;i++) {
		ni = iacpu[i+1]-iacpu[i];
		plists[i] = new int [ni*2];
		if (!plists[i]) MemoryFail (funcname);
	};

	int attr, k;

	for (i=0;i<_N1LevelCells;i++) {
		ind = _Index1LevelCells[i*2];
		attr = _Index1LevelCells[i*2+1];
		iproc = index2cpuini[i];
		k = iptrcpu[iproc]-iacpu[iproc];
		plists[iproc][k*2] = ind;
		plists[iproc][k*2+1] = attr;
		iptrcpu[iproc]++;
	};

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nproc;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char* [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 2*(iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) plists[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CDecomp::DistributeCellsOverCPUsNestedDessectionParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	int *piarr;

	int j;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*2];
			attr = piarr[j*2+1];
			attrinterval[ind-ibeg] = attr;
		};
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

	for (i=0;i<nproc;i++) delete [] plists[i];

// Create supernodes partitioning

	int *plist2 = mtrnew.GetList2 ();
	int nlistnew = mtrnew.GetNlist ();

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	nsuploc = 0;
	int isupprev = -1;

	sprndsloc[0] = 0;

	if (nlistnew > 0) {
		isupprev = plist2[0];
		listsp[nsuploc] = isupprev;
		nsuploc = 1;
		sprndsloc[nsuploc] = 1;
	};

	int isup;

	for (i=1;i<nlistnew;i++) {
		sprndsloc[nsuploc] = i+1;
		isup = plist2[i];
		if (isup != isupprev) {
			sprndsloc[nsuploc] = i;
			isupprev = isup;
			listsp[nsuploc] = isupprev;
			nsuploc++;
		};
	};

	sprndsloc[nsuploc] = nlistnew;

// Renumber the matrix data (condense pairs numbering into the index numbering)

	mtrnew.Transform2IndexIntoContinuous (_comm);

// Create local communicator

	int *listcpu;

	listcpu = new int [nprocpart];
	if (!listcpu) MemoryFail (funcname);

	for (i=0;i<nprocpart;i++) listcpu[i] = i;

	CMPIComm commloc;

	commloc = _comm.CreateComm (nprocpart, listcpu);

	delete [] listcpu;

// Distribute via ParMetis

	CTree treeloc;

	if (myid < nprocpart) {

// Distribute point data

//		int nlistloc = mtrnew.GetNlist ();

		int nprocdist = 1;

		while (nprocdist*2 <= nproc) nprocdist *= 2;

		int *orderloc;
		int *blksloc;

		mtrnew.OrderNDParMetis (commloc, nprocdist,
										treeloc, orderloc, _NHyperCells, blksloc, _HyperCell2CPU);

		mtrnew = mtrdummy;

// Compute supernodes ordering and corresponding blocks partitioning in terms of supernodes

		_HyperCells = new int [_NHyperCells+1];
		if (!_HyperCells) MemoryFail (funcname);

		SpPartitioning (commloc, 
								nsuploc, listsp, sprndsloc,
								orderloc, _NHyperCells, blksloc, _HyperCell2CPU,
								ordersp, _HyperCells);

		delete [] orderloc;
		delete [] blksloc;

		if (myid > 0) delete [] _HyperCells;
		if (myid > 0) delete [] _HyperCell2CPU;

	};

// Exchange tree and blocks information

	int isize = 0;
	char *ptree;

	if (myid == 0) {
		treeloc.PackTree (isize,ptree);
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &isize, &isize);

	if (myid > 0) {
		ptree = new char [isize];
		if (!ptree) MemoryFail (funcname);
		for (i=0;i<isize;i++) ptree[i] = 0;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, CHARVALUE, ADD, isize, ptree, ptree);

	treeloc.UnPackTree (isize,ptree);

	delete [] ptree;

	if (myid > 0) {
		_NHyperCells = 0;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &_NHyperCells, &_NHyperCells);

	if (myid > 0) {
		_HyperCells = new int [_NHyperCells+1];
		if (!_HyperCells) MemoryFail (funcname);
		_HyperCell2CPU = new int [_NHyperCells];
		if (!_HyperCell2CPU) MemoryFail (funcname);

		for (i=0;i<_NHyperCells+1;i++) _HyperCells[i] = 0;
		for (i=0;i<_NHyperCells;i++) _HyperCell2CPU[i] = 0;

	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells+1, _HyperCells, _HyperCells);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells, _HyperCell2CPU, _HyperCell2CPU);

// For each supernode find its block number

	int *spl2blk;

	spl2blk = new int [nsuploc];
	if (!spl2blk) MemoryFail (funcname);

	int isupnew, iblkbeg, iblkend, iblk;

	for (i=0;i<nsuploc;i++) {
		isup = ibeg+i;
		isupnew = ordersp[i];
		iblkbeg = 0;
		iblkend = _NHyperCells-1;
		while (iblkbeg != iblkend) {
			iblk = (iblkbeg+iblkend)/2;
			if (isupnew >= _HyperCells[iblk] && isupnew < _HyperCells[iblk+1]) {
				iblkbeg = iblk;
				iblkend = iblk;
			} else {
				if (isupnew < _HyperCells[iblk]) iblkend = iblk-1;
				if (isupnew >= _HyperCells[iblk+1]) iblkbeg = iblk+1;
			};
		};
		spl2blk[i] = iblkbeg;
	};

// Create the lists for exchange

	for (i=0;i<nproc+1;i++) iacpu[i] = 0;

	for (i=0;i<nsuploc;i++) {
		iblk = spl2blk[i];
		iproc = _HyperCell2CPU[iblk];
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	for (i=0;i<nproc;i++) {
		ni = iacpu[i+1]-iacpu[i];
		plists[i] = new int [ni*5];
		if (!plists[i]) MemoryFail (funcname);
	};

	for (i=0;i<nsuploc;i++) {
		iblk = spl2blk[i];
		isup = ibeg+i;
		attr = attrinterval[i];
		isupnew = ordersp[i];
		isize = sprndsloc[i+1]-sprndsloc[i];
		iproc = _HyperCell2CPU[iblk];
		k = iptrcpu[iproc];
		ind = k-iacpu[iproc];
		plists[iproc][ind*5] = iblk;
		plists[iproc][ind*5+1] = isup;
		plists[iproc][ind*5+2] = attr;
		plists[iproc][ind*5+3] = isupnew;
		plists[iproc][ind*5+4] = isize;
		iptrcpu[iproc]++;
	};

// Exchange the lists

	NObjSend = nproc;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char* [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 5*(iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) plists[i];
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CDecomp::DistributeCellsOverCPUsNestedDessectionParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	_N1LevelCellsNew = 0;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (5*sizeof(int));
		_N1LevelCellsNew += ni;
	};

	_1LevelCellsNew2HCells = new int [_N1LevelCellsNew];
	if(!_1LevelCellsNew2HCells) MemoryFail(funcname);
	_Index1LevelCellsNew = new int [_N1LevelCellsNew*2];
	if(!_Index1LevelCellsNew) MemoryFail(funcname);
	_IA1Level2RealNew = new int [_N1LevelCellsNew+1];
	if(!_IA1Level2RealNew) MemoryFail(funcname);

	iiarr = new CIntInt [_N1LevelCellsNew];
	if(!iiarr) MemoryFail(funcname);

	_N1LevelCellsNew = 0;

	_IA1Level2RealNew[0] = 0;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (5*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			iblk = piarr[j*5];
			isup = piarr[j*5+1];
			attr = piarr[j*5+2];
			isupnew = piarr[j*5+3];
			isize = piarr[j*5+4];
			_1LevelCellsNew2HCells[_N1LevelCellsNew] = iblk;
			_Index1LevelCellsNew[_N1LevelCellsNew*2] = isup;
			_Index1LevelCellsNew[_N1LevelCellsNew*2+1] = attr;
			_IA1Level2RealNew[_N1LevelCellsNew+1] = isize;
			iiarr[_N1LevelCellsNew].intvalue = isupnew;
			iiarr[_N1LevelCellsNew].int2value = _N1LevelCellsNew;
			_N1LevelCellsNew++;
		};
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

	qsort (iiarr, _N1LevelCellsNew, sizeof(CIntInt), CIntInt::CompareIntInt);

	int *imask;

	imask = new int [_N1LevelCellsNew];
	if(!imask) MemoryFail(funcname);

	for (i=0;i<_N1LevelCellsNew;i++) {
		j = iiarr[i].int2value;
		imask[i] = _1LevelCellsNew2HCells[j];
	};
	for (i=0;i<_N1LevelCellsNew;i++) _1LevelCellsNew2HCells[i] = imask[i];

	for (i=0;i<_N1LevelCellsNew;i++) {
		j = iiarr[i].int2value;
		imask[i] = _Index1LevelCellsNew[2*j];
	};
	for (i=0;i<_N1LevelCellsNew;i++) _Index1LevelCellsNew[2*i] = imask[i];

	for (i=0;i<_N1LevelCellsNew;i++) {
		j = iiarr[i].int2value;
		imask[i] = _Index1LevelCellsNew[2*j+1];
	};
	for (i=0;i<_N1LevelCellsNew;i++) _Index1LevelCellsNew[2*i+1] = imask[i];

	for (i=0;i<_N1LevelCellsNew;i++) {
		j = iiarr[i].int2value;
		imask[i] = _IA1Level2RealNew[j+1];
	};
	for (i=0;i<_N1LevelCellsNew;i++) _IA1Level2RealNew[i+1] = imask[i];

	for (i=0;i<_N1LevelCellsNew;i++) _IA1Level2RealNew[i+1] = _IA1Level2RealNew[i]+_IA1Level2RealNew[i+1];

// Free work arrays

	delete [] iiarr;
	delete [] listindex;
	delete [] attrindex;
	delete [] ordindexinv;
	delete [] sprnds;
	delete [] sprndsord;
	delete [] blksindex;
	delete [] listsp;
	delete [] ordersp;
	delete [] attrinterval;
	delete [] index2cpuini;
	delete [] sprndsloc;
	delete [] spl2blk;
	delete [] iacpu;
	delete [] iptrcpu;
	for (i=0;i<nproc;i++) delete [] plists[i];
	delete [] plists;
	delete [] imask;

// Init cputree

	_CpuTree.SetTree (treeloc);
	_CpuTree.SetTreeExt (treeloc);
	_CpuTree.SetTreeHc (treeloc);

	return 0;

};

// Author: Kharchenko S.A.
// Description: Compute partitioning of the indices for balanced load
// CDecomp::PartitionIndices()
//========================================================================================
void CDecomp::PartitionIndices (CMPIComm &_comm, // Compute partitioning of the indices for balanced load
											int _nprocpart, int _ncells, int *_listindex, int *_sprnds,
											int *_blksindex) {

	const char *funcname = "PartitionIndices";

//	int myid = _comm.GetMyid ();
//	int nproc = _comm.GetNproc ();

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","ChkPart_",myid,".dat");
//	ofstream fout (strbuff);

// Compute the number of rows

	int ncellstot = _ncells;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &ncellstot, &ncellstot);

	int indcellmax = 0;

	if (_ncells > 0) indcellmax = _listindex[_ncells-1];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, MAXIMUM, 1, &indcellmax, &indcellmax);

	indcellmax += 1;

	int nloc = _sprnds[_ncells];
	int ntot = nloc;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &ntot, &ntot);

// Compute approximate number of elements in each part

	int ni = ntot / _nprocpart;

// Main cycle over the parts

	int npartsmax = 200;

	int *indbegarr, *indendarr;
	int *blksloc;

	indbegarr = new int [npartsmax];
	if (!indbegarr) MemoryFail (funcname);
	indendarr = new int [npartsmax];
	if (!indendarr) MemoryFail (funcname);
	blksloc = new int [npartsmax+1];
	if (!blksloc) MemoryFail (funcname);

	int *blkslist;

	blkslist = new int [_nprocpart+1];
	if (!blkslist) MemoryFail (funcname);

	_blksindex[0] = 0;
	blkslist[0] = 0;

	int iproc, indexbeg, indexend, icellbeg, icellend, icellbeg0;
	int nsuploc, nparts, nsupi, i, j, ind, ivalue, ibeg, iend, niprev;

	for (iproc=0;iproc<_nprocpart-1;iproc++) {
		indexbeg = _blksindex[iproc];
		indexend = indcellmax-1;
		icellend = _ncells-1;
		ibeg = -1;
		for (j=0;j<=icellend;j++) {
			ind = _listindex[j];
			if (ibeg == -1 && ind >= indexbeg) {
				ibeg = j;
				break;
			};
		};
		if (ibeg >= 0) {
			icellbeg = ibeg;
		} else {
			icellbeg = _ncells;
		};
		icellbeg0 = icellbeg;
		while (indexbeg != indexend) {

// Compute intervals partitioning

			nsuploc = indexend-indexbeg+1;
			if (nsuploc < npartsmax) {
				nparts = nsuploc;
			} else {
				nparts = npartsmax;
			};
			nsupi = (nsuploc+nparts-1) / nparts;
			for (i=0;i<nparts;i++) {
				indbegarr[i] = -1;
				indendarr[i] = -1;
			};

// Fill local intervals data

			for (j=icellbeg;j<=icellend;j++) {
				ind = _listindex[j];
				if (ind >= indexbeg && ind <= indexend) {
					ivalue = (ind - indexbeg) / nsupi;
					if (indbegarr[ivalue] == -1) indbegarr[ivalue] = j;
					indendarr[ivalue] = j;
				};
			};

// Count the local numbers of rows in each interval

			for (i=0;i<=nparts;i++) blksloc[i] = 0;

			for (i=0;i<nparts;i++) {
				if (indbegarr[i] >= 0) {
					ibeg = indbegarr[i];
					iend = indendarr[i];
					blksloc[i+1] = _sprnds[iend+1]-_sprnds[ibeg];
				};
			};

// Exchange interval data

			CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nparts+1, blksloc, blksloc);

			for (i=0;i<nparts;i++) blksloc[i+1] = blksloc[i]+blksloc[i+1];

// Count previous volume

			ibeg = icellbeg0;
			iend = icellbeg-1;
			niprev = _sprnds[iend+1]-_sprnds[ibeg];

			CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &niprev, &niprev);

// Modify control values

			ivalue = -1;

			for (i=0;i<nparts;i++) {
				if (niprev+blksloc[i+1] < ni) ivalue = i;
			};

			if (ivalue < 0) {
				iend = indexbeg + nsupi-1;
				indexend = iend;
				if (indexend > indcellmax-1) indexend = indcellmax-1;
			} else if (ivalue == nparts-1) {
				ibeg = indexbeg + nsupi*ivalue;
				indexbeg = ibeg;
				if (indexbeg > indexend) indexbeg = indexend;
			} else {
				ibeg = indexbeg + nsupi*ivalue;
				iend = indexbeg + nsupi*(ivalue+1)-1;
				if (iend > indcellmax-1) iend = indcellmax-1;
				indexbeg = ibeg;
				indexend = iend;
				if (indexbeg > indexend) indexbeg = indexend;
			};

			ibeg = -1;
			for (j=0;j<=icellend;j++) {
				ind = _listindex[j];
				if (ibeg == -1 && ind >= indexbeg) {
					ibeg = j;
					break;
				};
			};
			if (ibeg >= 0) {
				icellbeg = ibeg;
			} else {
				icellbeg = _ncells;
			};
		};
		_blksindex[iproc+1] = indexbeg+1;
		blkslist[iproc+1] = icellbeg+1;
	};

	_blksindex[_nprocpart] = indcellmax;

// Check array blksindex

	for (i=0;i<_nprocpart;i++) {
		if (_blksindex[i] > _blksindex[i+1]) throw " CDecomp::PartitionIndices: wrong partitioning ";
	};

// Free work arrays

	delete [] indbegarr;
	delete [] indendarr;
	delete [] blksloc;
	delete [] blkslist;

};

// Author: Kharchenko S.A.
// Description: Create CSMatrix 2index structure that is based on pair {1levelcell,realcell} instead of {hypercell,realcell1}
// CDecomp::TransformSparsity()
//========================================================================================
CSMatrix CDecomp::TransformSparsity (CMPIComm &_comm, // Create CSMatrix 2index structure that is based on pair {1levelcell,realcell} instead of {hypercell,realcell1}
												int _N1LevelCells, int *_sprnds, int *_Index1LevelCells,
												int *_Real2Index,
												int *_IARealCells, int *_JAReal2Index) {

	const char *funcname = "TransformSparsity";

	int myid = _comm.GetMyid ();
	int nproc = _comm.GetNproc ();

// Count the number of hypercells

	int nblkstot = 0;

	int i, j, ihcell, icell;

	for (i=0;i<_N1LevelCells;i++) {
		for (j=_sprnds[i];j<_sprnds[i+1];j++) {
			ihcell = _Real2Index[j*2];
			if (ihcell > nblkstot) nblkstot = ihcell;
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, MAXIMUM, 1, &nblkstot, &nblkstot);

	nblkstot += 1;

// Create distribution of the hypercells to cpus and sizes of hypercells

	int *blks2cpu;
	int *blksini;

	blks2cpu = new int [nblkstot+1];
	if (!blks2cpu) MemoryFail (funcname);
	blksini = new int [nblkstot+1];
	if (!blksini) MemoryFail (funcname);

	for (i=0;i<=nblkstot;i++) blks2cpu[i] = 0;
	for (i=0;i<=nblkstot;i++) blksini[i] = 0;

	for (i=0;i<_N1LevelCells;i++) {
		for (j=_sprnds[i];j<_sprnds[i+1];j++) {
			ihcell = _Real2Index[j*2];
			icell = _Real2Index[j*2+1];
			blks2cpu[ihcell] = myid;
			if (icell+1 > blksini[ihcell+1]) blksini[ihcell+1] = icell+1;
		};
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nblkstot+1, blks2cpu, blks2cpu);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nblkstot+1, blksini, blksini);

	for (i=0;i<nblkstot;i++) blksini[i+1] = blksini[i] + blksini[i+1];

// Split sparsity into processors description

	int ncellsloc = _sprnds[_N1LevelCells];
	int nzjaloc = _IARealCells[ncellsloc];

	int *iacpu;
	int *iptrcpu;
	int *jaindcpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptrcpu = new int [nproc+1];
	if (!iptrcpu) MemoryFail (funcname);
	jaindcpu = new int [nzjaloc];
	if (!jaindcpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	int iproc, k, ibs;

	for (i=0;i<nzjaloc;i++) {
		ihcell = _JAReal2Index[i*2];
		iproc = blks2cpu[ihcell];
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i] + iacpu[i+1];
	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	for (i=0;i<nzjaloc;i++) {
		ihcell = _JAReal2Index[i*2];
		iproc = blks2cpu[ihcell];
		k = iptrcpu[iproc];
		jaindcpu[k] = i;
		iptrcpu[iproc]++;
	};

// Determine the maximal size of the mask array

	int *ibsblk;

	ibsblk = new int [nblkstot];
	if (!ibsblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) iptrcpu[i] = 0;

	for (i=0;i<nblkstot;i++) {
		iproc = blks2cpu[i];
		ibsblk[i] = iptrcpu[iproc];
		iptrcpu[iproc] += (blksini[i+1]-blksini[i]);
	};

	int nimax = 0;

	for (i=0;i<nproc;i++) {
		if (iptrcpu[i] > nimax) nimax = iptrcpu[i];
	};

// Register existing pairs

	int *i2index;

	i2index = new int [nimax*2];
	if (!i2index) MemoryFail (funcname);

	for (i=0;i<nimax*2;i++) i2index[i] = -1;

	int isup;

	for (i=0;i<_N1LevelCells;i++) {
		isup = _Index1LevelCells[i];
		for (j=_sprnds[i];j<_sprnds[i+1];j++) {
			ihcell = _Real2Index[j*2];
			icell = _Real2Index[j*2+1];
			ibs = ibsblk[ihcell]+icell;
			i2index[ibs*2] = isup;
			i2index[ibs*2+1] = j-_sprnds[i];
		};
	};

// Create the lists of pairs

	int *imask, *imask1;
	int *nlistarr;
	int **plists;

	imask = new int [nimax];
	if (!imask) MemoryFail (funcname);
	imask1 = new int [nimax];
	if (!imask1) MemoryFail (funcname);
	nlistarr = new int [nproc];
	if (!nlistarr) MemoryFail (funcname);
	plists = new int * [nproc];
	if (!plists) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nimax;i++) imask[i] = icycle;

	int nlistloc, ind;

	for (iproc=0;iproc<nproc;iproc++) {
		icycle++;
		nlistloc = 0;
		for (i=iacpu[iproc];i<iacpu[iproc+1];i++) {
			ind = jaindcpu[i];
			ihcell = _JAReal2Index[ind*2];
			icell = _JAReal2Index[ind*2+1];
			ibs = ibsblk[ihcell]+icell;
			if (imask[ibs] != icycle) {
				nlistloc++;
				imask[ibs] = icycle;
			};
		};
		nlistarr[iproc] = nlistloc;
		plists[iproc] = new int [4*nlistloc];
		if (!plists[iproc]) MemoryFail (funcname);
		icycle++;
		nlistloc = 0;
		for (i=iacpu[iproc];i<iacpu[iproc+1];i++) {
			ind = jaindcpu[i];
			ihcell = _JAReal2Index[ind*2];
			icell = _JAReal2Index[ind*2+1];
			ibs = ibsblk[ihcell]+icell;
			if (imask[ibs] != icycle) {
				plists[iproc][nlistloc*2] = ihcell;
				plists[iproc][nlistloc*2+1] = icell;
				nlistloc++;
				imask[ibs] = icycle;
			};
		};
	};

// Exchange the lists

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nproc;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char* [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 2*nlistarr[i]*sizeof(int);
		ObjSend[i] = (char *) plists[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CDecomp::TransformSparsity: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int ni;
	int ihcellnew, icellnew;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			ihcell = piarr[j*2];
			icell = piarr[j*2+1];
			ibs = ibsblk[ihcell]+icell;
			ihcellnew = i2index[ibs*2];
			icellnew = i2index[ibs*2+1];
			if (ihcellnew < 0) throw " CDecomp::TransformSparsity: Index not found";
			piarr[j*2] = ihcellnew;
			piarr[j*2+1] = icellnew;
		};
	};

// Exchange the lists back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CDecomp::TransformSparsity: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Register received data

	for (i=0;i<NObjSend;i++) {
		iproc = CpuIDSend[i];
		ni = ObjSizeSend[i] / (2*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ihcell = piarr[j*2];
			icell = piarr[j*2+1];
			plists[iproc][ni*2+j*2] = ihcell;
			plists[iproc][ni*2+j*2+1] = icell;
		};
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Modify the indices and pack resulting data into CSMatrix structure

	int nlisttot = ncellsloc;
	int nzjatot = nzjaloc;

	CSMatrix amatr (nlisttot,nlisttot,nzjatot,nzjatot);

	int *plist = amatr.GetList ();
	int *plist2 = amatr.GetList2 ();
	int *pia = amatr.GetIa ();
	int *pja = amatr.GetJa ();
	int *pja2 = amatr.GetJa2 ();

	for (i=0;i<=ncellsloc;i++) pia[i] = _IARealCells[i];

	int nlistproc, ihcell1, icell1;

	for (iproc=0;iproc<nproc;iproc++) {
		nlistproc = nlistarr[iproc];
		for (i=0;i<nlistproc;i++) {
			ihcell = plists[iproc][i*2];
			icell = plists[iproc][i*2+1];
			ihcell1 = plists[iproc][nlistproc*2+i*2];
			icell1 = plists[iproc][nlistproc*2+i*2+1];
			ibs = ibsblk[ihcell]+icell;
			imask[ibs] = ihcell1;
			imask1[ibs] = icell1;
		};
		if (iproc == myid) {
			for (i=0;i<ncellsloc;i++) {
				ihcell = _Real2Index[i*2];
				icell = _Real2Index[i*2+1];
				ibs = ibsblk[ihcell]+icell;
				ihcell1 = imask[ibs];
				icell1 = imask1[ibs];
				plist[i] = icell1;
				plist2[i] = ihcell1;
			};
		};
		for (i=iacpu[iproc];i<iacpu[iproc+1];i++) {
			ind = jaindcpu[i];
			ihcell = _JAReal2Index[ind*2];
			icell = _JAReal2Index[ind*2+1];
			ibs = ibsblk[ihcell]+icell;
			ihcell1 = imask[ibs];
			icell1 = imask1[ibs];
			pja[ind] = icell1;
			pja2[ind] = ihcell1;
		};
	};

// Free work arrays

	delete [] blks2cpu;
	delete [] blksini;
	delete [] iacpu;
	delete [] iptrcpu;
	delete [] jaindcpu;
	delete [] ibsblk;
	delete [] i2index;
	delete [] imask;
	delete [] imask1;
	delete [] nlistarr;

	for (i=0;i<nproc;i++) delete [] plists[i];

	delete [] plists;

// Sort matrix data

	amatr.SortListAndColumns2Index ();

	return amatr;

};

// Author: Kharchenko S.A.
// Description: Compute supernodes ordering and blocks partitioning
// CDecomp::SpPartitioning()
//========================================================================================
void CDecomp::SpPartitioning (CMPIComm &_comm, // Compute supernodes ordering and blocks partitioning
												int _nlistsp, int *_listsp, int *_sprndssp,
												int *_ordernd, int _nblks, int *_blksnd, int *_blks2cpu,
												int *_ordersp, int *_blkssp) {

	const char *funcname = "SpPartitioning";

	int myid = _comm.GetMyid ();
	int nproc = _comm.GetNproc ();

// Compute block number for each node element by binary search

	int nloc = _sprndssp[_nlistsp];

	int *nd2blk;

	nd2blk = new int [nloc];
	if (!nd2blk) MemoryFail (funcname);

	int iblk, ind, j, ii, i, jj, iblkbeg, iblkend;

	for (i=0;i<nloc;i++) {
		jj = _ordernd[i];
		iblkbeg = 0;
		iblkend = _nblks-1;
		while (iblkbeg != iblkend) {
			iblk = (iblkbeg+iblkend) / 2;
			if (jj >= _blksnd[iblk] && jj < _blksnd[iblk+1]) {
				iblkbeg = iblk;
				iblkend = iblk;
			} else {
				if (jj >= _blksnd[iblk+1]) iblkbeg = iblk+1;
				if (jj < _blksnd[iblk]) iblkend = iblk-1;
			};
		};
		nd2blk[i] = iblkbeg;
	};

// Compute block number and index for each supernode

	int *sp2blk, *sp2ind;

	sp2blk = new int [_nlistsp];
	if (!sp2blk) MemoryFail (funcname);
	sp2ind = new int [_nlistsp];
	if (!sp2ind) MemoryFail (funcname);

	int jblk;

	for (i=0;i<_nlistsp;i++) {
		ind = _sprndssp[i];
		ii = -1;
		iblk = nd2blk[ind];
		for (j=_sprndssp[i]+1;j<_sprndssp[i+1];j++) {
			jblk = nd2blk[j];
			if (jblk > iblk) iblk = jblk;
		};
		for (j=_sprndssp[i];j<_sprndssp[i+1];j++) {
			jblk = nd2blk[j];
			jj = _ordernd[j];
			if (jblk == iblk) {
				if (ii == -1) {
					ii = jj;
				} else {
					if (jj < ii) ii = jj;
				};
			};
		};
		sp2blk[i] = iblk;
		sp2ind[i] = ii;
	};

// Compute the new blocks in terms of supernodes

	for (i=0;i<=_nblks;i++) _blkssp[i] = 0;

	for (i=0;i<_nlistsp;i++) {
		iblk = sp2blk[i];
		_blkssp[iblk+1]++;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblks+1, _blkssp, _blkssp);

	for (i=0;i<_nblks;i++) _blkssp[i+1] = _blkssp[i]+_blkssp[i+1];

// Prepare for sends the lists

	int *iacpu;
	int *iptrcpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iacpu[i] = 0;

	int iproc;

	for (i=0;i<_nlistsp;i++) {
		iblk = sp2blk[i];
		iproc = _blks2cpu[iblk] % nproc;
		iacpu[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	int **plists;

	plists = new int * [nproc];
	if (!plists) MemoryFail (funcname);

	int ni;

	for (i=0;i<nproc;i++) {
		ni = iacpu[i+1] - iacpu[i];
		plists[i] = new int [ni*2];
		if (!plists[i]) MemoryFail (funcname);
	};

	int k;

	for (i=0;i<_nlistsp;i++) {
		iblk = sp2blk[i];
		ii = sp2ind[i];
		iproc = _blks2cpu[iblk] % nproc;
		k = iptrcpu[iproc]-iacpu[iproc];
		plists[iproc][k*2] = iblk;
		plists[iproc][k*2+1] = ii;
		iptrcpu[iproc]++;
	};

// Exchange the lists

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nproc;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char* [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 2*(iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) plists[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CDecomp::SpPartitioning: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Create local orderings and prepare reply

	int *ibsblk;
	int *iptrblk;

	ibsblk = new int [_nblks];
	if(!ibsblk) MemoryFail(funcname);
	iptrblk = new int [_nblks];
	if(!iptrblk) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) ibsblk[i] = -1;

	int nsupmax = 0;

	int nsupcpu = 0;

	for (i=0;i<_nblks;i++) {
		iproc = _blks2cpu[i];
		if (iproc%nproc == myid) {
			ibsblk[i] = nsupcpu;
			iptrblk[i] = nsupcpu;
			ni = _blkssp[i+1]-_blkssp[i];
			if (ni > nsupmax) nsupmax = ni;
			nsupcpu += ni;
		};
	};

	int *isuparr;

	isuparr = new int [nsupcpu*4];
	if(!isuparr) MemoryFail(funcname);

	CIntInt *iiarr;

	iiarr = new CIntInt [nsupmax];
	if(!iiarr) MemoryFail(funcname);

	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		iproc = CpuIDRecv[i];
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			iblk = piarr[j*2];
			ii = piarr[j*2+1];
			k = iptrblk[iblk];
			isuparr[k*4] = iblk;
			isuparr[k*4+1] = ii;
			isuparr[k*4+2] = j;
			isuparr[k*4+3] = iproc;
			iptrblk[iblk]++;
		};
	};

	int ibs, ind0, indold;

	for (i=0;i<_nblks;i++) {
		iproc = _blks2cpu[i];
		if (iproc%nproc == myid) {
			ni = _blkssp[i+1]-_blkssp[i];
			ibs = ibsblk[i];
			for (j=_blkssp[i];j<_blkssp[i+1];j++) {
				ind0 = j-_blkssp[i];
				ind = ibs+ind0;
				iiarr[ind0].intvalue = isuparr[ind*4+1];
				iiarr[ind0].int2value = ind0;
			};
			qsort (iiarr, ni, sizeof(CIntInt),CIntInt::CompareIntInt);
			for (j=_blkssp[i];j<_blkssp[i+1];j++) {
				ind0 = j-_blkssp[i];
//				ind = ibs+ind0;
				ind = _blkssp[i]+ind0;
				indold = ibs+iiarr[ind0].int2value;
				isuparr[indold*4+1] = ind;
			};
		};
	};

	for (i=0;i<_nblks;i++) iptrblk[i] = ibsblk[i];

	int *icpu2obj;

	icpu2obj = new int [nproc];
	if(!icpu2obj) MemoryFail(funcname);

	for (i=0;i<NObjRecv;i++) {
		iproc = CpuIDRecv[i];
		icpu2obj[iproc] = i;
	};

	int isupnew, iobj;

	for (i=0;i<nsupcpu;i++) {
		iblk    = isuparr[i*4];
		isupnew = isuparr[i*4+1];
		j       = isuparr[i*4+2];
		iproc   = isuparr[i*4+3];
		iobj = icpu2obj[iproc];
		piarr = (int *) ObjRecv[iobj];
		piarr[j*2+1] = isupnew;
	};

// Exchange the lists back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CDecomp::SpPartitioning: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Register received data

	for (i=0;i<NObjSend;i++) {
		iproc = CpuIDSend[i];
		ni = ObjSizeSend[i] / (2*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			iblk = piarr[j*2];
			isupnew = piarr[j*2+1];
			plists[iproc][j*2] = iblk;
			plists[iproc][j*2+1] = isupnew;
		};
	};

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Modify the indices

	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	for (i=0;i<_nlistsp;i++) {
		iblk = sp2blk[i];
		iproc = _blks2cpu[iblk] % nproc;
		k = iptrcpu[iproc]-iacpu[iproc];
		isupnew = plists[iproc][k*2+1];
		_ordersp[i] = isupnew;
//		_ordersp[i] = isupnew+_blkssp[iblk];
		iptrcpu[iproc]++;
	};

// Free work arrays

	delete [] nd2blk;
	delete [] sp2blk;
	delete [] sp2ind;
	delete [] iacpu;
	delete [] iptrcpu;
	for (i=0;i<nproc;i++) delete [] plists[i];
	delete [] plists;
	delete [] ibsblk;
	delete [] iptrblk;
	delete [] isuparr;
	delete [] iiarr;
	delete [] icpu2obj;

};

// Author: Kharchenko S.A.
// Description: Compute distributions of the cells over CPU's for regular 3D mesh
// CDecomp::DistributeRegular3DGrid()
//========================================================================================
int CDecomp::DistributeRegular3DGrid (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
													int _Nx, int _Ny, int _Nz,
													CSolverPar &_Param, 
													CCpuTree &_CpuTree,
													int &_NHyperCells, int *&_HyperCell2CPU,
													int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
													int *&_Index1LevelCellsNew) {

	const char *funcname = "DistributeRegular3DGrid";

	int myid = _comm.GetMyid ();

// Compute sparsity of the collapsed matrix

	CSMatrix axyz;

	if (myid == 0) {
		axyz = CSMatrix::Sparsity3D (_Nz, _Ny, _Nx);
	} else {
		axyz = CSMatrix::Sparsity3D (0, 0, 0);
	};

// Prepare the data

	int ntot = axyz.GetN ();
	int nztot = axyz.GetNzja ();
	int *ialoc = axyz.GetIa ();
	int *jaloc = axyz.GetJa ();

	int *index1LevelCells;
	int *ia1Level2Real;
	int *real2Index;
	int *iaRealCells;
	int *jaReal2Index;

	index1LevelCells = new int [ntot*2];
	if(!index1LevelCells) MemoryFail(funcname);
	ia1Level2Real = new int [ntot+1];
	if(!ia1Level2Real) MemoryFail(funcname);
	real2Index = new int [ntot*2];
	if(!real2Index) MemoryFail(funcname);
	iaRealCells = new int [ntot+1];
	if(!iaRealCells) MemoryFail(funcname);
	jaReal2Index = new int [nztot*2];
	if(!jaReal2Index) MemoryFail(funcname);

	ia1Level2Real[0] = 0;
	iaRealCells[0] = 0;

	int i;

	for (i=0;i<ntot;i++) {
		index1LevelCells[i*2] = i;
		index1LevelCells[i*2+1] = i;
		ia1Level2Real[i+1] = i+1;
		real2Index[i*2] = 0;
		real2Index[i*2+1] = i;
		iaRealCells[i+1] = ialoc[i+1];
	};

	for (i=0;i<nztot;i++) {
		jaReal2Index[i*2] = 0;
		jaReal2Index[i*2+1] = jaloc[i];
	};

// Compute cpu's partitioning

	int *ia1Level2RealNew;

	int info;

	info = DistributeCellsOverCPUsNestedDessectionParMetis (_comm, 
																				ntot, index1LevelCells,
																				ia1Level2Real, real2Index,
																				iaRealCells, jaReal2Index,
																				_Param, 
																				_CpuTree,
																				_NHyperCells, _HyperCell2CPU, _HyperCells, 
																				_N1LevelCellsNew, _1LevelCellsNew2HCells,
																				_Index1LevelCellsNew, ia1Level2RealNew);
	if (info != 0) return info;

// Free work arrays

	delete [] index1LevelCells;
	delete [] ia1Level2Real;
	delete [] real2Index;
	delete [] iaRealCells;
	delete [] jaReal2Index;
	delete [] ia1Level2RealNew;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Compute distributions of the cells over CPU's for regular 3D mesh
// CDecomp::DistributeRegular3DGridDirect()
//========================================================================================
int CDecomp::DistributeRegular3DGridDirect (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
													int _Nx, int _Ny, int _Nz,
													CSolverPar &_Param, 
													CCpuTree &_CpuTree,
													int &_NHyperCells, int *&_HyperCell2CPU,
													int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
													int *&_Index1LevelCellsNew) {

	const char *funcname = "DistributeRegular3DGridDirect";

	int myid = _comm.GetMyid ();
	int nproc = _comm.GetNproc ();

// Compute partitioning

	int nyz = _Ny * _Nz;

	int nproc_2 = 1;
	int nlev = 1;

	while (nproc_2*2 <= nproc) {
		nproc_2 *= 2;
		nlev += 1;
	};

	int nxloc = _Nx;
	int nyloc = _Ny;
	int nzloc = _Nz;

	int npartsx = 1;
	int npartsy = 1;
	int npartsz = 1;

	int ilev = 1;

	while (ilev != nlev) {
		if (nxloc >= nyloc && nxloc >= nzloc) {
			nxloc = nxloc / 2;
			if (nxloc < 1) nxloc = 1;
			npartsx *= 2;
		} else if (nyloc >= nxloc && nyloc >= nzloc) {
			nyloc = nyloc / 2;
			if (nyloc < 1) nyloc = 1;
			npartsy *= 2;
		} else if (nzloc >= nxloc && nzloc >= nyloc) {
			nzloc = nzloc / 2;
			if (nzloc < 1) nzloc = 1;
			npartsz *= 2;
		};
		ilev++;
	};

	int *sprndx;
	int *sprndy;
	int *sprndz;

	sprndx = new int [npartsx+1];
	if(!sprndx) MemoryFail(funcname);
	sprndy = new int [npartsy+1];
	if(!sprndy) MemoryFail(funcname);
	sprndz = new int [npartsz+1];
	if(!sprndz) MemoryFail(funcname);

	int i;

	sprndx[0] = 0;
	for (i=0;i<npartsx-1;i++) {
		sprndx[i+1] = sprndx[i]+nxloc;
		if (sprndx[i+1] > _Nx) sprndx[i+1] = _Nx;
	};
	sprndx[npartsx] = _Nx;
	sprndy[0] = 0;
	for (i=0;i<npartsy-1;i++) {
		sprndy[i+1] = sprndy[i]+nyloc;
		if (sprndy[i+1] > _Ny) sprndy[i+1] = _Ny;
	};
	sprndy[npartsy] = _Ny;
	sprndz[0] = 0;
	for (i=0;i<npartsz-1;i++) {
		sprndz[i+1] = sprndz[i]+nzloc;
		if (sprndz[i+1] > _Nz) sprndz[i+1] = _Nz;
	};
	sprndz[npartsz] = _Nz;

// Create the list for current processor

// Determine index number

	int nhcellsini = npartsx*npartsy*npartsz;

	int *hcells;

	hcells = new int [nhcellsini+1];
	if(!hcells) MemoryFail(funcname);

	int indx = -1;
	int indy = -1;
	int indz = -1;

	int iproc = 0;

	hcells[0] = 0;

	int ni, nj, nk, j, k;

	for (i=0;i<npartsx;i++) {
		ni = sprndx[i+1]-sprndx[i];
		for (j=0;j<npartsy;j++) {
			nj = sprndy[j+1]-sprndy[j];
			for (k=0;k<npartsz;k++) {
				nk = sprndz[k+1]-sprndz[k];
				if (iproc == myid) {
					indx = i;
					indy = j;
					indz = k;
				};
				hcells[iproc+1] = hcells[iproc] + ni*nj*nk;
				iproc++;
			};
		};
	};

// Create binary tree

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nproc_2+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nproc_2+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nproc_2;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nproc_2;i++) procweight[i] = 1.0e0;

	CTree tree (nproc_2, 2, procweight, memory2cpu);

	tree.SetMyid (myid);

	delete [] memory2cpu;
	delete [] procweight;

	tree.SortNodes ();

	int nnodesloc = tree.GetNnodes ();

	CNode *pnodes = tree.GetNode ();

	for (i=0;i<nnodesloc;i++) {
		pnodes[i].indbeg = i;
		pnodes[i].indend = i;
		pnodes[i].indbegtot = i;
		pnodes[i].indendtot = i;
	};

// Create HyperCell2CPU and HyperCells arrays

	_NHyperCells = nnodesloc;

	_HyperCell2CPU = new int [_NHyperCells];
	if(!_HyperCell2CPU) MemoryFail(funcname);
	_HyperCells = new int [_NHyperCells+1];
	if(!_HyperCells) MemoryFail(funcname);

	iproc = 0;

	int nchildsloc;

	_HyperCells[0] = 0;

	for (i=0;i<nnodesloc;i++) {
		_HyperCell2CPU[i] = pnodes[i].GetNodecpu ();
		nchildsloc = pnodes[i].GetNchilds ();
		if (nchildsloc == 1) {
			_HyperCells[i+1] = _HyperCells[i] + hcells[iproc+1]-hcells[iproc];
			iproc++;
		} else {
			_HyperCells[i+1] = _HyperCells[i];
		};
	};

	_N1LevelCellsNew = 0;
	if (indx >= 0) {
		_N1LevelCellsNew = hcells[myid+1]-hcells[myid];
	};

	_1LevelCellsNew2HCells = new int [_N1LevelCellsNew];
	if (!_1LevelCellsNew2HCells) MemoryFail (funcname);
	_Index1LevelCellsNew = new int [_N1LevelCellsNew*2];
	if (!_Index1LevelCellsNew) MemoryFail (funcname);

	int ihcell = -1;

	if (indx >= 0) {
		int *pcpuidend = tree.GetCpuidend ();
		ihcell = pcpuidend[myid];
	};

	for (i=0;i<_N1LevelCellsNew;i++) _1LevelCellsNew2HCells[i] = ihcell;

	if (indx >= 0) {
		int nz = 0;
		int ind;
		for (i=sprndx[indx];i<sprndx[indx+1];i++) {
			for (j=sprndy[indy];j<sprndy[indy+1];j++) {
				for (k=sprndz[indz];k<sprndz[indz+1];k++) {
					ind = i*nyz+j*_Nz+k;
					_Index1LevelCellsNew[nz*2] = ind;
					_Index1LevelCellsNew[nz*2+1] = ind;
					nz++;
				};
			};
		};
	};

// Modify the tree

	_CpuTree.SetTree (tree);
	_CpuTree.SetTreeExt (tree);
	_CpuTree.SetTreeHc (tree);

// Free work arrays

	delete [] sprndx;
	delete [] sprndy;
	delete [] sprndz;
	delete [] hcells;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the cells over CPU's for regular 3D mesh
// CDecomp::DistributeRegular3DGridWeights()
//========================================================================================
int CDecomp::DistributeRegular3DGridWeights (CMPIComm &_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
													int _Nx, int _Ny, int _Nz,
													CSolverPar &_Param, 
													CCpuTree &_CpuTree,
													int &_NHyperCells, int *&_HyperCell2CPU,
													int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
													int *&_Index1LevelCellsNew) {

	const char *funcname = "DistributeRegular3DGridWeights";

// Compute initial distribution

	int nhcells, *hcell2cpu, *hcells;

	DistributeRegular3DGridDirect (_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
												_Nx, _Ny, _Nz,
												_Param, 
												_CpuTree,
												nhcells, hcell2cpu, hcells, 
												_N1LevelCellsNew, _1LevelCellsNew2HCells, _Index1LevelCellsNew);

// Perform modification of the hcells according to the weights

	CSlvParam *pparam = _Param.GetParam ();

	int weightcellvolumeloc = pparam->weightcellvolume;

	int *hco2n;

	hco2n = new int [nhcells+1];
	if (!hco2n) MemoryFail (funcname);

	hco2n[0] = 0;

	_NHyperCells = 0;

	int i, j, nhi, ni, nj;

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		nhi = (ni+weightcellvolumeloc-1) / weightcellvolumeloc;
		if (nhi < 1) nhi = 1;
		_NHyperCells += nhi;
		hco2n[i+1] = _NHyperCells;
	};

	_HyperCell2CPU = new int [_NHyperCells];
	if(!_HyperCell2CPU) MemoryFail(funcname);
	_HyperCells = new int [_NHyperCells+1];
	if(!_HyperCells) MemoryFail(funcname);

	_HyperCells[0] = 0;

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		nhi = (ni+weightcellvolumeloc-1) / weightcellvolumeloc;
		if (nhi < 1) nhi = 1;
		for (j=hco2n[i];j<hco2n[i+1];j++) {
			_HyperCell2CPU[j] = hcell2cpu[i];
		};
		nj = ni / nhi;
		for (j=hco2n[i];j<hco2n[i+1]-1;j++) {
			_HyperCells[j+1] = _HyperCells[j]+nj;
		};
		_HyperCells[hco2n[i+1]] = hcells[i+1];
	};

// Modify array 1LevelCellsNew2HCells

	int ihcell, ibeg, iend, ibegn, ihcellnew;

	int index = 0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		if (iend-ibeg+1 != hcells[ihcell+1]-hcells[ihcell]) throw " CDecomp::DistributeRegular3DGridWeights: error in array 1LevelCellsNew2HCells";
		ibegn = hcells[ihcell];
		for (ihcellnew=hco2n[ihcell];ihcellnew<hco2n[ihcell+1];ihcellnew++) {
			for (i=_HyperCells[ihcellnew];i<_HyperCells[ihcellnew+1];i++) {
				_1LevelCellsNew2HCells[i-ibegn+ibeg] = ihcellnew;
			};
		};
		index = iend+1;
	};

// Register new hypercells in the tree's

	CTree *ptree = _CpuTree.GetTree ();
	CTree *ptreeext = _CpuTree.GetTreeext ();
	CTree *ptreeexthc = _CpuTree.GetTreeHc ();

	ptree->ReassignBlockNumbers (hco2n);
	ptreeext->ReassignBlockNumbers (hco2n);
	ptreeexthc->ReassignBlockNumbers (hco2n);

// Free work arrays

	delete [] hcell2cpu;
	delete [] hcells;
	delete [] hco2n;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
// CDecomp::DistributeCellsOverCPUsNestedDessectionParMetisWeights()
//========================================================================================
int CDecomp::DistributeCellsOverCPUsNestedDessectionParMetisWeights (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																					int _N1LevelCells, int *_Index1LevelCells,
																					int *_IA1Level2Real, int *_Real2Index,
																					int *_IARealCells, int *_JAReal2Index,
																					CSolverPar &_Param, 
																					CCpuTree &_CpuTree,
																					int &_NHyperCells, int *&_HyperCell2CPU,
																					int *&_HyperCells, int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
																					int *&_Index1LevelCellsNew) {

	const char *funcname = "DistributeCellsOverCPUsNestedDessectionParMetisWeights";

	CSlvParam *pparam = _Param.GetParam ();

	int weightcellvolumeloc = pparam->weightcellvolume;

// Compute initial distribution

	int nhcells, *hcell2cpu, *hcells, *ia1Level2RealNew;

	DistributeCellsOverCPUsNestedDessectionParMetis (_comm, // Compute distributions of the cells over CPU's for regular 3D mesh
																		_N1LevelCells, _Index1LevelCells,
																		_IA1Level2Real, _Real2Index,
																		_IARealCells, _JAReal2Index,
																		_Param, 
																		_CpuTree,
																		nhcells, hcell2cpu, hcells, 
																		_N1LevelCellsNew, _1LevelCellsNew2HCells, _Index1LevelCellsNew, ia1Level2RealNew);

// Count new diff array

	int *nRealCellsNewIn1Level;

	nRealCellsNewIn1Level = new int [_N1LevelCellsNew];
	if (!nRealCellsNewIn1Level) MemoryFail (funcname);

	int i;

	for (i=0;i<_N1LevelCellsNew;i++) nRealCellsNewIn1Level[i] = ia1Level2RealNew[i+1]-ia1Level2RealNew[i];

// Perform repartitioning of the cells

	int *hcweights;

	hcweights = new int [nhcells];
	if (!hcweights) MemoryFail (funcname);

	int ihcell;

	for (i=0;i<nhcells;i++) hcweights[i] = 0;

	for (i=0;i<_N1LevelCellsNew;i++) {
		ihcell = _1LevelCellsNew2HCells[i];
		hcweights[ihcell] += ia1Level2RealNew[i+1]-ia1Level2RealNew[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nhcells, hcweights, hcweights);

	int *hco2n;

	hco2n = new int [nhcells+1];
	if (!hco2n) MemoryFail (funcname);

	for (i=0;i<=nhcells;i++) hco2n[i] = 0;

	int j;

	int ibeg, iend, icell, icell1, icount, icount1, nhcellsloc;

	int index = 0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		nhcellsloc = 0;
		icell = ibeg;
		while (icell <= iend) {
			icount = ia1Level2RealNew[icell+1]-ia1Level2RealNew[icell];
			if (icount < weightcellvolumeloc) {
				while (icell+1 <= iend) {
					icell1 = icell+1;
					icount1 = icount + ia1Level2RealNew[icell1+1]-ia1Level2RealNew[icell1];
					if (icount1 <= weightcellvolumeloc) {
						icount = icount1;
						icell = icell1;
					} else {
						break;
					};
				};
			};
			nhcellsloc++;
			icell++;
		};
		hco2n[ihcell+1] = nhcellsloc;
		index = iend+1;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nhcells+1, hco2n, hco2n);

	for (i=0;i<nhcells;i++) hco2n[i+1] = hco2n[i]+hco2n[i+1];

	_NHyperCells = hco2n[nhcells];

	_HyperCell2CPU = new int [_NHyperCells];
	if(!_HyperCell2CPU) MemoryFail(funcname);

	for (i=0;i<nhcells;i++) {
		for (j=hco2n[i];j<hco2n[i+1];j++) {
			_HyperCell2CPU[j] = hcell2cpu[i];
		};
	};

	_HyperCells = new int [_NHyperCells+1];
	if(!_HyperCells) MemoryFail(funcname);

	_HyperCells[0] = 0;

	for (i=0;i<=_NHyperCells;i++) _HyperCells[i] = 0;

	index = 0;

	int ibegnew, iendnew, ihcellnew, ibeg0, iend0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		icell = ibeg;
		ibeg0 = ibeg;
		iend0 = ibeg-1;
		ibegnew = icell;
		iendnew = icell;
		icount = ia1Level2RealNew[icell+1]-ia1Level2RealNew[icell];
		ihcellnew = hco2n[ihcell];
		while (ihcellnew < hco2n[ihcell+1]-1) {
			while (icell+1 <= iend) {
				icell1 = icell+1;
				icount1 = icount + ia1Level2RealNew[icell1+1]-ia1Level2RealNew[icell1];
				if (icount1 <= weightcellvolumeloc) {
					icount = icount1;
					icell = icell1;
					iendnew = icell;
				} else {
					iendnew = icell;
					break;
				};
			};
			_HyperCells[ihcellnew+1] = iendnew-ibegnew+1;
			ihcellnew++;
			ibegnew = icell+1;
			iendnew = icell+1;
			iend0 = icell;
			icell = ibegnew;
			icount = ia1Level2RealNew[icell+1]-ia1Level2RealNew[icell];
		};
		_HyperCells[hco2n[ihcell+1]] = (hcells[ihcell+1]-hcells[ihcell])-(iend0-ibeg0+1);
		index = iend+1;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells+1, _HyperCells, _HyperCells);

	for (i=0;i<_NHyperCells;i++) _HyperCells[i+1] = _HyperCells[i]+_HyperCells[i+1];

// Modify array 1LevelCellsNew2HCells

	int ibegn;

	index = 0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		if (iend-ibeg+1 != hcells[ihcell+1]-hcells[ihcell]) throw " CDecomp::DistributeCellsOverCPUsNestedDessectionParMetisWeights: error in array 1LevelCellsNew2HCells";
		ibegn = hcells[ihcell];
		for (ihcellnew=hco2n[ihcell];ihcellnew<hco2n[ihcell+1];ihcellnew++) {
			for (i=_HyperCells[ihcellnew];i<_HyperCells[ihcellnew+1];i++) {
				_1LevelCellsNew2HCells[i-ibegn+ibeg] = ihcellnew;
			};
		};
		index = iend+1;
	};

// Count weights of the new hypercells

	int *hcweightsnew;

	hcweightsnew = new int [_NHyperCells];
	if (!hcweightsnew) MemoryFail (funcname);

	for (i=0;i<_NHyperCells;i++) hcweightsnew[i] = 0;

	for (i=0;i<_N1LevelCellsNew;i++) {
		ihcell = _1LevelCellsNew2HCells[i];
		hcweightsnew[ihcell] += ia1Level2RealNew[i+1]-ia1Level2RealNew[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells, hcweightsnew, hcweightsnew);

// Register new hypercells in the tree's

	CTree *ptree = _CpuTree.GetTree ();
	CTree *ptreeext = _CpuTree.GetTreeext ();
	CTree *ptreeexthc = _CpuTree.GetTreeHc ();

	ptree->ReassignBlockNumbers (hco2n);
	ptreeext->ReassignBlockNumbers (hco2n);
	ptreeexthc->ReassignBlockNumbers (hco2n);

// Free work arrays

	delete [] hcell2cpu;
	delete [] hcells;
	delete [] ia1Level2RealNew;
	delete [] nRealCellsNewIn1Level;
	delete [] hcweights;
	delete [] hco2n;
	delete [] hcweightsnew;

	return 0;

};

// Author: Kharchenko S.A.
// Description: Compute distributions of the 1 level cells via ParMetis via real cells sparsity, inter CPU's boundary is on the 1 level cells boundary
// CDecomp::DistributeCellsNDParMetisBoundary1Level()
//========================================================================================
int CDecomp::DistributeCellsNDParMetisBoundary1Level (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																			int _N1LevelCells, int *_Index1LevelCells,
																			int *_IA1Level2Real, int *_Real2Index,
																			int *_IARealCells, int *_JAReal2Index,
																			CSolverPar &_Param, 
																			CCpuTree &_CpuTree,
																			int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																			int &_NHyperCellsReal, int *&_HyperCell2CPUReal, int *&_HyperCellsInRealCells,
																			int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells, int *&_IA1Level2RealNew,
																			int *&_Index1LevelCellsNew, int *&_IndexRealCellsNew, int *&_IndexRealCellsIni) {

	const char *funcname = "DistributeCellsNDParMetisBoundary1Level";

	int myid = _comm.GetMyid ();
	int nproc = _comm.GetNproc ();

// Create local binary tree

	int nprocext = _Param.GetNcpuext ();

	if (nprocext <= 0) {
		nprocext = 1;
		while (nprocext*2 <= nproc) nprocext *= 2;
	};

	if (nprocext > nproc) {
		nprocext = 1;
		while (nprocext*2 <= nproc) nprocext *= 2;
	};

	int nchilds = 2;

	double *perf2cpu;
	double *memory2cpu;

	perf2cpu = new double [nprocext];
	if (!perf2cpu) MemoryFail (funcname);
	memory2cpu = new double [nprocext];
	if (!memory2cpu) MemoryFail (funcname);

	int i;

	for (i=0;i<nprocext;i++) {
		perf2cpu[i] = 1.0e0;
		memory2cpu[i] = 1.0e0;
	};

	CTree treebinary (nprocext,nchilds,perf2cpu,memory2cpu);
	treebinary.SortNodes ();

	delete [] perf2cpu;
	delete [] memory2cpu;

	treebinary.CpuListsForNodesSchur ();

// Sort the indices of 1 level cells

	CIntInt *iiarr;

	iiarr = new CIntInt [_N1LevelCells];
	if (!iiarr) MemoryFail (funcname);

	for (i=0;i<_N1LevelCells;i++) {
		iiarr[i].intvalue = _Index1LevelCells[i*2];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, _N1LevelCells, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *listindex, *attrindex, *ordindexinv;

	listindex = new int [_N1LevelCells];
	if (!listindex) MemoryFail (funcname);
	attrindex = new int [_N1LevelCells];
	if (!attrindex) MemoryFail (funcname);
	ordindexinv = new int [_N1LevelCells];
	if (!ordindexinv) MemoryFail (funcname);

	for (i=0;i<_N1LevelCells;i++) {
		listindex[i] = iiarr[i].intvalue;
		ordindexinv[i] = iiarr[i].int2value;
	};

	delete [] iiarr;

// Check that array Index1LevelCells is correct

	int indexmax = 0;
	if (_N1LevelCells > 0) indexmax = listindex[_N1LevelCells-1];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, MAXIMUM, 1, &indexmax, &indexmax);

	int ncellstot = _N1LevelCells;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &ncellstot, &ncellstot);

	if (indexmax != ncellstot-1) throw " CDecomp::DistributeCellsNDParMetisBoundary1Level: array Index1LevelCells on entry is incorrect ";

// Create ordered sprnds

	int *sprnds;
	int *sprndsord;

	sprnds = new int [_N1LevelCells+1];
	if (!sprnds) MemoryFail (funcname);
	sprndsord = new int [_N1LevelCells+1];
	if (!sprndsord) MemoryFail (funcname);

	sprnds[0] = 0;
	sprndsord[0] = 0;

	int ind;

	for (i=0;i<_N1LevelCells;i++) {
		ind = ordindexinv[i];
		sprnds[i+1] = sprnds[i] + (_IA1Level2Real[i+1]-_IA1Level2Real[i]);
		sprndsord[i+1] = sprndsord[i] + (_IA1Level2Real[ind+1]-_IA1Level2Real[ind]);
	};

// Compute the total size of the matrix

	int ntotmtr = sprndsord[_N1LevelCells];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &ntotmtr, &ntotmtr);

// Partition the indices into the parts for balanced load

	CSlvParam *pparam = _Param.GetParam ();

	int nprocpart = pparam->ncpupart;

	if (nprocpart <= 0) {
		nprocpart = 1;
		while (nprocpart*2 <= nproc) nprocpart *= 2;
	};

	if (nprocpart > nproc) {
		nprocpart = 1;
		while (nprocpart*2 <= nproc) nprocpart *= 2;
	};

	nprocpart = 1;

	int *blksindex;

	blksindex = new int [nprocpart+1];
	if (!blksindex) MemoryFail (funcname);

	PartitionIndices (_comm, 
							nprocpart, _N1LevelCells, listindex, sprndsord,
							blksindex);

// Create CSMatrix 2index structure that is based on pair {1levelcell,realcell} instead of {hypercell,realcell1}

	CSMatrix mtr;

	for (i=0;i<_N1LevelCells;i++) {
		attrindex[i] = _Index1LevelCells[2*i];
	};

	mtr = TransformSparsity (_comm, 
										_N1LevelCells, sprnds, attrindex,
										_Real2Index,
										_IARealCells, _JAReal2Index);

	mtr.SetM (ntotmtr);
	mtr.SetN (ntotmtr);
	mtr.SetNsupr (ntotmtr);
	mtr.SetNsupcBase (ntotmtr);

// Exchange matrix data

	int nint = 0;
	int ibeg = 0;
	int iend = -1;

	if (myid < nprocpart) {
		nint = 1;
		ibeg = blksindex[myid];
		iend = blksindex[myid+1]-1;
	};

	int nsuploc = iend-ibeg+1;

	int *listsp;
	int *attrinterval;

	listsp = new int [nsuploc];
	if (!listsp) MemoryFail (funcname);
	attrinterval = new int [nsuploc];
	if (!attrinterval) MemoryFail (funcname);

	for (i=0;i<nsuploc;i++) listsp[i] = ibeg+i;

	CSMatrix mtrnew;

	mtrnew = mtr.GetSubmatrixViaIntervals2Index (_comm,
																nint, &ibeg, &iend);

	CSMatrix mtrdummy;

	mtr = mtrdummy;

// Compute global attributes array for the condenced interval [ibeg,iend]

	int *index2cpuini;

	index2cpuini = new int [_N1LevelCells];
	if (!index2cpuini) MemoryFail (funcname);

	int *iacpu;
	int *iptrcpu;
	int **plists;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);
	plists = new int * [nproc];
	if (!plists) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iacpu[i] = 0;

	int iproc, ni, icpubeg, icpuend, icpu;

	for (i=0;i<_N1LevelCells;i++) {
		ind = _Index1LevelCells[i*2];
		icpubeg = 0;
		icpuend = nprocpart-1;
		while (icpubeg != icpuend) {
			icpu = (icpubeg+icpuend)/2;
			if (ind >= blksindex[icpu] && ind < blksindex[icpu+1]) {
				icpubeg = icpu;
				icpuend = icpu;
			} else if (ind < blksindex[icpu]) {
				icpuend = icpu-1;
			} else if (ind >= blksindex[icpu+1]) {
				icpubeg = icpu+1;
			};
		};
		index2cpuini[i] = icpubeg;
		iacpu[icpubeg+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];
	for (i=0;i<nproc;i++) iptrcpu[i] = iacpu[i];

	for (i=0;i<nproc;i++) {
		ni = iacpu[i+1]-iacpu[i];
		plists[i] = new int [ni*2];
		if (!plists[i]) MemoryFail (funcname);
	};

	int attr, k;

	for (i=0;i<_N1LevelCells;i++) {
		ind = _Index1LevelCells[i*2];
		attr = _Index1LevelCells[i*2+1];
		iproc = index2cpuini[i];
		k = iptrcpu[iproc]-iacpu[iproc];
		plists[iproc][k*2] = ind;
		plists[iproc][k*2+1] = attr;
		iptrcpu[iproc]++;
	};

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nproc;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char* [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	for(i = 0; i < NObjSend; i++)
	{
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 2*(iacpu[i+1]-iacpu[i])*sizeof(int);
		ObjSend[i] = (char *) plists[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CDecomp::DistributeCellsNDParMetisBoundary1Level: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	int *piarr;

	int j;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*2];
			attr = piarr[j*2+1];
			attrinterval[ind-ibeg] = attr;
		};
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

	for (i=0;i<nproc;i++) delete [] plists[i];

// Create supernodes partitioning

	int *plist2 = mtrnew.GetList2 ();
	int nlistnew = mtrnew.GetNlist ();

	int *sprndsloc;

	sprndsloc = new int [nsuploc+1];
	if (!sprndsloc) MemoryFail (funcname);

	nsuploc = 0;
	int isupprev = -1;

	sprndsloc[0] = 0;

	if (nlistnew > 0) {
		isupprev = plist2[0];
		listsp[nsuploc] = isupprev;
		nsuploc = 1;
		sprndsloc[nsuploc] = 1;
	};

	int isup;

	for (i=1;i<nlistnew;i++) {
		sprndsloc[nsuploc] = i+1;
		isup = plist2[i];
		if (isup != isupprev) {
			sprndsloc[nsuploc] = i;
			isupprev = isup;
			listsp[nsuploc] = isupprev;
			nsuploc++;
		};
	};

	sprndsloc[nsuploc] = nlistnew;

// Renumber the matrix data (condense pairs numbering into the index numbering)

	mtrnew.Transform2IndexIntoContinuous (_comm);

// Create local communicator

	int *listcpu;

	listcpu = new int [nprocpart];
	if (!listcpu) MemoryFail (funcname);

	for (i=0;i<nprocpart;i++) listcpu[i] = i;

	CMPIComm commloc;

	commloc = _comm.CreateComm (nprocpart, listcpu);

	delete [] listcpu;

// Distribute via ParMetis

	int nsuptot = nsuploc;
	int ntot = sprndsloc[nsuploc];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &nsuptot, &nsuptot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &ntot, &ntot);

	int *ordersp;
	int *orderloc;

	ordersp = new int [nsuptot];
	if (!ordersp) MemoryFail (funcname);
	orderloc = new int [ntot];
	if (!orderloc) MemoryFail (funcname);

	CTree treebinarynd;

	if (myid < nprocpart) {

// Distribute point data

		mtrnew.PartBinaryTreeNDSups (treebinary, treebinarynd, 
												nsuploc, sprndsloc, 
												ordersp, orderloc,
												_NHyperCells, _HyperCells, _HyperCell2CPU, 
												_NHyperCellsReal, _HyperCellsInRealCells, _HyperCell2CPUReal);

		CSMatrix *pablkstr = treebinarynd.GetAblktree ();

		*pablkstr = treebinarynd.MaximalBlockSparsityForSortedTree (_NHyperCellsReal);

		treebinary.StorePartitioning (_NHyperCells, _HyperCells, _HyperCell2CPU);
		treebinarynd.StorePartitioning (_NHyperCellsReal, _HyperCellsInRealCells, _HyperCell2CPUReal);

// Redistribute the boundary if necessary

		CSlvParam *pparam = _Param.GetParam ();

		if (pparam->schurtype > 0) {

			treebinarynd.ReassignBlockNumbers (_HyperCellsInRealCells);

			delete [] _HyperCellsInRealCells;
			delete [] _HyperCell2CPUReal;

			CSMatrix ao = mtrnew.OrdMtr (orderloc);

			ao.PartBinaryTreeNDSchurBnd_FixedInitialOrdering (treebinarynd, orderloc, *pparam,
																				_NHyperCellsReal, _HyperCellsInRealCells, _HyperCell2CPUReal);

		};

		if (myid > 0) delete [] _HyperCells;
		if (myid > 0) delete [] _HyperCell2CPU;
		if (myid > 0) delete [] _HyperCellsInRealCells;
		if (myid > 0) delete [] _HyperCell2CPUReal;

	};

	mtrnew = mtrdummy;

// Exchange tree and blocks information

	int isize = 0;
	char *ptree;

	if (myid == 0) {
		treebinary.PackTree (isize,ptree);
	};
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &isize, &isize);
	if (myid > 0) {
		ptree = new char [isize];
		if (!ptree) MemoryFail (funcname);
		for (i=0;i<isize;i++) ptree[i] = 0;
	};
	CMPIExchange::ExchangeArrayMPI (_comm, CHARVALUE, ADD, isize, ptree, ptree);
	treebinary.UnPackTree (isize,ptree);
	delete [] ptree;

	if (myid == 0) {
		treebinarynd.PackTree (isize,ptree);
	} else {
		isize = 0;
	};
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &isize, &isize);
	if (myid > 0) {
		ptree = new char [isize];
		if (!ptree) MemoryFail (funcname);
		for (i=0;i<isize;i++) ptree[i] = 0;
	};
	CMPIExchange::ExchangeArrayMPI (_comm, CHARVALUE, ADD, isize, ptree, ptree);
	treebinarynd.UnPackTree (isize,ptree);
	delete [] ptree;

	if (myid > 0) {
		_NHyperCells = 0;
	};
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &_NHyperCells, &_NHyperCells);

	if (myid > 0) {
		_NHyperCellsReal = 0;
	};
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 1, &_NHyperCellsReal, &_NHyperCellsReal);

	if (myid > 0) {

		_HyperCells = new int [_NHyperCells+1];
		if (!_HyperCells) MemoryFail (funcname);
		_HyperCell2CPU = new int [_NHyperCells];
		if (!_HyperCell2CPU) MemoryFail (funcname);
		_HyperCellsInRealCells = new int [_NHyperCellsReal+1];
		if (!_HyperCellsInRealCells) MemoryFail (funcname);
		_HyperCell2CPUReal = new int [_NHyperCellsReal];
		if (!_HyperCell2CPUReal) MemoryFail (funcname);

		for (i=0;i<_NHyperCells+1;i++) _HyperCells[i] = 0;
		for (i=0;i<_NHyperCells;i++) _HyperCell2CPU[i] = 0;
		for (i=0;i<_NHyperCellsReal+1;i++) _HyperCellsInRealCells[i] = 0;
		for (i=0;i<_NHyperCellsReal;i++) _HyperCell2CPUReal[i] = 0;

	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells+1, _HyperCells, _HyperCells);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells, _HyperCell2CPU, _HyperCell2CPU);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCellsReal+1, _HyperCellsInRealCells, _HyperCellsInRealCells);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCellsReal, _HyperCell2CPUReal, _HyperCell2CPUReal);

	if (myid >= nprocpart) {
		for (i=0;i<nsuptot;i++) ordersp[i] = 0;
		for (i=0;i<ntot;i++) orderloc[i] = 0;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsuptot, ordersp, ordersp);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, ntot, orderloc, orderloc);

// Store attributes of all 1 level cells

	int *attrtot;

	attrtot = new int [nsuptot];
	if (!attrtot) MemoryFail (funcname);

	for (i=0;i<nsuptot;i++) attrtot[i] = 0;

	int ind1, ind2;

	for (i=0;i<_N1LevelCells;i++) {
		ind1 = _Index1LevelCells[i*2];
		ind2 = _Index1LevelCells[i*2+1];
		attrtot[ind1] = ind2;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsuptot, attrtot, attrtot);

// Create array sprndstot

	int *sprndstot;

	sprndstot = new int [nsuptot+1];
	if (!sprndstot) MemoryFail (funcname);

	for (i=0;i<=nsuptot;i++) sprndstot[i] = 0;

	if (myid == 0) {
		for (i=0;i<=nsuptot;i++) sprndstot[i] = sprndsloc[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsuptot+1, sprndstot, sprndstot);

// Store cpu number for all 1 level cells

	int *cpuinitot;

	cpuinitot = new int [nsuptot];
	if (!cpuinitot) MemoryFail (funcname);

	for (i=0;i<nsuptot;i++) cpuinitot[i] = 0;

	for (i=0;i<_N1LevelCells;i++) {
		ind1 = _Index1LevelCells[i*2];
		cpuinitot[ind1] = myid;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nsuptot, cpuinitot, cpuinitot);

// Prepare return data

	int *iordersp;
	int *iorderloc;

	iordersp = new int [nsuptot];
	if (!iordersp) MemoryFail (funcname);
	iorderloc = new int [ntot];
	if (!iorderloc) MemoryFail (funcname);

	for (i=0;i<nsuptot;i++) iordersp[ordersp[i]] = i;
	for (i=0;i<ntot;i++) iorderloc[orderloc[i]] = i;

	_N1LevelCellsNew = 0;

	for (i=0;i<_NHyperCells;i++) {
		if (_HyperCell2CPU[i] == myid) {
			_N1LevelCellsNew += _HyperCells[i+1]-_HyperCells[i];
		};
	};

	_1LevelCellsNew2HCells = new int [_N1LevelCellsNew];
	if (!_1LevelCellsNew2HCells) MemoryFail (funcname);
	_Index1LevelCellsNew = new int [_N1LevelCellsNew*2];
	if (!_Index1LevelCellsNew) MemoryFail (funcname);
	_IA1Level2RealNew = new int [_N1LevelCellsNew+1];
	if (!_IA1Level2RealNew) MemoryFail (funcname);

	_N1LevelCellsNew = 0;
	_IA1Level2RealNew[0] = 0;
	int nrcells = 0;

	for (i=0;i<_NHyperCells;i++) {
		if (_HyperCell2CPU[i] == myid) {
			for (j=_HyperCells[i];j<_HyperCells[i+1];j++) {
				_1LevelCellsNew2HCells[_N1LevelCellsNew] = i;
				ind = iordersp[j];
				_Index1LevelCellsNew[_N1LevelCellsNew*2] = ind;
				_Index1LevelCellsNew[_N1LevelCellsNew*2+1] = attrtot[ind];
				_IA1Level2RealNew[_N1LevelCellsNew+1] = _IA1Level2RealNew[_N1LevelCellsNew]+(sprndstot[ind+1]-sprndstot[ind]);
				_N1LevelCellsNew++;
				nrcells += sprndstot[ind+1]-sprndstot[ind];
			};
		};
	};

	_IndexRealCellsNew = new int [nrcells];
	if (!_IndexRealCellsNew) MemoryFail (funcname);

	nrcells = 0;

	for (i=0;i<_NHyperCells;i++) {
		if (_HyperCell2CPU[i] == myid) {
			for (j=_HyperCells[i];j<_HyperCells[i+1];j++) {
				ind = iordersp[j];
				for (k=sprndstot[ind];k<sprndstot[ind+1];k++) {
					_IndexRealCellsNew[nrcells] = iorderloc[k];
					nrcells++;
				};
			};
		};
	};

// Create array IndexRealCellsIni

	int ncellsloc = _IA1Level2Real[_N1LevelCells];

	_IndexRealCellsIni = new int [ncellsloc];
	if (!_IndexRealCellsIni) MemoryFail (funcname);

	int nz = 0;

	for (i=0;i<_N1LevelCells;i++) {
		ind = _Index1LevelCells[i*2];
		for (k=sprndstot[ind];k<sprndstot[ind+1];k++) {
			_IndexRealCellsIni[nz] = orderloc[k];
			nz++;
		};
	};

// Modify the tree

	_CpuTree.SetTree (treebinary);
	_CpuTree.SetTreeExt (treebinarynd);
	_CpuTree.SetTreeHc (treebinary);

// Free work arrays

	delete [] listindex;
	delete [] attrindex;
	delete [] ordindexinv;
	delete [] sprnds;
	delete [] sprndsord;
	delete [] blksindex;
	delete [] listsp;
	delete [] attrinterval;
	delete [] index2cpuini;
	delete [] iacpu;
	delete [] iptrcpu;
	delete [] plists;
	delete [] sprndsloc;
	delete [] ordersp;
	delete [] orderloc;
	delete [] attrtot;
	delete [] cpuinitot;
	delete [] iordersp;
	delete [] iorderloc;
	delete [] sprndstot;

	return 0;

};

// Author: Kharchenko S.A. 
// Description: Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
// CDecomp::DistributeCellsNDParMetisBoundary1LevelWeights()
//========================================================================================
int CDecomp::DistributeCellsNDParMetisBoundary1LevelWeights (CMPIComm &_comm, // Compute distributions of the 1 level cells over CPU's via ParMetis according to real cells sparsity
																			int _N1LevelCells, int *_Index1LevelCells,
																			int *_IA1Level2Real, int *_Real2Index,
																			int *_IARealCells, int *_JAReal2Index,
																			CSolverPar &_Param, 
																			CCpuTree &_CpuTree,
																			int &_NHyperCells, int *&_HyperCell2CPU, int *&_HyperCells, 
																			int &_NHyperCellsReal, int *&_HyperCell2CPUReal, int *&_HyperCellsInRealCells,
																			int &_N1LevelCellsNew, int *&_1LevelCellsNew2HCells,
																			int *&_Index1LevelCellsNew, int *&_IndexRealCellsNew, int *&_IndexRealCellsIni) {

	const char *funcname = "DistributeCellsNDParMetisBoundary1LevelWeights";

	int myid = _comm.GetMyid ();

	CSlvParam *pparam = _Param.GetParam ();

	int weightcellvolumeloc = pparam->weightcellvolume;

// Compute initial distribution

	int nhcells, *hcell2cpu, *hcells, *ia1Level2RealNew;

	DistributeCellsNDParMetisBoundary1Level (_comm, 
															_N1LevelCells, _Index1LevelCells,
															_IA1Level2Real, _Real2Index,
															_IARealCells, _JAReal2Index,
															_Param, 
															_CpuTree,
															nhcells, hcell2cpu, hcells, 
															_NHyperCellsReal, _HyperCell2CPUReal, _HyperCellsInRealCells,
															_N1LevelCellsNew, _1LevelCellsNew2HCells, ia1Level2RealNew, _Index1LevelCellsNew, 
															_IndexRealCellsNew, _IndexRealCellsIni);

// Count new diff array

	int *nRealCellsNewIn1Level;

	nRealCellsNewIn1Level = new int [_N1LevelCellsNew];
	if (!nRealCellsNewIn1Level) MemoryFail (funcname);

	int i;

	for (i=0;i<_N1LevelCellsNew;i++) nRealCellsNewIn1Level[i] = ia1Level2RealNew[i+1]-ia1Level2RealNew[i];

// Perform repartitioning of the cells

	int *hcweights;

	hcweights = new int [nhcells];
	if (!hcweights) MemoryFail (funcname);

	int ihcell;

	for (i=0;i<nhcells;i++) hcweights[i] = 0;

	for (i=0;i<_N1LevelCellsNew;i++) {
		ihcell = _1LevelCellsNew2HCells[i];
		hcweights[ihcell] += ia1Level2RealNew[i+1]-ia1Level2RealNew[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nhcells, hcweights, hcweights);

	int *hco2n;

	hco2n = new int [nhcells+1];
	if (!hco2n) MemoryFail (funcname);

	for (i=0;i<=nhcells;i++) hco2n[i] = 0;

	int j;

	int ibeg, iend, icell, icell1, icount, icount1, nhcellsloc;

	int index = 0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		nhcellsloc = 0;
		icell = ibeg;
		while (icell <= iend) {
			icount = ia1Level2RealNew[icell+1]-ia1Level2RealNew[icell];
			if (icount < weightcellvolumeloc) {
				while (icell+1 <= iend) {
					icell1 = icell+1;
					icount1 = icount + ia1Level2RealNew[icell1+1]-ia1Level2RealNew[icell1];
					if (icount1 <= weightcellvolumeloc) {
						icount = icount1;
						icell = icell1;
					} else {
						break;
					};
				};
			};
			nhcellsloc++;
			icell++;
		};
		if (nhcellsloc == 0) nhcellsloc = 1;
		hco2n[ihcell+1] = nhcellsloc;
		index = iend+1;
	};

	for (i=0;i<nhcells;i++) {
		if (hcell2cpu[i] == myid && hco2n[i+1] == 0) hco2n[i+1] = 1;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nhcells+1, hco2n, hco2n);

	for (i=0;i<nhcells;i++) hco2n[i+1] = hco2n[i]+hco2n[i+1];

	_NHyperCells = hco2n[nhcells];

	_HyperCell2CPU = new int [_NHyperCells];
	if(!_HyperCell2CPU) MemoryFail(funcname);

	for (i=0;i<nhcells;i++) {
		for (j=hco2n[i];j<hco2n[i+1];j++) {
			_HyperCell2CPU[j] = hcell2cpu[i];
		};
	};

	_HyperCells = new int [_NHyperCells+1];
	if(!_HyperCells) MemoryFail(funcname);

	_HyperCells[0] = 0;

	for (i=0;i<=_NHyperCells;i++) _HyperCells[i] = 0;

	index = 0;

	int ibegnew, iendnew, ihcellnew, ibeg0, iend0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		icell = ibeg;
		ibeg0 = ibeg;
		iend0 = ibeg-1;
		ibegnew = icell;
		iendnew = icell;
		icount = ia1Level2RealNew[icell+1]-ia1Level2RealNew[icell];
		ihcellnew = hco2n[ihcell];
		while (ihcellnew < hco2n[ihcell+1]-1) {
			while (icell+1 <= iend) {
				icell1 = icell+1;
				icount1 = icount + ia1Level2RealNew[icell1+1]-ia1Level2RealNew[icell1];
				if (icount1 <= weightcellvolumeloc) {
					icount = icount1;
					icell = icell1;
					iendnew = icell;
				} else {
					iendnew = icell;
					break;
				};
			};
			_HyperCells[ihcellnew+1] = iendnew-ibegnew+1;
			ihcellnew++;
			ibegnew = icell+1;
			iendnew = icell+1;
			iend0 = icell;
			icell = ibegnew;
			icount = ia1Level2RealNew[icell+1]-ia1Level2RealNew[icell];
		};
		_HyperCells[hco2n[ihcell+1]] = (hcells[ihcell+1]-hcells[ihcell])-(iend0-ibeg0+1);
		index = iend+1;
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells+1, _HyperCells, _HyperCells);

	for (i=0;i<_NHyperCells;i++) _HyperCells[i+1] = _HyperCells[i]+_HyperCells[i+1];

// Modify array 1LevelCellsNew2HCells

	int ibegn;

	index = 0;

	while (index < _N1LevelCellsNew) {
		ihcell = _1LevelCellsNew2HCells[index];
		ibeg = index;
		while (index < _N1LevelCellsNew-1) {
			if (_1LevelCellsNew2HCells[index+1] == ihcell) {
				index++;
			} else {
				break;
			};
		};
		iend = index;
		if (iend-ibeg+1 != hcells[ihcell+1]-hcells[ihcell]) throw " CDecomp::DistributeCellsOverCPUsNestedDessectionParMetisWeights: error in array 1LevelCellsNew2HCells";
		ibegn = hcells[ihcell];
		for (ihcellnew=hco2n[ihcell];ihcellnew<hco2n[ihcell+1];ihcellnew++) {
			for (i=_HyperCells[ihcellnew];i<_HyperCells[ihcellnew+1];i++) {
				_1LevelCellsNew2HCells[i-ibegn+ibeg] = ihcellnew;
			};
		};
		index = iend+1;
	};

// Count weights of the new hypercells

	int *hcweightsnew;

	hcweightsnew = new int [_NHyperCells];
	if (!hcweightsnew) MemoryFail (funcname);

	for (i=0;i<_NHyperCells;i++) hcweightsnew[i] = 0;

	for (i=0;i<_N1LevelCellsNew;i++) {
		ihcell = _1LevelCellsNew2HCells[i];
		hcweightsnew[ihcell] += ia1Level2RealNew[i+1]-ia1Level2RealNew[i];
	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _NHyperCells, hcweightsnew, hcweightsnew);

// Register new hypercells in the tree's

	CTree *ptree = _CpuTree.GetTree ();

	ptree->ReassignBlockNumbers (hco2n);

	ptree->StorePartitioning (_NHyperCells, _HyperCells, _HyperCell2CPU);

// Free work arrays

	delete [] hcell2cpu;
	delete [] hcells;
	delete [] ia1Level2RealNew;
	delete [] nRealCellsNewIn1Level;
	delete [] hcweights;
	delete [] hco2n;
	delete [] hcweightsnew;

	return 0;

};

// Author: Kharchenko S.A.
// CDecompHCells: Memory allocation zero data constructor
//========================================================================================
CDecompHCells::CDecompHCells () { // Memory allocation zero data constructor

	const char *funcname = "CDecompHCells_00";

	nhcells = 0;
	hcells = new int [nhcells+1];
	if (!hcells) MemoryFail (funcname);
	hcells[0] = 0;
	cellschars = new char ** [nhcells];
	if (!cellschars) MemoryFail (funcname);
	sizeofbox = 0;
	cellsboxes = new char * [nhcells];
	if (!cellsboxes) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CDecompHCells: Memory allocation zero data constructor
//========================================================================================
CDecompHCells::CDecompHCells (int _nhcells, int *_hcells) { // Memory allocation zero data constructor

	const char *funcname = "CDecompHCells_01";

	nhcells = _nhcells;
	hcells = new int [nhcells+1];
	if (!hcells) MemoryFail (funcname);

	int i, j, ni;

	for (i=0;i<nhcells+1;i++) {
		hcells[i] = _hcells[i];
	};

	cellschars = new char ** [nhcells];
	if (!cellschars) MemoryFail (funcname);

	const char *chdummy = "";

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		cellschars[i] = new char * [ni];
		if (!cellschars[i]) MemoryFail (funcname);
		for (j=0;j<ni;j++) {
			cellschars[i][j] = new char [1];
			if (!cellschars[i][j]) MemoryFail (funcname);
			strcpy (cellschars[i][j],chdummy);
		};
	};

	sizeofbox = 0;
	cellsboxes = new char * [nhcells];
	if (!cellsboxes) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CDecompHCells: Equality operator
//========================================================================================
CDecompHCells &CDecompHCells::operator= (const CDecompHCells &_obj) { // Equality operator

	const char *funcname = "CNode_=";

	int i, j, ni;

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		for (j=0;j<ni;j++) {
			delete [] cellschars[i][j];
		};
		delete [] cellschars[i];
		delete [] cellsboxes[i];
	};
	delete [] cellschars;
	delete [] cellsboxes;

	delete [] hcells;

	nhcells = _obj.nhcells;
	hcells = new int [nhcells+1];
	if (!hcells) MemoryFail (funcname);

	for (i=0;i<nhcells+1;i++) {
		hcells[i] = _obj.hcells[i];
	};

	cellschars = new char ** [nhcells];
	if (!cellschars) MemoryFail (funcname);

	size_t isize;

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		cellschars[i] = new char * [ni];
		if (!cellschars[i]) MemoryFail (funcname);
		for (j=0;j<ni;j++) {
			isize = strlen(_obj.cellschars[i][j])+1;
			cellschars[i][j] = new char [isize];
			if (!cellschars[i][j]) MemoryFail (funcname);
			strcpy (cellschars[i][j],_obj.cellschars[i][j]);
		};
	};

	sizeofbox = _obj.sizeofbox;

	cellsboxes = new char * [nhcells];
	if (!cellsboxes) MemoryFail (funcname);

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		cellsboxes[i] = new char [ni*sizeofbox];
		if (!cellsboxes[i]) MemoryFail (funcname);
		memcpy (cellsboxes[i],_obj.cellsboxes[i],ni*sizeofbox);
	};

	return *this;

};

// Author: Kharchenko S.A.
// CDecompHCells: Fill the data that describe the boxes
//========================================================================================
void CDecompHCells::FillBoxesData (CMPIComm &_comm, int _nblks, int *_blks, // Fill the data that describe the boxes
												int _nlist, int *_list, void *_boxes, int _sizeofbox) {

	const char *funcname = "FillBoxesData";

// Allocate the memory

	int ntot = _blks[_nblks];

	char *boxestot;

	boxestot = new char [ntot*_sizeofbox];
	if (!boxestot) MemoryFail (funcname);

	int i, j;

	for (i=0;i<ntot*_sizeofbox;i++) boxestot[i] = 0;

	for (i=0;i<_nlist;i++) {
		j = _list[i];
		memcpy (boxestot+j*_sizeofbox, ((char *)_boxes)+i*_sizeofbox, _sizeofbox);
	};

	CMPIExchange::ExchangeArrayMPI (_comm, CHARVALUE, ADD, ntot*_sizeofbox, boxestot, boxestot);

// Create new structure

	CDecompHCells hc (_nblks, _blks);

	hc.SetSizeofbox (_sizeofbox);

	char *ptr;

	int ni, ibeg;

	for (i=0;i<_nblks;i++) {
		ni = _blks[i+1]-_blks[i];
		ibeg = _blks[i];
		ptr = new char [ntot*_sizeofbox];
		if (!ptr) MemoryFail (funcname);
		memcpy (ptr, boxestot+ibeg*_sizeofbox, ni*_sizeofbox);
		hc.SetCellsboxes (i, ptr);
	};

// Free work arrays

	delete [] boxestot;

	*this = hc;

};

// Author: Kharchenko S.A.
// CDecompHCells: Destructor
//========================================================================================
CDecompHCells::~CDecompHCells () { // Destructor

	int i, j, ni;

	for (i=0;i<nhcells;i++) {
		ni = hcells[i+1]-hcells[i];
		for (j=0;j<ni;j++) {
			delete [] cellschars[i][j];
		};
		delete [] cellschars[i];
		delete [] cellsboxes[i];
	};
	delete [] cellschars;
	delete [] cellsboxes;

	nhcells = 0;
	delete [] hcells;

};

// Author: Kharchenko S.A.
// CDecompHCells: Output hcells
//========================================================================================
ostream &operator<< (ostream &_stream, const CDecompHCells &_obj) { // Output hcells

	_stream << " CDecompHCells:" << endl;

	_stream << " Nhcells = " << _obj.nhcells << endl;
	OutArr (_stream, " HCells = ",_obj.nhcells+1,_obj.hcells);

	int i, j, ni;

	for (i=0;i<_obj.nhcells;i++) {
		ni = _obj.hcells[i+1]-_obj.hcells[i];
		_stream << "    IHCell = " << i << " Isize = " << ni << endl;
		for (j=0;j<ni;j++) {
			_stream << "        JCell = " << j << " Val = <" << _obj.cellschars[i][j] << ">" << endl;
		};
	};

	return _stream;

};
