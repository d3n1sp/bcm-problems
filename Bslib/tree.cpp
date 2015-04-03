//------------------------------------------------------------------------------------------------
// File: tree.cpp
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
#include <cstdlib>

#include "globals.h"
#include "smatrix.h"
#include "tree.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif

// Author: Kharchenko S.A.
// CNode: Memory allocation zero data constructor
//========================================================================================
CNode::CNode () { // Memory allocation zero data constructor

	const char *funcname = "CNode_00";

	nodetype = undefnode;
	nodeid = -1;
	nodeidext = -1;
	nodeidsm = -1;
	nodecpu = -1;
	nodelv = -1;
	nodeperf = 1.0e0;
	indbeg = -1;
	indend = -1;
	indbegtot = -1;
	indendtot = -1;
	fatherid = -1;
	fathercpu = -1;
	fatherid2 = -1;
	fathercpu2 = -1;
	nchilds = 0;
	nchilds2 = 0;
	childs = new int [nchilds+1];
	if (!childs) MemoryFail (funcname);
	childs2 = new int [nchilds2+1];
	if (!childs2) MemoryFail (funcname);
	int i;
	for (i=0;i<nchilds+1;i++) {
		childs[i] = -1;
	};
	for (i=0;i<nchilds2+1;i++) {
		childs2[i] = -1;
	};
	indtree = -1;
	nprocgl = 0;
	nprocnode = 0;
	cpuidarrnode = new int [nprocnode];
	if (!cpuidarrnode) MemoryFail (funcname);
	cpuidarrnodeloc = new int [nprocnode];
	if (!cpuidarrnodeloc) MemoryFail (funcname);
	cpuarrndgl2loc = new int [nprocgl];
	if (!cpuarrndgl2loc) MemoryFail (funcname);
};

// Author: Kharchenko S.A.
// CNode: Memory allocation zero data constructor
//========================================================================================
CNode::CNode (int _nchilds) { // Memory allocation zero data constructor

	const char *funcname = "CNode_01";

	nodetype = undefnode;
	nodeid = -1;
	nodeidext = -1;
	nodeidsm = -1;
	nodecpu = -1;
	nodelv = -1;
	nodeperf = 1.0e0;
	indbeg = -1;
	indend = -1;
	indbegtot = -1;
	indendtot = -1;
	fatherid = -1;
	fathercpu = -1;
	fatherid2 = -1;
	fathercpu2 = -1;
	nchilds = _nchilds;
	nchilds2 = 0;
	childs = new int [nchilds+1];
	if (!childs) MemoryFail (funcname);
	childs2 = new int [nchilds2+1];
	if (!childs2) MemoryFail (funcname);
	int i;
	for (i=0;i<nchilds+1;i++) {
		childs[i] = -1;
	};
	for (i=0;i<nchilds2+1;i++) {
		childs2[i] = -1;
	};
	indtree = -1;
	nprocgl = 0;
	nprocnode = 0;
	cpuidarrnode = new int [nprocnode];
	if (!cpuidarrnode) MemoryFail (funcname);
	cpuidarrnodeloc = new int [nprocnode];
	if (!cpuidarrnodeloc) MemoryFail (funcname);
	cpuarrndgl2loc = new int [nprocgl];
	if (!cpuarrndgl2loc) MemoryFail (funcname);
};

// Author: Kharchenko S.A.
// CNode: Memory allocation zero data constructor
//========================================================================================
CNode::CNode (int _nchilds, int _nchilds2) { // Memory allocation zero data constructor

	const char *funcname = "CNode_02";

	nodetype = undefnode;
	nodeid = -1;
	nodeidext = -1;
	nodeidsm = -1;
	nodecpu = -1;
	nodelv = -1;
	nodeperf = 1.0e0;
	indbeg = -1;
	indend = -1;
	indbegtot = -1;
	indendtot = -1;
	fatherid = -1;
	fathercpu = -1;
	fatherid2 = -1;
	fathercpu2 = -1;
	nchilds = _nchilds;
	nchilds2 = _nchilds2;
	childs = new int [nchilds+1];
	if (!childs) MemoryFail (funcname);
	childs2 = new int [nchilds2+1];
	if (!childs2) MemoryFail (funcname);
	int i;
	for (i=0;i<nchilds+1;i++) {
		childs[i] = -1;
	};
	for (i=0;i<nchilds2+1;i++) {
		childs2[i] = -1;
	};
	indtree = -1;
	nprocgl = 0;
	nprocnode = 0;
	cpuidarrnode = new int [nprocnode];
	if (!cpuidarrnode) MemoryFail (funcname);
	cpuidarrnodeloc = new int [nprocnode];
	if (!cpuidarrnodeloc) MemoryFail (funcname);
	cpuarrndgl2loc = new int [nprocgl];
	if (!cpuarrndgl2loc) MemoryFail (funcname);
};


// Author: Kharchenko S.A.
// CNode: Copy constructor
//========================================================================================
CNode::CNode (const CNode &_node) { // Copy constructor

	const char *funcname = "CNode_copy";

	nodetype = _node.nodetype;
	nodeid = _node.nodeid;
	nodeidext = _node.nodeidext;
	nodeidsm = _node.nodeidsm;
	nodecpu = _node.nodecpu;
	nodelv = _node.nodelv;
	nodeperf = _node.nodeperf;
	indbeg = _node.indbeg;
	indend = _node.indend;
	indbegtot = _node.indbegtot;
	indendtot = _node.indendtot;
	fatherid = _node.fatherid;
	fathercpu = _node.fathercpu;
	fatherid2 = _node.fatherid2;
	fathercpu2 = _node.fathercpu2;
	nchilds = _node.nchilds;
	nchilds2 = _node.nchilds2;
	childs = new int [nchilds+1];
	if (!childs) MemoryFail (funcname);
	int i;
	for (i=0;i<nchilds;i++) {
		childs[i] = _node.childs[i];
	};
	childs2 = new int [nchilds2+1];
	if (!childs2) MemoryFail (funcname);
	for (i=0;i<nchilds2;i++) {
		childs2[i] = _node.childs2[i];
	};
	indtree = _node.indtree;
	nprocgl = _node.nprocgl;
	nprocnode = _node.nprocnode;
	cpuidarrnode = new int [nprocnode];
	if (!cpuidarrnode) MemoryFail (funcname);
	for (i=0;i<nprocnode;i++) {
		cpuidarrnode[i] = _node.cpuidarrnode[i];
	};
	cpuidarrnodeloc = new int [nprocnode];
	if (!cpuidarrnodeloc) MemoryFail (funcname);
	for (i=0;i<nprocnode;i++) {
		cpuidarrnodeloc[i] = _node.cpuidarrnodeloc[i];
	};
	cpuarrndgl2loc = new int [nprocgl];
	if (!cpuarrndgl2loc) MemoryFail (funcname);
	for (i=0;i<nprocgl;i++) {
		cpuarrndgl2loc[i] = _node.cpuarrndgl2loc[i];
	};
};

// Author: Kharchenko S.A.
// CNode: Equality operator
//========================================================================================
CNode &CNode::operator= (const CNode &_node2) { // Equality operator

	const char *funcname = "CNode_=";

	delete [] childs;
	delete [] childs2;
	delete [] cpuidarrnode;
	delete [] cpuidarrnodeloc;
	delete [] cpuarrndgl2loc;

	nodetype = _node2.nodetype;
	nodeid = _node2.nodeid;
	nodeidext = _node2.nodeidext;
	nodeidsm = _node2.nodeidsm;
	nodecpu = _node2.nodecpu;
	nodelv = _node2.nodelv;
	nodeperf = _node2.nodeperf;
	indbeg = _node2.indbeg;
	indend = _node2.indend;
	indbegtot = _node2.indbegtot;
	indendtot = _node2.indendtot;
	fatherid = _node2.fatherid;
	fathercpu = _node2.fathercpu;
	fatherid2 = _node2.fatherid2;
	fathercpu2 = _node2.fathercpu2;
	nchilds = _node2.nchilds;
	nchilds2 = _node2.nchilds2;
	childs = new int [nchilds+1];
	if (!childs) MemoryFail (funcname);
	int i;
	for (i=0;i<nchilds;i++) {
		childs[i] = _node2.childs[i];
	};
	childs2 = new int [nchilds2+1];
	if (!childs2) MemoryFail (funcname);
	for (i=0;i<nchilds2;i++) {
		childs2[i] = _node2.childs2[i];
	};
	indtree = _node2.indtree;
	nprocgl = _node2.nprocgl;
	nprocnode = _node2.nprocnode;
	cpuidarrnode = new int [nprocnode];
	if (!cpuidarrnode) MemoryFail (funcname);
	for (i=0;i<nprocnode;i++) {
		cpuidarrnode[i] = _node2.cpuidarrnode[i];
	};
	cpuidarrnodeloc = new int [nprocnode];
	if (!cpuidarrnodeloc) MemoryFail (funcname);
	for (i=0;i<nprocnode;i++) {
		cpuidarrnodeloc[i] = _node2.cpuidarrnodeloc[i];
	};
	cpuarrndgl2loc = new int [nprocgl];
	if (!cpuarrndgl2loc) MemoryFail (funcname);
	for (i=0;i<nprocgl;i++) {
		cpuarrndgl2loc[i] = _node2.cpuarrndgl2loc[i];
	};

	return *this;

};

// Author: Kharchenko S.A.
// CNode: Destructor
//========================================================================================
CNode::~CNode () { // Destructor
//	cout << " On entry to CNode destructor " << endl;
	nodetype = undefnode;
	nodeid = -1;
	nodeidext = -1;
	nodeidsm = -1;
	nodecpu = -1;
	nodelv = -1;
	nodeperf = 0.0e0;
	indbeg = -1;
	indend = -1;
	indbegtot = -1;
	indendtot = -1;
	fatherid = -1;
	fathercpu = -1;
	fatherid2 = -1;
	fathercpu2 = -1;
	nchilds = 0;
	nchilds2 = 0;
	delete [] childs;
	delete [] childs2;
	nprocgl = 0;
	nprocnode = 0;
	delete [] cpuidarrnode;
	delete [] cpuidarrnodeloc;
	delete [] cpuarrndgl2loc;
//	cout << " On return from CNode destructor " << endl;
};

// Author: Kharchenko S.A.
// CNode: Output node
//========================================================================================
ostream &operator<< (ostream &_stream, const CNode &_node) { // Output node

	_stream << " CNode:" << endl;

	_stream << " NodeType = " << _node.nodetype << endl;
	_stream << " NodeId = " << _node.nodeid << " NodeIdExt = " << _node.nodeidext << " NodeIdSm = " << _node.nodeidsm << endl;
	_stream << " NodeCpu = " << _node.nodecpu << " NodeLv = " << _node.nodelv << endl;
	_stream << " NodePerf = " << _node.nodeperf << endl;
	_stream << " IndBeg = " << _node.indbeg << " IndEnd = " << _node.indend << endl;
	_stream << " IndBegTot = " << _node.indbegtot << " IndEndTot = " << _node.indendtot << endl;
	_stream << " FatherId = " << _node.fatherid << " FatherCpu = " << _node.fathercpu << endl;
	_stream << " FatherId2 = " << _node.fatherid2 << " FatherCpu2 = " << _node.fathercpu2 << endl;
	_stream << " Nchilds = " << _node.nchilds << endl;
	_stream << " Nchilds2 = " << _node.nchilds2 << endl;

	OutArr (_stream, " Childs =", _node.nchilds, _node.childs);
	OutArr (_stream, " Childs2 =", _node.nchilds2, _node.childs2);

	_stream << " NprocGl = " << _node.nprocgl << endl;
	_stream << " NprocNode = " << _node.nprocnode << endl;
	OutArr (_stream, " CpuIdArrNode =", _node.nprocnode, _node.cpuidarrnode);
	OutArr (_stream, " CpuIdArrNodeLoc =", _node.nprocnode, _node.cpuidarrnodeloc);
	OutArr (_stream, " CpuArrNdGl2Loc =", _node.nprocgl, _node.cpuarrndgl2loc);
//	_stream << " CommNode = " << _node.commnode << endl;

	return _stream;

};

// Author: Kharchenko S.A.
// CNode: Order node indices inside the node
//========================================================================================
void CNode::OrderIndices (int *_order) { // Order node indices inside the node

	int i, ind;

	if (nodeid >= 0) {
		ind = _order[nodeid];
		nodeid = ind;
	};
	if (fatherid >= 0) {
		ind = _order[fatherid];
		fatherid = ind;
	};
	if (fatherid2 >= 0) {
		ind = _order[fatherid2];
		fatherid2 = ind;
	};
	for (i=0;i<nchilds;i++) {
		if (childs[i] >= 0) {
			ind = _order[childs[i]];
			childs[i] = ind;
		};
	};
	for (i=0;i<nchilds2;i++) {
		if (childs2[i] >= 0) {
			ind = _order[childs2[i]];
			childs2[i] = ind;
		};
	};

};

// Author: Kharchenko S.A.
// CNode: Pack node data
//========================================================================================
void CNode::PackNode (int &_length, char *&_obj) { // Pack node data

	const char* funcname = "PackNode";

// Allocate the structure with the head and init pointers

	_length = 19 * sizeof(int) + sizeof(double) + (nchilds+nchilds2+2*nprocnode+nprocgl) * sizeof(int);

	_obj = new char [_length];
	if (!_obj) MemoryFail(funcname);

	char* pLoc;

	pLoc = _obj;

	int* pHead;
	double* pdouble;
	int* pchilds;
	int* pchilds2;
	int* pcpuidarrnode;
	int* pcpuidarrnodeloc;
	int* pcpuarrndgl2loc;

	pHead = (int *) pLoc;
	pLoc += 19 * sizeof(int);

	pdouble = (double *) pLoc;
	pLoc += sizeof(double);

	pchilds = (int *) pLoc;
	pLoc += nchilds * sizeof(int);

	pchilds2 = (int *) pLoc;
	pLoc += nchilds2 * sizeof(int);

	pcpuidarrnode = (int *) pLoc;
	pLoc += nprocnode * sizeof(int);

	pcpuidarrnodeloc = (int *) pLoc;
	pLoc += nprocnode * sizeof(int);

	pcpuarrndgl2loc = (int *) pLoc;
	pLoc += nprocgl * sizeof(int);

// Store the data

	pHead[0]  = nodetype;
	pHead[1]  = nodeid;
	pHead[2]  = nodeidext;
	pHead[3]  = nodeidsm;
	pHead[4]  = nodecpu;
	pHead[5]  = nodelv;
	pHead[6]  = indbeg;
	pHead[7]  = indend;
	pHead[8]  = indbegtot;
	pHead[9]  = indendtot;
	pHead[10] = fatherid;
	pHead[11] = fathercpu;
	pHead[12] = fatherid2;
	pHead[13] = fathercpu2;
	pHead[14] = nchilds;
	pHead[15] = nchilds2;
	pHead[16] = indtree;
	pHead[17] = nprocgl;
	pHead[18] = nprocnode;

	pdouble[0] = nodeperf;

	int j;

	for(j = 0; j < nchilds; j++) pchilds[j] = childs[j];
	for(j = 0; j < nchilds2; j++) pchilds2[j] = childs2[j];
	for(j = 0; j < nprocnode; j++) pcpuidarrnode[j] = cpuidarrnode[j];
	for(j = 0; j < nprocnode; j++) pcpuidarrnodeloc[j] = cpuidarrnodeloc[j];
	for(j = 0; j < nprocgl; j++) pcpuarrndgl2loc[j] = cpuarrndgl2loc[j];

};

// Author: Kharchenko S.A.
// CNode: UnPack node data
//========================================================================================
void CNode::UnPackNode (int _length, char *_obj) { // UnPack node data

	const char* funcname = "UnPackNode";

// Delete previous arrays

	delete [] childs;
	delete [] childs2;
	delete [] cpuidarrnode;
	delete [] cpuidarrnodeloc;
	delete [] cpuarrndgl2loc;

// Open the head

	char* pLoc;

	pLoc = _obj;

	int* pHead;
	double* pdouble;
	int* pchilds;
	int* pchilds2;
	int* pcpuidarrnode;
	int* pcpuidarrnodeloc;
	int* pcpuarrndgl2loc;

	pHead = (int *) pLoc;
	pLoc += 19 * sizeof(int);

	pdouble = (double *) pLoc;
	pLoc += sizeof(double);

// Get head data

	nodetype   = pHead[0];
	nodeid     = pHead[1];
	nodeidext  = pHead[2];
	nodeidsm   = pHead[3];
	nodecpu    = pHead[4];
	nodelv     = pHead[5];
	indbeg     = pHead[6];
	indend     = pHead[7];
	indbegtot  = pHead[8];
	indendtot  = pHead[9];
	fatherid   = pHead[10];
	fathercpu  = pHead[11];
	fatherid2  = pHead[12];
	fathercpu2 = pHead[13];
	nchilds    = pHead[14];
	nchilds2   = pHead[15];
	indtree    = pHead[16];
	nprocgl    = pHead[17];
	nprocnode  = pHead[18];

	nodeperf = pdouble[0];

// Get childs data

	childs = new int [nchilds];
	if (!childs) MemoryFail(funcname);
	childs2 = new int [nchilds2];
	if (!childs2) MemoryFail(funcname);
	cpuidarrnode = new int [nprocnode];
	if (!cpuidarrnode) MemoryFail(funcname);
	cpuidarrnodeloc = new int [nprocnode];
	if (!cpuidarrnodeloc) MemoryFail(funcname);
	cpuarrndgl2loc = new int [nprocgl];
	if (!cpuarrndgl2loc) MemoryFail(funcname);

	pchilds = (int *) pLoc;
	pLoc += nchilds * sizeof(int);

	pchilds2 = (int *) pLoc;
	pLoc += nchilds2 * sizeof(int);

	pcpuidarrnode = (int *) pLoc;
	pLoc += nprocnode * sizeof(int);

	pcpuidarrnodeloc = (int *) pLoc;
	pLoc += nprocnode * sizeof(int);

	pcpuarrndgl2loc = (int *) pLoc;
	pLoc += nprocgl * sizeof(int);

	int j;

	for(j = 0; j < nchilds;  j++) childs[j] = pchilds[j];
	for(j = 0; j < nchilds2; j++) childs2[j] = pchilds2[j];
	for(j = 0; j < nprocnode;j++) cpuidarrnode[j] = pcpuidarrnode[j];
	for(j = 0; j < nprocnode;j++) cpuidarrnodeloc[j] = pcpuidarrnodeloc[j];
	for(j = 0; j < nprocgl;  j++) cpuarrndgl2loc[j] = pcpuarrndgl2loc[j];

};

// Author: Kharchenko S.A.
// CTree: Memory allocation zero data constructor
//========================================================================================
CTree::CTree () { // Memory allocation zero data constructor

	const char *funcname = "CTree_00";

	myid = 0;
	nlev = -1;
	nnodes = 0;
	nproc = 0;
	nproctot = 0;
	rootid = -1;
	rootcpu = -1;
	nodes = new CNode [nnodes+1];
	if (!nodes) MemoryFail (funcname);
	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	int i;
	for (i=0;i<nproctot;i++) {
		cpuidarr[i] = -1;
	};
	cpuidbeg = new int [nproc+1];
	if (!cpuidbeg) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidbeg[i] = -1;
	};
	cpuidend = new int [nproc+1];
	if (!cpuidend) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidend[i] = -1;
	};
	perf2cpu = new double [nproc+1];
	if (!perf2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		perf2cpu[i] = 1.0e0;
	};
	memory2cpu = new double [nproc+1];
	if (!memory2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		memory2cpu[i] = 1.0e0;
	};
	nsubtree = 0;
	subtreearr = 0;
	nblkstree = 0;
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	blkstree[0] = 0;

};

// Author: Kharchenko S.A.
// CTree: Memory allocation zero data constructor
//========================================================================================
CTree::CTree (int _nproc) { // Memory allocation zero data constructor

	const char *funcname = "CTree_01";

	myid = 0;
	nlev = -1;
	nnodes = 12*_nproc+10;
	nproc = _nproc;
	nproctot = _nproc;
	rootid = -1;
	rootcpu = -1;
	nodes = new CNode [nnodes+1];
	if (!nodes) MemoryFail (funcname);
	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	int i;
	for (i=0;i<nproctot;i++) {
		cpuidarr[i] = -1;
	};
	cpuidbeg = new int [nproc+1];
	if (!cpuidbeg) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidbeg[i] = -1;
	};
	cpuidend = new int [nproc+1];
	if (!cpuidend) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidend[i] = -1;
	};
	perf2cpu = new double [nproc+1];
	if (!perf2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		perf2cpu[i] = 1.0e0;
	};
	memory2cpu = new double [nproc+1];
	for (i=0;i<nproc;i++) {
		memory2cpu[i] = 1.0e0;
	};
	nsubtree = 0;
	subtreearr = 0;
	nblkstree = 0;
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	blkstree[0] = 0;

};

// Author: Kharchenko S.A.
// CTree: Create the tree
//========================================================================================
CTree::CTree (int _nproc, int _nchilds) { // Create the tree

	const char *funcname = "CTree_02";

	myid = 0;
	nlev = -1;
	nnodes = 4*_nproc+10;
	nproc = _nproc;
	nproctot = _nproc;
	rootid = -1;
	rootcpu = -1;
	nodes = new CNode [nnodes+1];
	if (!nodes) MemoryFail (funcname);
	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	int i;
	for (i=0;i<nproctot;i++) {
		cpuidarr[i] = -1;
	};
	cpuidbeg = new int [nproc+1];
	if (!cpuidbeg) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidbeg[i] = -1;
	};
	cpuidend = new int [nproc+1];
	if (!cpuidend) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidend[i] = -1;
	};
	perf2cpu = new double [nproc+1];
	if (!perf2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		perf2cpu[i] = 1.0e0;
	};
	memory2cpu = new double [nproc+1];
	for (i=0;i<nproc;i++) {
		memory2cpu[i] = 1.0e0;
	};
	nsubtree = 0;
	subtreearr = 0;
	nblkstree = 0;
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	blkstree[0] = 0;

// Create 1 level tree

	if (_nproc == 1) {

		nnodes = 1;
		nlev = 1;
		rootid = 0;
		rootcpu = 0;

		CNode temp (1);

		temp.nodetype = compnode;
		temp.nodeid = 0;
		temp.nodecpu = 0;
		temp.nodelv = 0;
		temp.fatherid = 0;
		temp.fathercpu = 0;
		temp.fatherid2 = 0;
		temp.fathercpu2 = 0;
		temp.nchilds = 1;
		temp.childs[0] = 0;

		nodes[0] = temp;

	} else {

// Create 2 level tree

		if (_nchilds == 1) {

			nnodes = nproc+1;
			nlev = 2;
			rootid = 0;
			rootcpu = 0;

// Root node

			CNode temp (nproc);

			temp.nodetype = compnode;
			temp.nodeid = 0;
			temp.nodecpu = 0;
			temp.nodelv = 0;
			temp.fatherid = 0;
			temp.fathercpu = 0;
			temp.fatherid2 = 0;
			temp.fathercpu2 = 0;
			temp.nchilds = nproc;
			for (int i=0;i<nproc;i++) {
				temp.childs[i] = i+1;
			};

			nodes[0] = temp;

// Child nodes

			for (int iproc=0;iproc<nproc;iproc++) {

				CNode temp2 (1);

				temp2.nodetype = compnode;
				temp2.nodeid = iproc+1;
				temp2.nodecpu = iproc;
				temp2.nodelv = 1;
				temp2.fatherid = 0;
				temp2.fathercpu = 0;
				temp2.fatherid2 = 0;
				temp2.fathercpu2 = 0;
				temp2.nchilds = 1;
				temp2.childs[0] = iproc+1;

				nodes[iproc+1] = temp2;

			};

		} else {

			if (nproc%_nchilds != 0) {
				throw " Error in the number of processes used ";
			};

			int nphead = nproc / _nchilds;

			if (true) {

// Create 3 level tree

				nnodes = nproc+nphead+1;
				nlev = 3;
				rootid = 0;
				rootcpu = 0;

// Root node

				CNode temp (nphead);

				temp.nodetype = compnode;
				temp.nodeid = 0;
				temp.nodecpu = 0;
				temp.nodelv = 0;
				temp.fatherid = 0;
				temp.fathercpu = 0;
				temp.fatherid2 = 0;
				temp.fathercpu2 = 0;
				temp.nchilds = nphead;
				for (int i=0;i<nphead;i++) {
					temp.childs[i] = i+1;
				};

				nodes[0] = temp;

// Child nodes

				int iproc;
				for (iproc=0;iproc<nphead;iproc++) {

					CNode temp2 (_nchilds);

					temp2.nodetype = compnode;
					temp2.nodeid = iproc+1;
					temp2.nodecpu = iproc*_nchilds;
					temp2.nodelv = 1;
					temp2.fatherid = 0;
					temp2.fathercpu = 0;
					temp2.fatherid2 = 0;
					temp2.fathercpu2 = 0;
					temp2.nchilds = _nchilds;

					for (int j=0;j<_nchilds;j++) {
						temp2.childs[j] = 1+nphead+iproc*_nchilds+j;
					};

					nodes[iproc+1] = temp2;

				};

// Child child nodes

				for (iproc=0;iproc<nphead;iproc++) {

					for (int jproc=0;jproc<_nchilds;jproc++) {

						CNode temp3 (1);

						temp3.nodetype = compnode;
						temp3.nodeid = 1+nphead+iproc*_nchilds+jproc;
						temp3.nodecpu = iproc*_nchilds+jproc;
						temp3.nodelv = 2;
						temp3.fatherid = iproc+1;
						temp3.fathercpu = iproc*_nchilds;
						temp3.fatherid2 = iproc+1;
						temp3.fathercpu2 = iproc*_nchilds;
						temp3.nchilds = 1;

						temp3.childs[0] = 1+nphead+iproc*_nchilds+jproc;

						nodes[1+nphead+iproc*_nchilds+jproc] = temp3;

					};

				};

			};
		};

	};

// Assign cpuidbeg and cpuidend

	int iproc;
	for (iproc=0;iproc<nproc;iproc++) {
		cpuidbeg[iproc] = -1;
		cpuidend[iproc] = -1;
	};
	int inode;
	for (inode=nnodes-1;inode>=0;inode--) {
		iproc = nodes[inode].nodecpu;
		cpuidbeg[iproc] = inode;
	};
	for (inode=0;inode<nnodes;inode++) {
		iproc = nodes[inode].nodecpu;
		cpuidend[iproc] = inode;
	};

};

// Author: Kharchenko S.A.
// CTree: Create the tree
//========================================================================================
CTree::CTree (int _nproc, int _nchilds, double *_perf2cpu, double *_memory2cpu) { // Create the tree

	const char *funcname = "CTree_03";

	myid = 0;
	nlev = -1;
	nnodes = 4*_nproc+10;
	nproc = _nproc;
	nproctot = _nproc;
	rootid = -1;
	rootcpu = -1;
	nodes = new CNode [nnodes+1];
	if (!nodes) MemoryFail (funcname);
	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	int i;
	for (i=0;i<nproctot;i++) {
		cpuidarr[i] = -1;
	};
	cpuidbeg = new int [nproc+1];
	if (!cpuidbeg) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidbeg[i] = -1;
	};
	cpuidend = new int [nproc+1];
	if (!cpuidend) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidend[i] = -1;
	};
	perf2cpu = new double [nproc+1];
	if (!perf2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		perf2cpu[i] = _perf2cpu[i];
	};
	memory2cpu = new double [nproc+1];
	for (i=0;i<nproc;i++) {
		memory2cpu[i] = _memory2cpu[i];
	};
	nsubtree = 0;
	subtreearr = 0;
	nblkstree = 0;
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	blkstree[0] = 0;

// Create temporary node

	int nchildsloc = _nchilds+1;

	CNode temp (nchildsloc);

// Init the lowest level nodes

	int *sprndsini;

	sprndsini = new int [nproc+2];
	if (!sprndsini) MemoryFail (funcname);

	nlev = 0;

	sprndsini[0] = 0;

	int iproc;
	for (iproc=0;iproc<_nproc;iproc++) {

		temp.nodetype = compnode;
		temp.nodeid = iproc;
		temp.nodecpu = iproc;
		temp.nodelv = nlev;
		temp.fatherid = 0;
		temp.fathercpu = 0;
		temp.fatherid2 = 0;
		temp.fathercpu2 = 0;
		temp.nchilds = 1;
		temp.childs[0] = iproc;

		nodes[iproc] = temp;

	};

	int inodebeg = 0;
	int inodeend = _nproc-1;

	int nnodesloc = _nproc;

	sprndsini[nlev+1] = _nproc;

// Main cycle over the levels

	int inodeendsave;

	while (inodeend > inodebeg) {

		nlev++;

// Combine the sets of cpu's into the groups

		int inodeendl = inodebeg-1;

		inodeendsave = inodeendl;

		while (inodeend > inodeendl) {

// Detect the local group

			int inodebegl = inodeendl+1;
			inodeendl = inodebegl+_nchilds-1;
			if (inodeendl>inodeend) inodeendl = inodeend;
			if (inodeend == inodeendl+1) inodeendl = inodeend;

			int nchildsloc = inodeendl-inodebegl+1;

// Assign data of the group

			inodeendsave = inodeendl;

			if (nchildsloc > 1) {

// Check the memory

				if (nnodesloc >= nnodes-1) throw " Insufficient local number of nodes";

// Create and register new node

				temp.nodetype = compnode;
				temp.nodeid = nnodesloc;
				temp.nodecpu = -1;
				temp.nodelv = nlev;
				temp.fatherid = nnodesloc;
				temp.fathercpu = -1;
				temp.fatherid2 = nnodesloc;
				temp.fathercpu2 = -1;
				temp.nchilds = nchildsloc;
				for (int ichild=0;ichild<nchildsloc;ichild++) {
					int inodeloc = inodebegl+ichild;
					nodes[inodeloc].fatherid = nnodesloc;
					nodes[inodeloc].fatherid2 = nnodesloc;
					temp.childs[ichild] = inodeloc;
				};

				nodes[nnodesloc] = temp;

				nnodesloc++;

			} else {

				throw " Error in combining the childs";

			};

		};

		sprndsini[nlev+1] = nnodesloc;

		inodebeg = inodeendsave+1;
		inodeend = nnodesloc-1;

	};

// Renumber the nodes and the levels of the tree in the inverse order

	nlev++;

	int *order;
	int *orderlv;
	int *iorderlv;

	order = new int [nnodesloc+1];
	if (!order) MemoryFail (funcname);
	orderlv = new int [nlev+1];
	if (!orderlv) MemoryFail (funcname);
	iorderlv = new int [nlev+1];
	if (!iorderlv) MemoryFail (funcname);

	for (i=0;i<nlev;i++) orderlv[i] = nlev-i-1;
	for (i=0;i<nlev;i++) iorderlv[orderlv[i]] = i;

	int nz = 0;

	for (int ilevnew=0;ilevnew<nlev;ilevnew++) {
		int ilevold = iorderlv[ilevnew];
		for (int kii=sprndsini[ilevold];kii<sprndsini[ilevold+1];kii++) {
			order[kii] = nz++;
		};
	};

// Reassign the data of the tree

	CNode *nodesnew;

	nodesnew = new CNode [nnodesloc+1];
	if (!nodesnew) MemoryFail (funcname);

	int inode;
	for (inode=0;inode<nnodesloc;inode++) {
		int inodeloc = nodes[inode].nodeid;
		int inodenew = order[inodeloc];
		nodes[inode].nodeid = inodenew;
		inodeloc = nodes[inode].fatherid;
		inodenew = order[inodeloc];
		nodes[inode].fatherid = inodenew;
		nodes[inode].fatherid2 = inodenew;
		int ilevloc = nodes[inode].nodelv;
		int ilevnew = orderlv[ilevloc];
		nodes[inode].nodelv = ilevnew;
		for (int ichild=0;ichild<nodes[inode].nchilds;ichild++) {
			inodeloc = nodes[inode].childs[ichild];
			inodenew = order[inodeloc];
			nodes[inode].childs[ichild] = inodenew;
		};
	};

	for (inode=0;inode<nnodesloc;inode++) {
		int inodenew = order[inode];
		nodesnew[inodenew] = nodes[inode];
	};

	for (inode=0;inode<nnodesloc;inode++) nodes[inode] = nodesnew[inode];

	delete [] order;
	delete [] orderlv;
	delete [] iorderlv;
	delete [] sprndsini;

	delete [] nodesnew;

// Assign the cpu numbers

	for (iproc=0;iproc<_nproc;iproc++) {
		inodeend = nnodesloc-_nproc+iproc;
		inode = inodeend;
		while (inode != -1) {
			int fathernode = nodes[inode].fatherid;
			if (inode != fathernode) {
				if (nodes[fathernode].childs[0] == inode && nodes[inode].nodecpu == iproc) {
					nodes[fathernode].nodecpu = iproc;
				};
				inode = fathernode;
			} else {
				inode = -1;
			};
		};
	};

	rootid = 0;
	rootcpu = 0;

	nodes[rootid].nodecpu = rootcpu;

	for (inode=0;inode<nnodesloc;inode++) {
		int fathernode = nodes[inode].fatherid;
		int iproc = nodes[fathernode].nodecpu;
		nodes[inode].fathercpu = iproc;
		nodes[inode].fathercpu2 = iproc;
	};

	nnodes = nnodesloc;

// Assign cpuidbeg and cpuidend

	for (iproc=0;iproc<nproc;iproc++) {
		cpuidbeg[iproc] = -1;
		cpuidend[iproc] = -1;
	};
	for (inode=nnodes-1;inode>=0;inode--) {
		iproc = nodes[inode].nodecpu;
		cpuidbeg[iproc] = inode;
	};
	for (inode=0;inode<nnodes;inode++) {
		iproc = nodes[inode].nodecpu;
		cpuidend[iproc] = inode;
	};

// Compute the total performance of each computing node

	for (i=0;i<nnodes;i++) {
		nodes[i].indbeg = -1;
		nodes[i].nodeperf = 0.0e0;
	};

	for (i=0;i<nnodes;i++) {
		if (nodes[i].nchilds == 1) {
			iproc = nodes[i].nodecpu;
			nodes[i].nodeperf = perf2cpu[iproc];
			nodes[i].indbeg = 0;
		};
	};

	int icount = 1;

	while (icount > 0) {
		icount = 0;
		for (i=0;i<nnodes;i++) {
			if (nodes[i].indbeg == -1) {
				icount++;
				int jcount = 0;
				for (int ichild=0;ichild<nodes[i].nchilds;ichild++) {
					int ichildnode = nodes[i].childs[ichild];
					if (nodes[ichildnode].indbeg == -1) jcount++;
				};
				if (jcount == 0) {
					for (int ichild=0;ichild<nodes[i].nchilds;ichild++) {
						int ichildnode = nodes[i].childs[ichild];
						nodes[i].nodeperf += nodes[ichildnode].nodeperf;
					};
					nodes[i].indbeg = 0;
				};
			};
		};
	};

	for (i=0;i<nnodes;i++) {
		nodes[i].indbeg = -1;
	};

};

// Author: Kharchenko S.A.
// CTree: Copy constructor
//========================================================================================
CTree::CTree (const CTree &_tree) { // Copy constructor

	const char *funcname = "CTree_copy";

	myid = _tree.myid;
	nlev = _tree.nlev;
	nnodes = _tree.nnodes;
	nproc = _tree.nproc;
	nproctot = _tree.nproctot;
	rootid = _tree.rootid;
	rootcpu = _tree.rootcpu;
	nsubtree = _tree.nsubtree;
	nodes = new CNode [nnodes+1];
	if (!nodes) MemoryFail (funcname);
	int i;
	for (i=0;i<nnodes;i++) {
		nodes[i] = _tree.nodes[i];
	};
	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	for (i=0;i<nproctot;i++) {
		cpuidarr[i] = _tree.cpuidarr[i];
	};
	cpuidbeg = new int [nproc+1];
	if (!cpuidbeg) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidbeg[i] = _tree.cpuidbeg[i];
	};
	cpuidend = new int [nproc+1];
	if (!cpuidend) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidend[i] = _tree.cpuidend[i];
	};
	perf2cpu = new double [nproc+1];
	if (!perf2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		perf2cpu[i] = _tree.perf2cpu[i];
	};
	memory2cpu = new double [nproc+1];
	for (i=0;i<nproc;i++) {
		memory2cpu[i] = _tree.memory2cpu[i];
	};
	if (nsubtree != 0) {
		subtreearr = new CTree [nsubtree];
		if (!subtreearr) MemoryFail (funcname);
		for (i=0;i<nsubtree;i++) {
			subtreearr[i] = _tree.subtreearr[i];
		};
	} else {
		subtreearr = 0;
	};

	nblkstree = _tree.nblkstree;
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	for (i=0;i<nblkstree+1;i++) {
		blkstree[i] = _tree.blkstree[i];
	};
	for (i=0;i<nblkstree;i++) {
		blk2cputree[i] = _tree.blk2cputree[i];
	};

	ablktree = _tree.ablktree;

	comm = _tree.comm;

};

// Author: Kharchenko S.A.
// CTree: Equality operator
//========================================================================================
CTree &CTree::operator= (const CTree &_tree2) { // Equality operator

	const char *funcname = "CTree_=";

	delete [] nodes;
	delete [] cpuidarr;
	delete [] cpuidbeg;
	delete [] cpuidend;
	delete [] perf2cpu;
	delete [] memory2cpu;
	if (subtreearr != 0) delete [] subtreearr;
	delete [] blkstree;
	delete [] blk2cputree;

	myid = _tree2.myid;
	nlev = _tree2.nlev;
	nnodes = _tree2.nnodes;
	nproc = _tree2.nproc;
	nproctot = _tree2.nproctot;
	rootid = _tree2.rootid;
	rootcpu = _tree2.rootcpu;
	nsubtree = _tree2.nsubtree;
	nodes = new CNode [nnodes+1];
	if (!nodes) MemoryFail (funcname);
	int i;
	for (i=0;i<nnodes;i++) {
		nodes[i] = _tree2.nodes[i];
	};
	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	for (i=0;i<nproctot;i++) {
		cpuidarr[i] = _tree2.cpuidarr[i];
	};
	cpuidbeg = new int [nproc+1];
	if (!cpuidbeg) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidbeg[i] = _tree2.cpuidbeg[i];
	};
	cpuidend = new int [nproc+1];
	if (!cpuidend) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		cpuidend[i] = _tree2.cpuidend[i];
	};
	perf2cpu = new double [nproc+1];
	if (!perf2cpu) MemoryFail (funcname);
	for (i=0;i<nproc;i++) {
		perf2cpu[i] = _tree2.perf2cpu[i];
	};
	memory2cpu = new double [nproc+1];
	for (i=0;i<nproc;i++) {
		memory2cpu[i] = _tree2.memory2cpu[i];
	};
	if (nsubtree != 0) {
		subtreearr = new CTree [nsubtree];
		if (!subtreearr) MemoryFail (funcname);
		for (i=0;i<nsubtree;i++) {
			subtreearr[i] = _tree2.subtreearr[i];
		};
	} else {
		subtreearr = 0;
	};

	nblkstree = _tree2.nblkstree;
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	for (i=0;i<nblkstree+1;i++) {
		blkstree[i] = _tree2.blkstree[i];
	};
	for (i=0;i<nblkstree;i++) {
		blk2cputree[i] = _tree2.blk2cputree[i];
	};

	ablktree = _tree2.ablktree;

	comm = _tree2.comm;

	return *this;

};

// Author: Kharchenko S.A.
// CTree: Destructor
//========================================================================================
CTree::~CTree () { // Destructor
//	cout << " On entry to CTree destructor " << endl;
	myid = -1;
	nlev = 0;
	nnodes = 0;
	nproc = 0;
	nproctot = 0;
	rootid = -1;
	rootcpu = -1;
	delete [] nodes;
	delete [] cpuidarr;
	delete [] cpuidbeg;
	delete [] cpuidend;
	delete [] perf2cpu;
	delete [] memory2cpu;
	if (subtreearr != 0) delete [] subtreearr;
	delete [] blkstree;
	delete [] blk2cputree;
//	cout << " On return from CTree destructor " << endl;
};

// Author: Kharchenko S.A.
// CTree: Set comm for the tree and for the nodes of the tree
//========================================================================================
void CTree::SetCommTreeAndNodes (CMPIComm &_comm) { // Set comm for the tree and for the nodes of the tree

	SetComm (_comm);
	SetCommNodes (_comm);

};

// Author: Kharchenko S.A.
// CTree: Set comm for the tree
//========================================================================================
void CTree::SetCommSchur (CMPIComm &_comm) { // Set comm

	SetCommTreeAndNodes (_comm);

	int myidloc = myid;

	int inode, indtreeloc, i;

	for (inode=0;inode<nnodes;inode++) {
		indtreeloc = nodes[inode].indtree;
		if (indtreeloc >= 0) {

			CMPIComm commloc;
			commloc = _comm.CreateComm (subtreearr[indtreeloc].nproc, subtreearr[indtreeloc].cpuidarr);

			bool is_current_cpu = false;

			for (i=0;i<subtreearr[indtreeloc].nproc;i++) {
				if (subtreearr[indtreeloc].cpuidarr[i] == myidloc) is_current_cpu = true;
			};

			if (is_current_cpu) {
				subtreearr[indtreeloc].SetCommTreeAndNodes (commloc);
			};

		};

	};

};

// Author: Kharchenko S.A.
// CTree: Compute the list of processors for each node of the tree
//========================================================================================
void CTree::CpuListsForNodes () { // Compute the list of processors for each node of the tree

	const char *funcname = "CpuListsForNodes";

// Allocate work arrays

	int *imasknode, *imask, *listloc;

	imasknode = new int [nnodes];
	if (!imasknode) MemoryFail (funcname);
	imask = new int [nproc];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nproc];
	if (!listloc) MemoryFail (funcname);

	int i;

	for (i=0;i<nnodes;i++) imasknode[i] = -1;

	int icycle = -1;

	for (i=0;i<nproc;i++) {
		imask[i] = icycle;
	};

// Main cycle over the tree nodes

	int nlistloc, inode, j, k, ichild, iproc, jproc;

	int nmarked = 0;

	bool ready;

	while (nmarked != nnodes) {
		for (inode=0;inode<nnodes;inode++) {
			if (imasknode[inode] == -1) {
				ready = true;
				if (nodes[inode].nchilds > 1) {
					for (j=0;j<nodes[inode].nchilds;j++) {
						ichild = nodes[inode].childs[j];
						if (imasknode[ichild] == -1) ready = false;
					};
				};
				if (ready) {
					icycle++;
					nlistloc = 0;
					iproc = nodes[inode].nodecpu;
					imask[iproc] = icycle;
					listloc[nlistloc] = iproc;
					nlistloc++;
					if (nodes[inode].nchilds > 1) {
						for (j=0;j<nodes[inode].nchilds;j++) {
							ichild = nodes[inode].childs[j];
							for (k=0;k<nodes[ichild].nprocnode;k++) {
								jproc = nodes[ichild].cpuidarrnode[k];
								if (imask[jproc] != icycle) {
									imask[jproc] = icycle;
									listloc[nlistloc] = jproc;
									nlistloc++;
								};
							};
						};
					};
					qsort (listloc, nlistloc, sizeof(int), compint);
					delete [] nodes[inode].cpuidarrnode;
					delete [] nodes[inode].cpuidarrnodeloc;
					nodes[inode].nprocnode = nlistloc;
					nodes[inode].cpuidarrnode = new int [nlistloc];
					for (j=0;j<nlistloc;j++) nodes[inode].cpuidarrnode[j] = listloc[j];
					nodes[inode].cpuidarrnodeloc = new int [nlistloc];
					for (j=0;j<nlistloc;j++) nodes[inode].cpuidarrnodeloc[j] = listloc[j];
					delete [] nodes[inode].cpuarrndgl2loc;
					nodes[inode].nprocgl = nproc;
					nodes[inode].cpuarrndgl2loc = new int [nproc];
					for (j=0;j<nproc;j++) nodes[inode].cpuarrndgl2loc[j] = -1;
					for (j=0;j<nlistloc;j++) {
						iproc = listloc[j];
						nodes[inode].cpuarrndgl2loc[iproc] = j;
					};
					imasknode[inode] = 1;
					nmarked++;
				};
			};
		};
	};

// Free work arrays

	delete [] imasknode;
	delete [] imask;
	delete [] listloc;

};

// Author: Kharchenko S.A.
// CTree: Compute the list of processors for each node of the tree including Schur subtree's
//========================================================================================
void CTree::CpuListsForNodesSchur () { // Compute the list of processors for each node of the tree including Schur subtree's

//	const char *funcname = "CpuListsForNodesSchur";

	CpuListsForNodes ();

	int inode, indtreeloc;

	for (inode=0;inode<nnodes;inode++) {
		indtreeloc = nodes[inode].indtree;
		if (indtreeloc >= 0) {

			subtreearr[indtreeloc].CpuListsForNodes ();
			subtreearr[indtreeloc].ModifyCpuLists (nodes[inode].cpuidarrnode);

		};
	};

};

// Author: Kharchenko S.A.
// CTree: Modify the proc numbers in the nodes lists
//========================================================================================
void CTree::ModifyCpuLists (int *_cpuarr) { // Modify the proc numbers in the nodes lists

//	const char *funcname = "ModifyCpuLists";

	int inode, j, iproc, iprocnew;

	for (inode=0;inode<nnodes;inode++) {
		for (j=0;j<nodes[inode].nprocnode;j++) {
			iproc = nodes[inode].cpuidarrnode[j];
			iprocnew = _cpuarr[iproc];
			nodes[inode].cpuidarrnode[j] = iprocnew;
		};
	};

};

// Author: Kharchenko S.A.
// CTree: Set comm for the nodes of the tree
//========================================================================================
void CTree::SetCommNodes (CMPIComm &_comm) { // Set comm for the nodes of the tree

	int inode;

	for (inode=0;inode<nnodes;inode++) {
		CMPIComm commloc;
		commloc = _comm.CreateComm (nodes[inode].nprocnode, nodes[inode].cpuidarrnodeloc);
		nodes[inode].SetComm (commloc);
	};

};

// Author: Kharchenko S.A.
// CTree: Extend the tree
//========================================================================================
CTree CTree::ExtendTree () const { // Extend the tree

	const char *funcname = "ExtendTree";

// Create dummy tree

	CTree tree (nproc);

	tree.comm = comm;

// Copy original nodes

	int i;
	for (i=0;i<nnodes;i++) tree.nodes[i] = nodes[i];

// Scan the set of original nodes and create the new ones

	int nnodesloc = nnodes;

	int inode = 0;

	while (inode >= 0) {

		int nchildsloc = nodes[inode].nchilds;

		if (nchildsloc > 1) {

// Create a copy of the subtree

			CNode temp (nchildsloc,nchildsloc);

			temp.nodetype   = exchangenode;
			temp.nodeid     = nnodesloc;
			temp.nodecpu    = nodes[inode].nodecpu;
			temp.fatherid   = inode;
			temp.fathercpu  = nodes[inode].nodecpu;
			temp.fatherid2  = inode;
			temp.fathercpu2 = nodes[inode].nodecpu;

			int ichild;
			for (ichild=0;ichild<nchildsloc;ichild++) {
				int childnodeid = nodes[inode].childs[ichild];
				temp.childs[ichild] = nodes[inode].childs[ichild];
				temp.childs2[ichild] = nnodesloc+ichild+1;

				tree.nodes[childnodeid].fatherid = nnodesloc;
				tree.nodes[childnodeid].fatherid2 = nnodesloc;

			};

			tree.nodes[nnodesloc] = temp;

			int inode2 = nnodesloc;

			nnodesloc++;

			for (ichild=0;ichild<nchildsloc;ichild++) {

				int childnodeid = nodes[inode].childs[ichild];
				int childcpuid  = nodes[childnodeid].nodecpu;

				CNode temp2 (0,0);

				temp2.nodetype   = compnode2;
				temp2.nodeid     = nnodesloc;
				temp2.nodecpu    = childcpuid;
				temp2.fatherid   = inode;
				temp2.fathercpu  = nodes[inode].nodecpu;
				temp2.fatherid2  = inode2;
				temp2.fathercpu2 = nodes[inode].nodecpu;

				tree.nodes[nnodesloc] = temp2;

				tree.nodes[inode].childs[ichild] = nnodesloc;

				nnodesloc++;

			};

		};

		inode++;

		if (inode >= nnodes) inode = -1;

	};

// Reassing the tree levels

	for (i=0;i<nnodesloc;i++) {
		tree.nodes[i].nodelv = -1;
	};

	tree.nodes[rootid].nodelv = 0;

	int icount = 1;

	while (icount > 0) {
		icount = 0;
		for (i=0;i<nnodesloc;i++) {
			if (tree.nodes[i].nodelv == -1) {
				icount++;
				if (tree.nodes[i].nodetype != exchangenode) {
					int fathernodeid = tree.nodes[i].fatherid;
					if (tree.nodes[fathernodeid].nodelv != -1) {
						tree.nodes[i].nodelv = tree.nodes[fathernodeid].nodelv+1;
					};
				} else {
					int fathernodeid = tree.nodes[i].childs2[0];
					if (tree.nodes[fathernodeid].nodelv != -1) {
						tree.nodes[i].nodelv = tree.nodes[fathernodeid].nodelv+1;
					};
				};
			};
		};
	};

	int nlevloc = 0;

	for (i=0;i<nnodesloc;i++) {
		if (tree.nodes[i].nodelv > nlevloc) nlevloc = tree.nodes[i].nodelv;
	};

	nlevloc++;

// Compute cpuidbeg and cpuidend

	for (i=0;i<nproc;i++) tree.cpuidend[i] = cpuidend[i];
	for (i=0;i<nproc;i++) tree.cpuidbeg[i] = cpuidend[i];

	int *ilevarr;

	ilevarr = new int [nproc];
	if (!ilevarr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) {
		inode = cpuidend[i];
		ilevarr[i] = tree.nodes[inode].nodelv;
	};

	for (i=0;i<nnodesloc;i++) {
		int ilevloc = tree.nodes[i].nodelv;
		int iprocloc = tree.nodes[i].nodecpu;
		if (ilevloc < ilevarr[iprocloc]) {
			tree.cpuidbeg[iprocloc] = i;
			ilevarr[iprocloc] = ilevloc;
		};
	};

	delete [] ilevarr;

// Finalize the tree

	tree.myid = myid;
	tree.nlev = nlevloc;
	tree.nnodes = nnodesloc;
	tree.nproc = nproc;
	tree.rootid = rootid;
	tree.rootcpu = rootcpu;

	return tree;

};

// Author: Kharchenko S.A.
// CTree: Collapse back general extended tree
//========================================================================================
CTree CTree::CollapseGeneralTree () const { // Collapse back general extended tree

	const char *funcname = "CollapseGeneralTree";

// Create dummy tree

	CTree tree;

	tree = *this;

// Count the number of nodes

	int nnodesloc = 0;

	int *order;
	int *orderinv;

	order = new int [nnodes+2];
	if (!order) MemoryFail (funcname);
	orderinv = new int [nnodes+2];
	if (!orderinv) MemoryFail (funcname);

	int i;

	for (i=0;i<nnodes;i++) order[i] = -1;
	for (i=0;i<nnodes;i++) orderinv[i] = -1;

	int inode;

	for (inode=0;inode<nnodes;inode++) {
		int inodetype = nodes[inode].nodetype;
		if (inodetype != compnode2 && inodetype != exchangenode) {
			order[inode] = nnodesloc;
			orderinv[nnodesloc] = inode;
			nnodesloc++;
		};
	};

// Scan the set of original nodes redirect necessary ones

	for (inode=0;inode<nnodes;inode++) {
		int inodetype = nodes[inode].nodetype;
		if (inodetype == compnode2) {
			int fathernodeid = nodes[inode].fatherid;
			int fathernodeid2 = nodes[inode].fatherid2;
			int nchildsloc = nodes[fathernodeid].nchilds;
			tree.nodes[fathernodeid].nodetype = compnode;
			for (int ichild=0;ichild<nchildsloc;ichild++) {
				int childnodeid = nodes[fathernodeid2].childs[ichild];
				tree.nodes[childnodeid].nodetype = compnode;
				tree.nodes[childnodeid].fatherid = fathernodeid;
				tree.nodes[childnodeid].fatherid2 = fathernodeid;
				tree.nodes[fathernodeid].childs[ichild] = childnodeid;
			};
		};
	};

// Reassing the tree levels

	for (i=0;i<nnodes;i++) {
		tree.nodes[i].nodelv = -1;
	};

	tree.nodes[rootid].nodelv = 0;

	int icount = nnodesloc-2;

	while (icount >= 0) {
		for (i=0;i<nnodes;i++) {
			int inodetype = nodes[i].nodetype;
			if (tree.nodes[i].nodelv == -1 && inodetype != compnode2 && inodetype != exchangenode) {
				int fathernodeid = tree.nodes[i].fatherid;
				if (tree.nodes[fathernodeid].nodelv != -1) {
					tree.nodes[i].nodelv = tree.nodes[fathernodeid].nodelv+1;
					icount--;
				};
			};
		};
	};

	int nlevloc = 0;

	for (i=0;i<nnodes;i++) {
		if (tree.nodes[i].nodelv > nlevloc) nlevloc = tree.nodes[i].nodelv;
	};

	nlevloc++;

// Recompute new nodes array

	CNode *nodesnew;

	nodesnew = new CNode [nnodesloc+2];
	if (!nodesnew) MemoryFail (funcname);

	int iold;

	for (i=0;i<nnodesloc;i++) {
		iold = orderinv[i];
		nodesnew[i] = tree.nodes[iold];
	};

// Replace node registration data inside

	for (i=0;i<nnodesloc;i++) {
		nodesnew[i].OrderIndices (order);
	};

// Replace nodes array

	delete [] tree.nodes;

	tree.nodes = nodesnew;

// Change tree control information

	int ind, ind0;

	ind = order[rootid];
	tree.rootid = ind;

	for (i=0;i<nproc;i++) {
		ind0 = tree.cpuidend[i];
		ind = order[ind0];
		tree.cpuidend[i] = ind;
		tree.cpuidbeg[i] = ind;
	};

// Recompute cpuidbeg

	int *ilevarr;

	ilevarr = new int [nproc];
	if (!ilevarr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) {
		inode = tree.cpuidend[i];
		ilevarr[i] = tree.nodes[inode].nodelv;
	};

	for (i=0;i<nnodesloc;i++) {
		int ilevloc = tree.nodes[i].nodelv;
		int iprocloc = tree.nodes[i].nodecpu;
		if (ilevloc < ilevarr[iprocloc]) {
			tree.cpuidbeg[iprocloc] = i;
			ilevarr[iprocloc] = ilevloc;
		};
	};

	delete [] ilevarr;

// Finalize the tree

	tree.nlev = nlevloc;
	tree.nnodes = nnodesloc;

// Free work arrays

	delete [] order;
	delete [] orderinv;

	return tree;

};

// Author: Kharchenko S.A.
// CTree: Collapse back extended tree
//========================================================================================
CTree CTree::CollapseTree () const { // Collapse back extended tree

	const char *funcname = "CollapseTree";

// Create dummy tree

	CTree tree;

	tree = *this;

// Determine the number of original nodes

	int nnodesloc = 0;

	int i;
	for (i=0;i<nproc;i++) {
		if (cpuidend[i] > nnodesloc) nnodesloc = cpuidend[i];
	};

	nnodesloc++;

// Scan the set of original nodes redirect necessary ones

	int inode;
	for (inode=0;inode<nnodes;inode++) {
		int inodetype = nodes[inode].nodetype;
		if (inodetype == compnode2) {
			int fathernodeid = nodes[inode].fatherid;
			int fathernodeid2 = nodes[inode].fatherid2;
			int nchildsloc = nodes[fathernodeid].nchilds;
			tree.nodes[fathernodeid].nodetype = compnode;
			for (int ichild=0;ichild<nchildsloc;ichild++) {
				int childnodeid = nodes[fathernodeid2].childs[ichild];
				tree.nodes[childnodeid].nodetype = compnode;
				tree.nodes[childnodeid].fatherid = fathernodeid;
				tree.nodes[childnodeid].fatherid2 = fathernodeid;
				tree.nodes[fathernodeid].childs[ichild] = childnodeid;
			};
		};
	};

// Reassing the tree levels

	for (i=0;i<nnodesloc;i++) {
		tree.nodes[i].nodelv = -1;
	};

	tree.nodes[rootid].nodelv = 0;

	int icount = 1;

	while (icount > 0) {
		icount = 0;
		for (i=0;i<nnodesloc;i++) {
			if (tree.nodes[i].nodelv == -1) {
				icount++;
				if (tree.nodes[i].nodetype != exchangenode) {
					int fathernodeid = tree.nodes[i].fatherid;
					if (tree.nodes[fathernodeid].nodelv != -1) {
						tree.nodes[i].nodelv = tree.nodes[fathernodeid].nodelv+1;
					};
				} else {
					int fathernodeid = tree.nodes[i].childs2[0];
					if (tree.nodes[fathernodeid].nodelv != -1) {
						tree.nodes[i].nodelv = tree.nodes[fathernodeid].nodelv+1;
					};
				};
			};
		};
	};

	int nlevloc = 0;

	for (i=0;i<nnodesloc;i++) {
		if (tree.nodes[i].nodelv > nlevloc) nlevloc = tree.nodes[i].nodelv;
	};

	nlevloc++;

// Compute cpuidbeg and cpuidend

	for (i=0;i<nproc;i++) tree.cpuidend[i] = cpuidend[i];
	for (i=0;i<nproc;i++) tree.cpuidbeg[i] = cpuidend[i];

	int *ilevarr;

	ilevarr = new int [nproc];
	if (!ilevarr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) {
		inode = cpuidend[i];
		ilevarr[i] = tree.nodes[inode].nodelv;
	};

	for (i=0;i<nnodesloc;i++) {
		int ilevloc = tree.nodes[i].nodelv;
		int iprocloc = tree.nodes[i].nodecpu;
		if (ilevloc < ilevarr[iprocloc]) {
			tree.cpuidbeg[iprocloc] = i;
			ilevarr[iprocloc] = ilevloc;
		};
	};

	delete [] ilevarr;

// Finalize the tree

	tree.myid = myid;
	tree.nlev = nlevloc;
	tree.nnodes = nnodesloc;
	tree.nproc = nproc;
	tree.rootid = rootid;
	tree.rootcpu = rootcpu;

	return tree;

};

// Author: Kharchenko S.A.
// CTree: Sort nodes of the tree according to the block ordering assumed in the matrix
//========================================================================================
void CTree::SortNodes () { // Sort nodes of the tree according to the block ordering assumed in the matrix

	const char *funcname = "SortNodes";

// Compute new nodes ordering and its inverse

	int *order;
	int *orderinv;

	order = new int [nnodes+2];
	if (!order) MemoryFail (funcname);
	orderinv = new int [nnodes+2];
	if (!orderinv) MemoryFail (funcname);

	int i;

	for (i=0;i<nnodes;i++) order[i] = -1;
	for (i=0;i<nnodes;i++) orderinv[i] = -1;

// Scan the tree

	int *inodelv, *ichildlv, *nchildslv;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);

// Prepare root node data

	int ilevcurr = 0;

	int rootnode = rootid;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

	int inumber = nnodes;

// Search the tree

	int inodecurr, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichild, nodetypecurr, nchilds2curr, fathernodeid2, ichildnode2, inodeloc;
	int *pchildscurr, *pchilds2curr;

	CTree treenew = TransformTree ();

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		nodetypecurr = treenew.nodes[inodecurr].nodetype;
		nchildscurr = treenew.nodes[inodecurr].nchilds;
		nchilds2curr = treenew.nodes[inodecurr].nchilds2;
		pchildscurr = treenew.nodes[inodecurr].childs;
		pchilds2curr = treenew.nodes[inodecurr].childs2;

		if (nodetypecurr != compnode) {
			throw " Internal error in SortNodes routine";
		};

// Register current node if not yet registered

		if (order[inodecurr] < 0) {
			inumber--;
			order[inodecurr] = inumber;
		};

// Register all secondary nodes (exchange and comp2)

		if (nchilds2curr > 0) {
			for (ichild=nchilds2curr-1;ichild>=0;ichild--) {
				ichildnode2 = pchilds2curr[ichild];
				inumber--;
				if (order[ichildnode2] >= 0) throw " Internal error in SortNodes routine";
				order[ichildnode2] = inumber;
			};
			ichildnode = pchildscurr[0];
			fathernodeid2 = treenew.nodes[ichildnode].fatherid2;
			inumber--;
			if (order[fathernodeid2] >= 0) throw " Internal error in SortNodes routine";
			order[fathernodeid2] = inumber;
		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = pchildscurr[nchildscurr-1];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = nchildscurr-1;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr > 0) {
			if (ilevcurr == 0) {
				throw " Internal error in SortNodes routine";
			};
			ifathernode = treenew.nodes[inodecurr].fatherid;
			ichildcurr--;
			inodeloc = treenew.nodes[ifathernode].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == 0) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Check obtained ordering

	for (i=0;i<nnodes;i++) if (order[i] < 0) throw " Internal error in SortNodes routine";

	for (i=0;i<nnodes;i++) orderinv[order[i]] = i;

// Compute new nodes array

	CNode *nodesnew;

	nodesnew = new CNode [nnodes+2];
	if (!nodesnew) MemoryFail (funcname);

	int iold;

	for (i=0;i<nnodes;i++) {
		iold = orderinv[i];
		nodesnew[i] = nodes[iold];
	};

// Replace node registration data inside

	for (i=0;i<nnodes;i++) {
		nodesnew[i].OrderIndices (order);
	};

// Replace nodes array

	delete [] nodes;

	nodes = nodesnew;

// Change tree control information

	int ind, ind0;

	ind = order[rootid];
	rootid = ind;

	for (i=0;i<nproc;i++) {
		ind0 = cpuidbeg[i];
		ind = order[ind0];
		cpuidbeg[i] = ind;
		ind0 = cpuidend[i];
		ind = order[ind0];
		cpuidend[i] = ind;
	};

// Free work arrays

	delete [] order;
	delete [] orderinv;
	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;

};

// Author: Kharchenko S.A.
// CTree: Transform the tree for easy tree scan
//========================================================================================
CTree CTree::TransformTree () const { // Transform the tree for easy tree scan

	const char *funcname = "TransformTree";

// Create copy of the tree

	CTree treenew = *this;

// Replace tree data in the new tree according to the old tree if necessary

	int i, j, inode, nchildsloc, ichildnode, inode2, ind;

	for (i=0;i<nnodes;i++) {
		inode = i;
		nchildsloc = nodes[inode].nchilds;
		if (nchildsloc > 0) {
			ichildnode = nodes[inode].childs[0];
			if (nodes[ichildnode].nodetype == compnode2) {
				inode2 = nodes[ichildnode].fatherid2;
				delete [] treenew.nodes[inode].childs2;
				treenew.nodes[inode].childs2 = new int [nchildsloc];
				if (!treenew.nodes[inode].childs2) MemoryFail (funcname);
				for (j=0;j<nchildsloc;j++) {
					treenew.nodes[inode].childs[j] = nodes[inode2].childs[j];
					treenew.nodes[inode].childs2[j] = nodes[inode2].childs2[j];
					ind = nodes[inode2].childs[j];
					treenew.nodes[ind].fatherid = inode;
					ind = nodes[inode2].childs2[j];
					treenew.nodes[ind].fatherid = inode;
				};
				treenew.nodes[inode].nchilds2 = nchildsloc;
			};
		};
	};

	return treenew;

};

// Author: Kharchenko S.A.
// CTree: Block partitioning according to the tree
//========================================================================================
void CTree::PartitionTreeSchur (int &_nblks, int *&_blks, int *&_blk2cpu) { // Block partitioning according to the tree

	const char *funcname = "PartitionTreeSchur";

// Estimate from above the total number of blocks

	int nblksext;
	int *blksext, *blk2cpuext;
	int nblksloc;
	int *blksloc, *blk2cpuloc;
	int *nblksarr;
	int **pblksloc, **pblk2cpuloc;

	_nblks = 0;

	PartitionTree (nblksext, blksext, blk2cpuext);

	_nblks += nblksext;

	nblksarr = new int [nsubtree];
	if (!nblksarr) MemoryFail (funcname);
	pblksloc = new int * [nsubtree];
	if (!pblksloc) MemoryFail (funcname);
	pblk2cpuloc = new int * [nsubtree];
	if (!pblk2cpuloc) MemoryFail (funcname);

	int i, j;

	for (i=0;i<nsubtree;i++) {

		subtreearr[i].PartitionTree (nblksloc, blksloc, blk2cpuloc);

		_nblks += nblksloc;

		nblksarr[i] = nblksloc;
		pblksloc[i] = blksloc;
		pblk2cpuloc[i] = blk2cpuloc;

	};

// Scan the tree and put the corresponding description into the arrays

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);
	_blk2cpu = new int [_nblks];
	if (!_blk2cpu) MemoryFail (funcname);

	_blks[0] = 0;
	_nblks = 0;

	int inode, indtreeloc, nblks0;

	for (inode=0;inode<nnodes;inode++) {
		indtreeloc = nodes[inode].indtree;
		if (indtreeloc < 0) {
			_blks[_nblks+1] = _blks[_nblks] + (blksext[inode+1]-blksext[inode]);
			nodes[inode].indbeg = _nblks;
			nodes[inode].indend = _nblks;
			nodes[inode].indbegtot = _nblks;
			nodes[inode].indendtot = _nblks;
			_blk2cpu[_nblks] = nodes[inode].nodecpu;
			_nblks++;
		} else {
			nblks0 = _nblks;
			for (j=0;j<nblksarr[indtreeloc];j++) {
				_blks[_nblks+1] = _blks[_nblks] + (pblksloc[indtreeloc][j+1]-pblksloc[indtreeloc][j]);
				subtreearr[indtreeloc].nodes[j].indbeg = _nblks;
				subtreearr[indtreeloc].nodes[j].indend = _nblks;
				subtreearr[indtreeloc].nodes[j].indbegtot = _nblks;
				subtreearr[indtreeloc].nodes[j].indendtot = _nblks;
				_blk2cpu[_nblks] = pblk2cpuloc[indtreeloc][j];
				_nblks++;
			};
			nodes[inode].indbeg = nblks0;
			nodes[inode].indend = _nblks-1;
			nodes[inode].indbegtot = nblks0;
			nodes[inode].indendtot = _nblks-1;
		};
	};

// Free work arrays

	delete [] blksext;
	delete [] blk2cpuext;

	for (i=0;i<nsubtree;i++) {

		delete [] pblksloc[i];
		delete [] pblk2cpuloc[i];

	};

	delete [] nblksarr;
	delete [] pblksloc;
	delete [] pblk2cpuloc;

};

// Author: Kharchenko S.A.
// CTree: Block partitioning according to the tree
//========================================================================================
void CTree::PartitionTree (int &_nblks, int *&_blks, int *&_blk2cpu) { // Block partitioning according to the tree

	const char *funcname = "PartitionTree";

	_blks = new int [nnodes+2];
	if (!_blks) MemoryFail (funcname);
	_blk2cpu = new int [nnodes+2];
	if (!_blk2cpu) MemoryFail (funcname);

	_blks[0] = 0;
	_nblks = 0;

	int inode;

	for (inode=0;inode<nnodes;inode++) {
		if (nodes[inode].nodetype==compnode || nodes[inode].nodetype==compnode2) {
			_blks[_nblks+1] = nodes[inode].indend+1;
			nodes[inode].indbeg = _nblks;
			nodes[inode].indend = _nblks;
			_blk2cpu[_nblks] = nodes[inode].nodecpu;
			_nblks++;
		};
	};

	int fathernode;
	int ichildnode;

	for (inode=0;inode<nnodes;inode++) {
		if (nodes[inode].nodetype==compnode || nodes[inode].nodetype==compnode2) {
			fathernode = nodes[inode].fatherid;
			ichildnode = nodes[inode].childs[0];
//			ichild2node = nodes[ichildnode].childs[0];
//			nodes[inode].indbegtot = nodes[ichildnode].indbeg;
//			nodes[inode].indendtot = nodes[fathernode].indend;
//			nodes[inode].indbegtot = nodes[ichildnode].indbeg;
//			nodes[inode].indendtot = nodes[inode].indend;
			nodes[inode].indbegtot = nodes[inode].indbeg;
			nodes[inode].indendtot = nodes[inode].indend;
		};
	};

};

// Author: Kharchenko S.A.
// CTree: Store block partitioning in a tree
//========================================================================================
void CTree::StorePartitioning (int _nblks, int *_blks, int *_blk2cpu) { // Store block partitioning in a tree

	const char *funcname = "PartitionTree";

	delete [] blkstree;
	delete [] blk2cputree;

	nblkstree = _nblks;

	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	int i;

	for (i=0;i<=nblkstree;i++) blkstree[i] = _blks[i];
	for (i=0;i<nblkstree;i++) blk2cputree[i] = _blk2cpu[i];

};

// Author: Kharchenko S.A.
// CTree: Include one extended binary tree into smaller one with inclusion of the subtrees
//========================================================================================
void CTree::IncludeBinaryTreeSchur (CTree &_treesm, int _nblks, int *_blks, CSMatrix &_amatr) { // Include one extended binary tree into smaller one

	const char *funcname = "IncludeBinaryTreeSchur";

// Perform upper level inclusion

	IncludeBinaryTree (_treesm);

// Count the number of subtrees to be generated and stored

	int nsubtreeloc = 0;

	int inode, inodeext;

	for (inode=0;inode<_treesm.nnodes;inode++) {
		if (_treesm.nodes[inode].nchilds > 1) {
			inodeext = _treesm.nodes[inode].nodeidext;
			if (nodes[inodeext].indtree >= 0) nsubtreeloc++;
		};
	};

	CTree *psubtreeloc = 0;

	if (nsubtreeloc > 0) {
		psubtreeloc = new CTree [nsubtreeloc];
		if (!psubtreeloc) MemoryFail (funcname);
		_treesm.SetNsubtree (nsubtreeloc);
		_treesm.SetSubtreearr (psubtreeloc);
	};

// Allocate list of processors numbers

	int nproctot = _treesm.GetNproc ();

	int *imaskcpu;
	int *listcpu;

	imaskcpu = new int [nproctot];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproctot];
	if (!listcpu) MemoryFail (funcname);

	int icycle = -1;

	int i;

	for (i=0;i<nproctot;i++) imaskcpu[i] = icycle;

// Prepare arrays for scanning the tree

	int nlevloc = _treesm.GetNlev ();

	int *inodelv, *ichildlv, *nchildslv;

	inodelv = new int [nlevloc+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlevloc+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlevloc+1];
	if (!nchildslv) MemoryFail (funcname);

// Scan all nodes of the smaller tree and generate and store if necessary the corresponding subtrees

	int indsubtreeloc, inodeloc, ilev, jnode;
	int iproc;

	nsubtreeloc = 0;

	for (inode=0;inode<_treesm.nnodes;inode++) {
		if (_treesm.nodes[inode].nchilds > 1) {
			inodeext = _treesm.nodes[inode].nodeidext;
			indsubtreeloc = nodes[inodeext].indtree;
			if (indsubtreeloc >= 0) {

// Count the number of free processors available for the subtree and store all of them in the list

				int nprocfreetot = 0;

				icycle++;

// Find the list

				bool is_node_found;

				for (jnode=0;jnode<_treesm.nnodes;jnode++) {
					is_node_found = false;
					ilev = _treesm.nodes[jnode].nodelv;
					inodeloc = jnode;
					if (inodeloc == inode) is_node_found = true;
					while (ilev > 0) {
						if (inodeloc == inode) is_node_found = true;
						inodeloc = _treesm.nodes[inodeloc].fatherid;
						ilev = _treesm.nodes[inodeloc].nodelv;
					};
					if (inodeloc == inode) is_node_found = true;
					if (is_node_found) {
						iproc = _treesm.nodes[jnode].nodecpu;
						if (imaskcpu[iproc] != icycle) {
							listcpu[nprocfreetot] = iproc;
							nprocfreetot++;
							imaskcpu[iproc] = icycle;
						};
					};
				};

				qsort (listcpu,nprocfreetot,sizeof(int),compint);

// Check the number of processors in the extended subtree

				int nprocfree = nprocfreetot;

				int nproclocext = subtreearr[indsubtreeloc].nproc;

				if (nprocfree > nproclocext) nprocfree = nproclocext;

// Create new local binary tree

				int nchilds = 2;

				double *memory2cpu;
				double *procweight;

				memory2cpu = new double [nprocfree+1];
				if (!memory2cpu) MemoryFail (funcname);
				procweight = new double [nprocfree+1];
				if (!procweight) MemoryFail (funcname);

				for (i=0;i<nprocfree;i++) memory2cpu[i] = 1.0e0;
				for (i=0;i<nprocfree;i++) procweight[i] = 1.0e0;

				CTree treeloc (nprocfree, nchilds, procweight, memory2cpu);

				delete [] memory2cpu;
				delete [] procweight;

				treeloc.SortNodes ();

// Include binary tree into the extended one

				subtreearr[indsubtreeloc].IncludeBinaryTree (treeloc);

// Store CPU ID numbers in transformation arrays

				int *pcpuidarr = treeloc.GetCpuidarr ();

				delete [] pcpuidarr;

				pcpuidarr = new int [nprocfreetot];
				if (!pcpuidarr) MemoryFail (funcname);

				for (i=0;i<nprocfreetot;i++) {
					pcpuidarr[i] = listcpu[i];
				};

				treeloc.SetCpuidarr (pcpuidarr);

				treeloc.SetNproctot (nprocfreetot);

// Store and register new subtree

				_treesm.nodes[inode].indtree = nsubtreeloc;
				psubtreeloc[nsubtreeloc] = treeloc;

				nsubtreeloc++;

			};
		};
	};

// Free work arrays

	delete [] imaskcpu;
	delete [] listcpu;

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;

};

// Author: Kharchenko S.A.
// CTree: Include one extended binary tree into smaller one
//========================================================================================
void CTree::IncludeBinaryTree (CTree &_treesm) { // Include one extended binary tree into smaller one

	const char *funcname = "IncludeBinaryTree";

// Prepare root node data

	int nlevtot = nlev;

	int *inodelv, *inodelvext, *ichildlv, *nchildslv;

	inodelv = new int [nlevtot+1];
	if (!inodelv) MemoryFail (funcname);
	inodelvext = new int [nlevtot+1];
	if (!inodelvext) MemoryFail (funcname);
	ichildlv = new int [nlevtot+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlevtot+1];
	if (!nchildslv) MemoryFail (funcname);

	int ilevcurr = 0;

	int rootnodeext = rootid;
	int rootnodesm = _treesm.rootid;

	inodelv[ilevcurr] = rootnodesm;
	inodelvext[ilevcurr] = rootnodeext;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	int inodecurr, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichildnodeext, ifathernodeext, inode, inodecurrext, inodelocext;
	int ichild, inodeloc;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];
		inodecurrext = inodelvext[ilevcurr];

		nchildscurr = _treesm.nodes[inodecurr].nchilds;

// Register tree data

		if (nchildscurr != 1) {
			_treesm.nodes[inodecurr].indbeg = nodes[inodecurrext].indbeg;
			_treesm.nodes[inodecurr].indend = nodes[inodecurrext].indend;
			_treesm.nodes[inodecurr].nodeidext = inodecurrext;
			nodes[inodecurrext].nodeidsm = inodecurr;
		} else {
			ichild = nodes[inodecurrext].childs[0];
			nodes[ichild].nodeidsm = inodecurr;
			while (nodes[ichild].nchilds > 1) {
				ichild = nodes[ichild].childs[0];
				nodes[ichild].nodeidsm = inodecurr;
			};
			_treesm.nodes[inodecurr].indbeg = nodes[ichild].indbeg;
			_treesm.nodes[inodecurr].indend = nodes[inodecurrext].indend;
			_treesm.nodes[inodecurr].nodeidext = inodecurrext;
		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _treesm.nodes[inodecurr].childs[0];
			ichildnodeext = nodes[inodecurrext].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			inodelvext[ilevcurr] = ichildnodeext;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in CTree::IncludeBinaryTree routine";
			};
			ifathernode = _treesm.nodes[inodecurr].fatherid;
			ifathernodeext = nodes[inodecurrext].fatherid;
			ichildcurr ++;
			inodeloc = _treesm.nodes[ifathernode].childs[ichildcurr];
			inodelocext = nodes[ifathernodeext].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			inodelvext[ilevcurr] = inodelocext;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		inodecurrext = inodelvext[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

	int fathernode;

	for (inode=0;inode<_treesm.nnodes;inode++) {
		if (_treesm.nodes[inode].nodetype==compnode || _treesm.nodes[inode].nodetype==compnode2) {
			fathernode = _treesm.nodes[inode].fatherid;
			ichildnode = _treesm.nodes[inode].childs[0];
//			_treesm.nodes[inode].indbegtot = _treesm.nodes[ichildnode].indbeg;
//			_treesm.nodes[inode].indendtot = _treesm.nodes[inode].indend;
			_treesm.nodes[inode].indbegtot = _treesm.nodes[inode].indbeg;
			_treesm.nodes[inode].indendtot = _treesm.nodes[inode].indend;
		};
	};

// Free work arrays

	delete [] inodelv;
	delete [] inodelvext;
	delete [] ichildlv;
	delete [] nchildslv;

};

// Author: Kharchenko S.A.
// CTree: Include one extended binary tree into smaller one
//========================================================================================
void CTree::IncludeBinaryTree (int _nproc, CTree &_treeext) { // Include one extended binary tree into smaller one

	const char *funcname = "IncludeBinaryTree";

// Create local binary tree

	int nprocloc = _nproc;
	int nchilds = 2;

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nprocloc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocloc+1];
	if (!procweight) MemoryFail (funcname);

	int i;

	for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

	CTree treeloc (nprocloc, nchilds, procweight, memory2cpu);

	delete [] memory2cpu;
	delete [] procweight;

	treeloc.SortNodes ();

	treeloc.CpuListsForNodes ();

// Copy major data

	int nlevtot = treeloc.GetNlev ();
	int nnodesloc = treeloc.GetNnodes ();

	int *inodelv, *inodelvext, *ichildlv, *nchildslv;
	int *inoden2o;

	inodelv = new int [nlevtot+1];
	if (!inodelv) MemoryFail (funcname);
	inodelvext = new int [nlevtot+1];
	if (!inodelvext) MemoryFail (funcname);
	ichildlv = new int [nlevtot+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlevtot+1];
	if (!nchildslv) MemoryFail (funcname);
	inoden2o = new int [nnodesloc];
	if (!inoden2o) MemoryFail (funcname);

	int ilevcurr = 0;

	int rootnodeext = _treeext.rootid;
	int rootnodesm = treeloc.rootid;

	inodelv[ilevcurr] = rootnodesm;
	inodelvext[ilevcurr] = rootnodeext;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

	for (i=0;i<nnodesloc;i++) inoden2o[i] = -1;

	inoden2o[rootnodesm] = rootnodeext;

// Search the tree

	int inodecurr, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichildnodeext, ifathernodeext, inodecurrext, inodelocext;
	int ichild, inodeloc, ichildext;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];
		inodecurrext = inoden2o[inodecurr];

		nchildscurr = treeloc.nodes[inodecurr].nchilds;

// Register tree data

		if (nchildscurr != 1) {
			treeloc.nodes[inodecurr].indbeg = _treeext.nodes[inodecurrext].indbegtot;
			treeloc.nodes[inodecurr].indend = _treeext.nodes[inodecurrext].indendtot;
			treeloc.nodes[inodecurr].indbegtot = _treeext.nodes[inodecurrext].indbegtot;
			treeloc.nodes[inodecurr].indendtot = _treeext.nodes[inodecurrext].indendtot;
			for (i=0;i<nchildscurr;i++) {
				ichild = treeloc.nodes[inodecurr].childs[i];
				ichildext = _treeext.nodes[inodecurrext].childs[i];
				inoden2o[ichild] = ichildext;
			};
		} else {
			treeloc.nodes[inodecurr].indend = _treeext.nodes[inodecurrext].indendtot;
			treeloc.nodes[inodecurr].indendtot = _treeext.nodes[inodecurrext].indendtot;
			int inodeextloc = inodecurrext;
			while (_treeext.nodes[inodeextloc].nchilds > 1) {
				int ichildextloc = _treeext.nodes[inodeextloc].childs[0];
				inodeextloc = ichildextloc;
			};
			treeloc.nodes[inodecurr].indbeg = _treeext.nodes[inodeextloc].indbeg;
//			treeloc.nodes[inodecurr].indbegtot = _treeext.nodes[inodecurrext].indbegtot;
		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = treeloc.nodes[inodecurr].childs[0];
			ichildnodeext = inoden2o[ichildnode];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			inodelvext[ilevcurr] = ichildnodeext;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in CTree::IncludeBinaryTreeSchur routine";
			};
			ifathernode = treeloc.nodes[inodecurr].fatherid;
			ifathernodeext = _treeext.nodes[inodecurrext].fatherid;
			ichildcurr++;
			inodeloc = treeloc.nodes[ifathernode].childs[ichildcurr];
			inodelocext = _treeext.nodes[ifathernodeext].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			inodelvext[ilevcurr] = inodelocext;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		inodecurrext = inodelvext[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Reassign indbegtot

	for (inodeloc=0;inodeloc<treeloc.nnodes;inodeloc++) {
		int inode1loc = inodeloc;
		while (treeloc.nodes[inode1loc].nchilds > 1) {
			int ichild1loc = treeloc.nodes[inode1loc].childs[0];
			inode1loc = ichild1loc;
		};
		treeloc.nodes[inodeloc].indbegtot = treeloc.nodes[inode1loc].indbeg;
	};

// Free work arrays

	delete [] inodelv;
	delete [] inodelvext;
	delete [] ichildlv;
	delete [] nchildslv;
	delete [] inoden2o;

// Copy the block sparsity matrix

	CSMatrix *pablkstrext = _treeext.GetAblktree ();
	CSMatrix *pablkstrloc = treeloc.GetAblktree ();

	*pablkstrloc = *pablkstrext;

// Copy blocks partitioning

	int nblksext = _treeext.GetNblkstree ();
	int *pblksext = _treeext.GetBlkstree ();
	int *pblksloc = treeloc.GetBlkstree ();

	delete [] pblksloc;

	pblksloc = new int [nblksext+1];
	if (!pblksloc) MemoryFail (funcname);

	for (i=0;i<=nblksext;i++) pblksloc[i] = pblksext[i];

// Recompute new data distribution

	int nprocext = _treeext.GetNproc ();

	int isize = nprocext / _nproc;

	if (isize <= 0) throw " CTree::IncludeBinaryTreeSchur: wrong number of processors ";

	if (nprocext%_nproc != 0) throw " CTree::IncludeBinaryTreeSchur: wrong number of processors ";

	int *cpuo2n;

	cpuo2n = new int [nprocext];
	if (!cpuo2n) MemoryFail (funcname);

	int nz = 0;

	int j;

	for (i=0;i<_nproc;i++) {
		for (j=0;j<isize;j++) {
			cpuo2n[nz] = i;
			nz++;
		};
	};

	int *pblks2cpuext = _treeext.GetBlk2cputree ();
	int *pblks2cpuloc = treeloc.GetBlk2cputree ();

	delete [] pblks2cpuloc;

	pblks2cpuloc = new int [nblksext];
	if (!pblks2cpuloc) MemoryFail (funcname);

	int iproc, iprocnew;

	for (i=0;i<nblksext;i++) {
		iproc = pblks2cpuext[i];
		iprocnew = cpuo2n[iproc];
		pblks2cpuloc[i] = iprocnew;
	};

	treeloc.SetNblkstree (nblksext);
	treeloc.SetBlkstree (pblksloc);
	treeloc.SetBlk2cputree (pblks2cpuloc);

// Copy the tree

	*this = treeloc;

// Free work arrays

	delete [] cpuo2n;

};

// Author: Kharchenko S.A.
// CTree: Include one extended binary tree into smaller one
//========================================================================================
void CTree::IncludeBinaryTreeSchur (int _nproc, CTree &_treeext) { // Include one extended binary tree into smaller one

	const char *funcname = "IncludeBinaryTreeSchur";

// Create local binary tree

	int nprocloc = _nproc;
	int nchilds = 2;

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nprocloc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocloc+1];
	if (!procweight) MemoryFail (funcname);

	int i;

	for (i=0;i<nprocloc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocloc;i++) procweight[i] = 1.0e0;

	CTree treeloc (nprocloc, nchilds, procweight, memory2cpu);

	delete [] memory2cpu;
	delete [] procweight;

	treeloc.SortNodes ();

	treeloc.CpuListsForNodes ();

// Copy major data

	int nlevtot = treeloc.GetNlev ();
	int nnodesloc = treeloc.GetNnodes ();

	int *inodelv, *inodelvext, *ichildlv, *nchildslv;
	int *inoden2o;

	inodelv = new int [nlevtot+1];
	if (!inodelv) MemoryFail (funcname);
	inodelvext = new int [nlevtot+1];
	if (!inodelvext) MemoryFail (funcname);
	ichildlv = new int [nlevtot+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlevtot+1];
	if (!nchildslv) MemoryFail (funcname);
	inoden2o = new int [nnodesloc];
	if (!inoden2o) MemoryFail (funcname);

	int ilevcurr = 0;

	int rootnodeext = _treeext.rootid;
	int rootnodesm = treeloc.rootid;

	inodelv[ilevcurr] = rootnodesm;
	inodelvext[ilevcurr] = rootnodeext;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

	for (i=0;i<nnodesloc;i++) inoden2o[i] = -1;

	inoden2o[rootnodesm] = rootnodeext;

// Search the tree

	int inodecurr, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int ichildnodeext, ifathernodeext, inodecurrext, inodelocext;
	int ichild, inodeloc, ichildext;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];
		inodecurrext = inoden2o[inodecurr];

		nchildscurr = treeloc.nodes[inodecurr].nchilds;

// Register tree data

		if (nchildscurr != 1) {
			treeloc.nodes[inodecurr].indbeg = _treeext.nodes[inodecurrext].indbeg;
			treeloc.nodes[inodecurr].indend = _treeext.nodes[inodecurrext].indend;
			treeloc.nodes[inodecurr].indbegtot = _treeext.nodes[inodecurrext].indbegtot;
			treeloc.nodes[inodecurr].indendtot = _treeext.nodes[inodecurrext].indendtot;
			for (i=0;i<nchildscurr;i++) {
				ichild = treeloc.nodes[inodecurr].childs[i];
				ichildext = _treeext.nodes[inodecurrext].childs[i];
				inoden2o[ichild] = ichildext;
			};
		} else {
			treeloc.nodes[inodecurr].indbeg = _treeext.nodes[inodecurrext].indbegtot;
			treeloc.nodes[inodecurr].indend = _treeext.nodes[inodecurrext].indendtot;
			treeloc.nodes[inodecurr].indbegtot = _treeext.nodes[inodecurrext].indbegtot;
			treeloc.nodes[inodecurr].indendtot = _treeext.nodes[inodecurrext].indendtot;
		};

// Go down level

		if (nchildscurr > 1) {

			ichildnode = treeloc.nodes[inodecurr].childs[0];
			ichildnodeext = inoden2o[ichildnode];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			inodelvext[ilevcurr] = ichildnodeext;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " Internal error in CTree::IncludeBinaryTreeSchur routine";
			};
			ifathernode = treeloc.nodes[inodecurr].fatherid;
			ifathernodeext = _treeext.nodes[inodecurrext].fatherid;
			ichildcurr++;
			inodeloc = treeloc.nodes[ifathernode].childs[ichildcurr];
			inodelocext = _treeext.nodes[ifathernodeext].childs[ichildcurr];
			inodelv[ilevcurr] = inodeloc;
			inodelvext[ilevcurr] = inodelocext;
			ichildlv[ilevcurr] = ichildcurr;
			goto exitcycle;
		};

// Go up level if necessary

uplevel:;

		ilevcurr--;

		if (ilevcurr<0) goto exitcycle;

// Search the next child at the current level

		inodecurr = inodelv[ilevcurr];
		inodecurrext = inodelvext[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Free work arrays

	delete [] inodelv;
	delete [] inodelvext;
	delete [] ichildlv;
	delete [] nchildslv;
	delete [] inoden2o;

// Copy the block sparsity matrix

	CSMatrix *pablkstrext = _treeext.GetAblktree ();
	CSMatrix *pablkstrloc = treeloc.GetAblktree ();

	*pablkstrloc = *pablkstrext;

// Copy blocks partitioning

	int nblksext = _treeext.GetNblkstree ();
	int *pblksext = _treeext.GetBlkstree ();
	int *pblksloc = treeloc.GetBlkstree ();

	delete [] pblksloc;

	pblksloc = new int [nblksext+1];
	if (!pblksloc) MemoryFail (funcname);

	for (i=0;i<=nblksext;i++) pblksloc[i] = pblksext[i];

// Recompute new data distribution

	int nprocext = _treeext.GetNproc ();

	int isize = nprocext / _nproc;

	if (isize <= 0) throw " CTree::IncludeBinaryTreeSchur: wrong number of processors ";

	if (nprocext%_nproc != 0) throw " CTree::IncludeBinaryTreeSchur: wrong number of processors ";

	int *cpuo2n;

	cpuo2n = new int [nprocext];
	if (!cpuo2n) MemoryFail (funcname);

	int nz = 0;

	int j;

	for (i=0;i<_nproc;i++) {
		for (j=0;j<isize;j++) {
			cpuo2n[nz] = i;
			nz++;
		};
	};

	int *pblks2cpuext = _treeext.GetBlk2cputree ();
	int *pblks2cpuloc = treeloc.GetBlk2cputree ();

	delete [] pblks2cpuloc;

	pblks2cpuloc = new int [nblksext];
	if (!pblks2cpuloc) MemoryFail (funcname);

	int iproc, iprocnew;

	for (i=0;i<nblksext;i++) {
		iproc = pblks2cpuext[i];
		iprocnew = cpuo2n[iproc];
		pblks2cpuloc[i] = iprocnew;
	};

	treeloc.SetNblkstree (nblksext);
	treeloc.SetBlkstree (pblksloc);
	treeloc.SetBlk2cputree (pblks2cpuloc);

// Copy the tree

	*this = treeloc;

// Free work arrays

	delete [] cpuo2n;

};

// Author: Kharchenko S.A.
// CTree: Block partitioning according to the tree
//========================================================================================
void CTree::Partition (int &_nblks, int *&_blks) const { // Block partitioning according to the tree

	const char *funcname = "Partition";

	_blks = new int [nnodes+2];
	if (!_blks) MemoryFail (funcname);

	_blks[0] = 0;
	_nblks = 0;
	for (int inode=0;inode<nnodes;inode++) {
		if (nodes[inode].nodetype==compnode || nodes[inode].nodetype==compnode2) {
			_blks[_nblks+1] = nodes[inode].indend+1;
			_nblks++;
		};
	};

	qsort (_blks, _nblks+1, sizeof(int), compint);

};

// Author: Kharchenko S.A.
// CTree: Modify block partitioning, return modified block partitioning and return tree in terms of blocks
//========================================================================================
void CTree::TreeInBlocks (int _nsupmx, int &_nblks, int *&_blks) { // Modify block partitioning, return modified block partitioning and return tree in terms of blocks

	const char *funcname = "TreeInBlocks";

// Take current blocks partitioning

	int nblksl;
	int *blksl;

	Partition (nblksl, blksl);

	int nsuptotal = blksl[nblksl];

// Modify block partitioning

	int *blksnew;

	blksnew = new int [nsuptotal+1];
	if (!blksnew) MemoryFail (funcname);

	_nblks = 0;

	blksnew[0] = 0;

	int iblk;
	for (iblk=0;iblk<nblksl;iblk++) {
		int isupbeg = blksl[iblk];
		int isupend = blksl[iblk+1]-1;

		int nsupl = isupend-isupbeg+1;

		int nblki = (nsupl+_nsupmx-1)/_nsupmx;

//		int ni = nsupl / nblki;
		int ni = _nsupmx;

		for (int iblkl=0;iblkl<nblki-1;iblkl++) {
			blksnew[_nblks+1]=blksnew[_nblks]+ni;
			_nblks++;
		};
		blksnew[_nblks+1] = isupend+1;
		_nblks++;
	};

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);

	for (iblk=0;iblk<=_nblks;iblk++) _blks[iblk] = blksnew[iblk];

// Check the result

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blks[iblk]>=_blks[iblk+1]) {
			throw " Error in modified blocks partitioning";
		};
	};

// Recompute data in the tree

	for (int inode=0;inode<nnodes;inode++) {
		int isupbeg = nodes[inode].indbeg;
		int isupend = nodes[inode].indend;
		int iblkbeg = -1;
		int iblkend = -1;
		int isupbegtot = nodes[inode].indbegtot;
		int isupendtot = nodes[inode].indendtot;
		int iblkbegtot = -1;
		int iblkendtot = -1;
		if (nodes[inode].nodetype != exchangenode) {
			for (iblk=0;iblk<_nblks;iblk++) {
				if (_blks[iblk] == isupbeg) iblkbeg = iblk;
				if (_blks[iblk+1]-1 == isupend) iblkend = iblk;
				if (_blks[iblk] == isupbegtot) iblkbegtot = iblk;
				if (_blks[iblk+1]-1 == isupendtot) iblkendtot = iblk;
			};
			if (iblkbeg == -1 || iblkend == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			if (iblkbegtot == -1 || iblkendtot == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			nodes[inode].indbeg = iblkbeg;
			nodes[inode].indend = iblkend;
			nodes[inode].indbegtot = iblkbegtot;
			nodes[inode].indendtot = iblkendtot;
		} else {
			for (iblk=0;iblk<_nblks;iblk++) {
				if (_blks[iblk] == isupbegtot) iblkbegtot = iblk;
				if (_blks[iblk+1]-1 == isupendtot) iblkendtot = iblk;
			};
			if (iblkbegtot == -1 || iblkendtot == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			nodes[inode].indbegtot = iblkbegtot;
			nodes[inode].indendtot = iblkendtot;
		};
	};

// Free work arrays

	delete [] blksl;
	delete [] blksnew;

};

// Author: Kharchenko S.A.
// CTree: Modify block partitioning introduced on entry, return modified block partitioning and return tree in terms of blocks
//========================================================================================
void CTree::TreeInBlocksExt (int _nsupmx, int &_nblks, int *&_blks) { // Modify block partitioning introduced on entry, return modified block partitioning and return tree in terms of blocks

	const char *funcname = "TreeInBlocksExt";

// Take current blocks partitioning

	int nblksl;
	int *blksl;

	nblksl = _nblks;
	blksl = new int [nblksl+1];
	if (!blksl) MemoryFail (funcname);

	int iblk;
	for (iblk=0;iblk<=nblksl;iblk++) blksl[iblk] = _blks[iblk];

	delete [] _blks;

	int nsuptotal = blksl[nblksl];

// Modify block partitioning

	int *blksnew;

	blksnew = new int [nsuptotal+1];
	if (!blksnew) MemoryFail (funcname);

	_nblks = 0;

	blksnew[0] = 0;

	for (iblk=0;iblk<nblksl;iblk++) {
		int isupbeg = blksl[iblk];
		int isupend = blksl[iblk+1]-1;

		int nsupl = isupend-isupbeg+1;

		int nblki = (nsupl+_nsupmx-1)/_nsupmx;

//		int ni = nsupl / nblki;
		int ni = _nsupmx;

		for (int iblkl=0;iblkl<nblki-1;iblkl++) {
			blksnew[_nblks+1]=blksnew[_nblks]+ni;
			_nblks++;
		};
		blksnew[_nblks+1] = isupend+1;
		_nblks++;
	};

	_blks = new int [_nblks+1];
	if (!_blks) MemoryFail (funcname);

	for (iblk=0;iblk<=_nblks;iblk++) _blks[iblk] = blksnew[iblk];

// Check the result

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blks[iblk]>=_blks[iblk+1]) {
			throw " Error in modified blocks partitioning";
		};
	};

// Recompute data in the tree

	for (int inode=0;inode<nnodes;inode++) {
		int isupbeg = nodes[inode].indbeg;
		int isupend = nodes[inode].indend;
		int iblkbeg = -1;
		int iblkend = -1;
		int isupbegtot = nodes[inode].indbegtot;
		int isupendtot = nodes[inode].indendtot;
		int iblkbegtot = -1;
		int iblkendtot = -1;
		if (nodes[inode].nodetype != exchangenode) {
			for (iblk=0;iblk<_nblks;iblk++) {
				if (_blks[iblk] == isupbeg) iblkbeg = iblk;
				if (_blks[iblk+1]-1 == isupend) iblkend = iblk;
				if (_blks[iblk] == isupbegtot) iblkbegtot = iblk;
				if (_blks[iblk+1]-1 == isupendtot) iblkendtot = iblk;
			};
			if (iblkbeg == -1 || iblkend == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			if (iblkbegtot == -1 || iblkendtot == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			nodes[inode].indbeg = iblkbeg;
			nodes[inode].indend = iblkend;
			nodes[inode].indbegtot = iblkbegtot;
			nodes[inode].indendtot = iblkendtot;
		} else {
			for (iblk=0;iblk<_nblks;iblk++) {
				if (_blks[iblk] == isupbegtot) iblkbegtot = iblk;
				if (_blks[iblk+1]-1 == isupendtot) iblkendtot = iblk;
			};
			if (iblkbegtot == -1 || iblkendtot == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			nodes[inode].indbegtot = iblkbegtot;
			nodes[inode].indendtot = iblkendtot;
		};
	};

// Free work arrays

	delete [] blksl;
	delete [] blksnew;

};

// Author: Kharchenko S.A.
// CTree: Return tree in terms of blocks
//========================================================================================
void CTree::TreeInBlocks (int _nblks, int *_blks) { // Return tree in terms of blocks

//	const char *funcname = "TreeInBlocks";

// Recompute data in the tree

	int iblk;

	for (int inode=0;inode<nnodes;inode++) {
		int isupbeg = nodes[inode].indbeg;
		int isupend = nodes[inode].indend;
		int iblkbeg = -1;
		int iblkend = -1;
		int isupbegtot = nodes[inode].indbegtot;
		int isupendtot = nodes[inode].indendtot;
		int iblkbegtot = -1;
		int iblkendtot = -1;
		if (nodes[inode].nodetype != exchangenode) {
			for (iblk=0;iblk<_nblks;iblk++) {
				if (_blks[iblk] == isupbeg) iblkbeg = iblk;
				if (_blks[iblk+1]-1 == isupend) iblkend = iblk;
				if (_blks[iblk] == isupbegtot) iblkbegtot = iblk;
				if (_blks[iblk+1]-1 == isupendtot) iblkendtot = iblk;
			};
			if (iblkbeg == -1 || iblkend == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			if (iblkbegtot == -1 || iblkendtot == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			nodes[inode].indbeg = iblkbeg;
			nodes[inode].indend = iblkend;
			nodes[inode].indbegtot = iblkbegtot;
			nodes[inode].indendtot = iblkendtot;
		} else {
			for (iblk=0;iblk<_nblks;iblk++) {
				if (_blks[iblk] == isupbegtot) iblkbegtot = iblk;
				if (_blks[iblk+1]-1 == isupendtot) iblkendtot = iblk;
			};
			if (iblkbegtot == -1 || iblkendtot == -1) {
				throw " Error in modified tree blocks partitioning";
			};
			nodes[inode].indbegtot = iblkbegtot;
			nodes[inode].indendtot = iblkendtot;
		};
	};

};

// Author: Kharchenko S.A.
// CTree: For each block compute the processor that should compute it (including the Schur complement)
//========================================================================================
void CTree::Block2CpuSchur (int *_bl2cpu) const { // For each block compute the processor that should compute it (including the Schur complement)

//	const char *funcname = "Block2CpuSchur";

// Scan recursively data in the tree

	int inode, iblk, iproc, iblkbeg, iblkend, indtreeloc;
	int inode1, iproc1, iproc2, iblkbeg1, iblkend1;

	for (inode=0;inode<nnodes;inode++) {
		iblkbeg = nodes[inode].indbeg;
		iblkend = nodes[inode].indend;
		iproc = nodes[inode].nodecpu;
		indtreeloc = nodes[inode].indtree;
		if (indtreeloc < 0) {
			for (iblk=iblkbeg;iblk<=iblkend;iblk++) _bl2cpu[iblk] = iproc;
		} else {
			for (inode1=0;inode1<subtreearr[indtreeloc].nnodes;inode1++) {
				iblkbeg1 = subtreearr[indtreeloc].nodes[inode1].indbeg;
				iblkend1 = subtreearr[indtreeloc].nodes[inode1].indend;
				iproc1 = subtreearr[indtreeloc].nodes[inode1].nodecpu;
//				iproc2 = subtreearr[indtreeloc].cpuidarr[iproc1];
				iproc2 = nodes[inode].cpuidarrnode[iproc1];
				for (iblk=iblkbeg1;iblk<=iblkend1;iblk++) _bl2cpu[iblk] = iproc2;
			};
		};
	};

};

// Author: Kharchenko S.A.
// CTree: For each block compute the processor that should compute it
//========================================================================================
void CTree::Block2Cpu (int *_bl2cpu) const { // For each block compute the processor that should compute it

//	const char *funcname = "Block2Cpu";

// Scan data in the tree

	for (int inode=0;inode<nnodes;inode++) {
		int iblkbeg = nodes[inode].indbeg;
		int iblkend = nodes[inode].indend;
		int iproc = nodes[inode].nodecpu;
		if (nodes[inode].nodetype != exchangenode) {
			for (int iblk=iblkbeg;iblk<=iblkend;iblk++) _bl2cpu[iblk] = iproc;
		};
	};

};

// Author: Kharchenko S.A.
// CTree: Check if the tree is sorted
//========================================================================================
bool CTree::IsTreeSorted () const { // Check if the tree is sorted

//	const char *funcname = "IsTreeSorted";

// Scan data in the tree

	bool check = true;

	int index = -1;

	for (int inode=0;inode<nnodes;inode++) {
		int iblkbeg = nodes[inode].indbeg;
		if (nodes[inode].nodetype != exchangenode) {
			if (iblkbeg < index) check = false;
			index = iblkbeg;
		};
	};

	return check;

};

// Author: Kharchenko S.A.
// CTree: Return maximal possible block sparsity
//========================================================================================
CSMatrix CTree::MaximalBlockSparsityForSortedTree (int _nblks) { // Return maximal possible block sparsity

	const char* funcname = "MaximalBlockSparsityForSortedTree";

	bool check = IsTreeSorted ();

	if (!check) {
		CSMatrix amatr = MaximalBlockSparsity (_nblks);
		return amatr;
	};

// Allocate the memory

	int *ordnd, *ordndinv;

	ordnd = new int [nnodes];
	if (!ordnd) MemoryFail (funcname);
	ordndinv = new int [nnodes];
	if (!ordndinv) MemoryFail (funcname);

	int i;

	for (i=0;i<nnodes;i++) ordnd[i] = -1;
	for (i=0;i<nnodes;i++) ordndinv[i] = -1;

// Initially sort nodes of the tree according to the indices

	int nnodes0 = 0;

	for (i=0;i<nnodes;i++) {
		int inodetype = nodes[i].nodetype;
		if (inodetype != exchangenode) {
			ordnd[i] = nnodes0;
			ordndinv[nnodes0] = i;
			nnodes0++;
		};
	};

// Create block sparsity in terms of node numbers 
// (upper triangular part of the matrix by rows)

	int *ianodes;
	int *janodes;

	ianodes = new int [nnodes0+1];
	if (!ianodes) MemoryFail (funcname);
	janodes = new int [nnodes0*nnodes0];
	if (!janodes) MemoryFail (funcname);

	ianodes[0] = 0;
	int nz0, nz = 0;

	for (i=0;i<nnodes0;i++) {

		nz0 = nz;

		janodes[nz] = i;
		nz++;

		int inodecurr = ordndinv[i];

		int inode = inodecurr;

		while (inode >= 0) {

//			int ilev = nodes[inode].nodelv;
			int inodetype = nodes[inode].nodetype;

			int fathernode = nodes[inode].fatherid;

// Go up level

			if (inode == fathernode) {
				inode = -1;
			} else {
				if (inodetype != exchangenode) {
					inode = fathernode;
					int jloc = ordnd[inode];
					if (jloc >= i) {
						janodes[nz] = jloc;
						nz++;
					};
				} else {
					int childnode = -1;
					for (int ichild=0;ichild<nodes[inode].nchilds2;ichild++) {
						int childnodeloc = nodes[inode].childs2[ichild];
						int jloc = ordnd[childnodeloc];
						if (jloc >= i) {
							janodes[nz] = jloc;
							nz++;
						};
						if (ichild == 0) {
							childnode = childnodeloc;
						};
					};
					if (childnode != -1) {
						inode = childnode;
					} else {
						inode = nodes[inode].childs2[0];
					};
				};
			};

		};

		qsort (janodes+nz0,nz-nz0, sizeof(int), compint);

		ianodes[i+1] = nz;

	};

// Symmetrize the sparsity

	CSMatrix anode (nnodes0, nz);

	int *listloc, *ialoc, *jaloc;

	anode.GetSparsity (listloc, ialoc, jaloc);

	for (i=0;i<nnodes0;i++) listloc[i] = 0;
	for (i=0;i<=nnodes0;i++) ialoc[i] = ianodes[i];
	for (i=0;i<nz;i++) jaloc[i] = janodes[i];

	CSMatrix as;

	as = anode.SymmMtr ();

// Create supernodes partitioning

	int *sprndsloc;

	sprndsloc = new int [nnodes0+1];
	if (!sprndsloc) MemoryFail (funcname);

	sprndsloc[0] = 0;

	int inode;

	for (i=0;i<nnodes0;i++) {
		inode = ordndinv[i];
		sprndsloc[i+1] = sprndsloc[i] + (nodes[inode].indend-nodes[inode].indbeg+1);
	};

	if (_nblks != sprndsloc[nnodes0]) {
		throw " CTree::MaximalBlockSparsity: wrong number of blocks ";
	};

// Count the number of elements in the matrix

	as.GetSparsity (listloc, ialoc, jaloc);

	nz = 0;

	int ni, j, jj, nj;

	for (i=0;i<nnodes0;i++) {
		ni = sprndsloc[i+1]-sprndsloc[i];
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			nj = sprndsloc[jj+1]-sprndsloc[jj];
			nz += ni * nj;
		};
	};

	CSMatrix anew (_nblks,nz);

	int *listnew, *ianew, *janew;

	anew.GetSparsity (listnew, ianew, janew);

	for (i=0;i<_nblks;i++) listnew[i] = i;

	nz = 0;
	ianew[0] = 0;

	int irow, icol;

	for (i=0;i<nnodes0;i++) {
		for (irow=sprndsloc[i];irow<sprndsloc[i+1];irow++) {
			for (j=ialoc[i];j<ialoc[i+1];j++) {
				jj = jaloc[j];
				for (icol=sprndsloc[jj];icol<sprndsloc[jj+1];icol++) {
					janew[nz] = icol;
					nz++;
				};
			};
			ianew[irow+1] = nz;
		};
	};

// Free work memory

	delete [] ordnd;
	delete [] ordndinv;
	delete [] ianodes;
	delete [] janodes;
	delete [] sprndsloc;

// Return the result

	return anew;

};

// Author: Kharchenko S.A.
// CTree: Return maximal possible block sparsity
//========================================================================================
CSMatrix CTree::MaximalBlockSparsity (int _nblks) { // Return maximal possible block sparsity

	const char* funcname = "MaximalBlockSparsity";

// Allocate the memory

	int *ordnd, *ordndinv;

	ordnd = new int [nnodes];
	if (!ordnd) MemoryFail (funcname);
	ordndinv = new int [nnodes];
	if (!ordndinv) MemoryFail (funcname);

	CIntInt *srtarr;

	srtarr = new CIntInt [nnodes];
	if (!srtarr) MemoryFail (funcname);

// Initially sort nodes of the tree according to the indices values

	int i;

	int nnodes0 = 0;

	for (i=0;i<nnodes;i++) {
		int inodetype = nodes[i].nodetype;
		if (inodetype != exchangenode) {
			srtarr[nnodes0].intvalue = nodes[i].indbeg;
			srtarr[nnodes0].int2value = i;
			nnodes0++;
		};
	};

	qsort (srtarr, nnodes0, sizeof(CIntInt), CIntInt::CompareIntInt);

	for (i=0;i<nnodes0-1;i++) {
		if (srtarr[i].intvalue == srtarr[i+1].intvalue) {
			throw " CTree::MaximalBlockSparsity: something is wrong with the tree ";
		};
	};

	for (i=0;i<nnodes0;i++) {
		ordndinv[i] = srtarr[i].int2value;
	};
	for (i=nnodes0;i<nnodes;i++) {
		ordndinv[i] = -1;
	};

	for (i=0;i<nnodes;i++) ordnd[i] = -1;

	for (i=0;i<nnodes0;i++) ordnd[ordndinv[i]] = i;

// Create block sparsity in terms of node numbers 
// (upper triangular part of the matrix by rows)

	if (_nblks < nnodes0) {
		throw " CTree::MaximalBlockSparsity: too small number of blocks on entry ";
	};

	int *ianodes;
	int *janodes;

	ianodes = new int [nnodes0+1];
	if (!ianodes) MemoryFail (funcname);
	janodes = new int [nnodes0*nnodes0];
	if (!janodes) MemoryFail (funcname);

	ianodes[0] = 0;
	int nz0, nz = 0;

	for (i=0;i<nnodes0;i++) {

		nz0 = nz;

		janodes[nz] = i;
		nz++;

		int inodecurr = ordndinv[i];

		int inode = inodecurr;

		while (inode >= 0) {

//			int ilev = nodes[inode].nodelv;
			int inodetype = nodes[inode].nodetype;

			int fathernode = nodes[inode].fatherid;

// Go up level

			if (inode == fathernode) {
				inode = -1;
			} else {
				if (inodetype != exchangenode) {
					inode = fathernode;
					int jloc = ordnd[inode];
					if (jloc >= i) {
						janodes[nz] = jloc;
						nz++;
					};
				} else {
					int childnode = -1;
					for (int ichild=0;ichild<nodes[inode].nchilds2;ichild++) {
						int childnodeloc = nodes[inode].childs2[ichild];
						int jloc = ordnd[childnodeloc];
						if (jloc >= i) {
							janodes[nz] = jloc;
							nz++;
						};
						if (ichild == 0) {
							childnode = childnodeloc;
						};
					};
					if (childnode != -1) {
						inode = childnode;
					} else {
						inode = nodes[inode].childs2[0];
					};
				};
			};

		};

		qsort (janodes+nz0,nz-nz0, sizeof(int), compint);

		ianodes[i+1] = nz;

	};

// Symmetrize the sparsity

	CSMatrix anode (nnodes0, nz);

	int *listloc, *ialoc, *jaloc;

	anode.GetSparsity (listloc, ialoc, jaloc);

	for (i=0;i<nnodes0;i++) listloc[i] = 0;
	for (i=0;i<=nnodes0;i++) ialoc[i] = ianodes[i];
	for (i=0;i<nz;i++) jaloc[i] = janodes[i];

	CSMatrix as;

	as = anode.SymmMtr ();

// Create supernodes partitioning

	int *sprndsloc;

	sprndsloc = new int [nnodes0+1];
	if (!sprndsloc) MemoryFail (funcname);

	sprndsloc[0] = 0;

	int inode;

	for (i=0;i<nnodes0;i++) {
		inode = ordndinv[i];
		sprndsloc[i+1] = sprndsloc[i] + (nodes[inode].indend-nodes[inode].indbeg+1);
	};

	if (_nblks != sprndsloc[nnodes0]) {
		throw " CTree::MaximalBlockSparsity: wrong number of blocks ";
	};

// Count the number of elements in the matrix

	as.GetSparsity (listloc, ialoc, jaloc);

	nz = 0;

	int ni, j, jj, nj;

	for (i=0;i<nnodes0;i++) {
		ni = sprndsloc[i+1]-sprndsloc[i];
		for (j=ialoc[i];j<ialoc[i+1];j++) {
			jj = jaloc[j];
			nj = sprndsloc[jj+1]-sprndsloc[jj];
			nz += ni * nj;
		};
	};

	CSMatrix anew (_nblks,nz);

	int *listnew, *ianew, *janew;

	anew.GetSparsity (listnew, ianew, janew);

	for (i=0;i<_nblks;i++) listnew[i] = i;

	nz = 0;
	ianew[0] = 0;

	int irow, icol;

	for (i=0;i<nnodes0;i++) {
		for (irow=sprndsloc[i];irow<sprndsloc[i+1];irow++) {
			for (j=ialoc[i];j<ialoc[i+1];j++) {
				jj = jaloc[j];
				for (icol=sprndsloc[jj];icol<sprndsloc[jj+1];icol++) {
					janew[nz] = icol;
					nz++;
				};
			};
			ianew[irow+1] = nz;
		};
	};

// Free work memory

	delete [] ordnd;
	delete [] ordndinv;
	delete [] srtarr;
	delete [] ianodes;
	delete [] janodes;
	delete [] sprndsloc;

// Return the result

	return anew;

};

// Author: Kharchenko S.A.
// CTree: Pack tree data
//========================================================================================
void CTree::PackTree (int &_length, char *&_obj) { // Pack tree data

	const char* funcname = "PackTree";

// Count the total size of the object

	_length = 10 * sizeof(int) + (2*nproc+nproctot+nnodes+nsubtree+nblkstree*2+1) * sizeof(int) + (2*nproc) * sizeof(double);

	int j, lenloc;
	char *objloc;

	for (j=0;j<nnodes;j++) {
		nodes[j].PackNode (lenloc,objloc);
		_length += lenloc;
		delete [] objloc;
	};
	for (j=0;j<nsubtree;j++) {
		subtreearr[j].PackTree (lenloc,objloc);
		_length += lenloc;
		delete [] objloc;
	};

	ablktree.PackMatrix (lenloc, objloc);

	_length += lenloc;
	delete [] objloc;

	int isizeamtr = lenloc;

// Allocate the structure with the head and init pointers

	_obj = new char [_length];
	if (!_obj) MemoryFail(funcname);

	char* pLoc;

	pLoc = _obj;

	int* pHead;
	int* pcpuidarr;
	int* pcpubeg;
	int* pcpuend;
	int* psize;
	int* psizesubtree;
	int* pblkstree;
	int* pblk2cputree;
	double* pperf;
	double* pmemo;

	pHead = (int *) pLoc;
	pLoc += 10 * sizeof(int);

	pcpuidarr = (int *) pLoc;
	pLoc += nproctot * sizeof(int);

	pcpubeg = (int *) pLoc;
	pLoc += nproc * sizeof(int);

	pcpuend = (int *) pLoc;
	pLoc += nproc * sizeof(int);

	psize = (int *) pLoc;
	pLoc += nnodes * sizeof(int);

	psizesubtree = (int *) pLoc;
	pLoc += nsubtree * sizeof(int);

	pblkstree = (int *) pLoc;
	pLoc += (nblkstree+1) * sizeof(int);

	pblk2cputree = (int *) pLoc;
	pLoc += (nblkstree) * sizeof(int);

	pperf = (double *) pLoc;
	pLoc += nproc * sizeof(double);

	pmemo = (double *) pLoc;
	pLoc += nproc * sizeof(double);

	pHead[0]  = myid;
	pHead[1]  = nlev;
	pHead[2]  = nnodes;
	pHead[3]  = nproc;
	pHead[4]  = nproctot;
	pHead[5]  = rootid;
	pHead[6]  = rootcpu;
	pHead[7]  = nsubtree;
	pHead[8]  = nblkstree;
	pHead[9]  = isizeamtr;

	for(j = 0; j < nproctot; j++) pcpuidarr[j] = cpuidarr[j];
	for(j = 0; j < nproc; j++) pcpubeg[j] = cpuidbeg[j];
	for(j = 0; j < nproc; j++) pcpuend[j] = cpuidend[j];
	for(j = 0; j < nproc; j++) pperf[j] = perf2cpu[j];
	for(j = 0; j < nproc; j++) pmemo[j] = memory2cpu[j];
	for(j = 0; j <= nblkstree; j++) pblkstree[j] = blkstree[j];
	for(j = 0; j <  nblkstree; j++) pblk2cputree[j] = blk2cputree[j];

	int k;

	for (j=0;j<nnodes;j++) {
		nodes[j].PackNode (lenloc,objloc);

		for(k = 0; k < lenloc; k++) pLoc[k] = objloc[k];

		psize[j] = lenloc;

		pLoc += lenloc;
		delete [] objloc;
	};

	for (j=0;j<nsubtree;j++) {
		subtreearr[j].PackTree (lenloc,objloc);

		for(k = 0; k < lenloc; k++) pLoc[k] = objloc[k];

		psizesubtree[j] = lenloc;

		pLoc += lenloc;
		delete [] objloc;
	};

	ablktree.PackMatrix (lenloc,objloc);

	for(k = 0; k < lenloc; k++) pLoc[k] = objloc[k];

	pLoc += lenloc;
	delete [] objloc;

};

// Author: Kharchenko S.A.
// CTree: UnPack tree data
//========================================================================================
void CTree::UnPackTree (int _length, char *_obj) { // UnPack tree data

	const char* funcname = "UnPackTree";

// Free previous arrays

	delete [] nodes;
	delete [] cpuidarr;
	delete [] cpuidbeg;
	delete [] cpuidend;
	delete [] perf2cpu;
	delete [] memory2cpu;
	if (subtreearr != 0) delete [] subtreearr;
	delete [] blkstree;
	delete [] blk2cputree;

// Get head data

	char* pLoc;

	pLoc = _obj;

	int* pHead;
	int* pcpuidarr;
	int* pcpubeg;
	int* pcpuend;
	int* psize;
	int* psizesubtree;
	int* pblkstree;
	int* pblk2cputree;
	double* pperf;
	double* pmemo;
	int isizeablk;

	pHead = (int *) pLoc;
	pLoc += 10 * sizeof(int);

	myid      = pHead[0];
	nlev      = pHead[1];
	nnodes    = pHead[2];
	nproc     = pHead[3];
	nproctot  = pHead[4];
	rootid    = pHead[5];
	rootcpu   = pHead[6];
	nsubtree  = pHead[7];
	nblkstree = pHead[8];
	isizeablk = pHead[9];

// Allocate and get the arrays

	cpuidarr = new int [nproctot];
	if (!cpuidarr) MemoryFail (funcname);
	cpuidbeg = new int [nproc];
	if (!cpuidbeg) MemoryFail(funcname);
	cpuidend = new int [nproc];
	if (!cpuidend) MemoryFail(funcname);
	perf2cpu = new double [nproc];
	if (!perf2cpu) MemoryFail(funcname);
	memory2cpu = new double [nproc];
	if (!memory2cpu) MemoryFail(funcname);
	nodes = new CNode [nnodes];
	if (!nodes) MemoryFail(funcname);
	if (nsubtree > 0) {
		subtreearr = new CTree [nsubtree];
		if (!subtreearr) MemoryFail(funcname);
	};
	blkstree = new int [nblkstree+1];
	if (!blkstree) MemoryFail (funcname);
	blk2cputree = new int [nblkstree];
	if (!blk2cputree) MemoryFail (funcname);

	pcpuidarr = (int *) pLoc;
	pLoc += nproctot * sizeof(int);

	pcpubeg = (int *) pLoc;
	pLoc += nproc * sizeof(int);

	pcpuend = (int *) pLoc;
	pLoc += nproc * sizeof(int);

	psize = (int *) pLoc;
	pLoc += nnodes * sizeof(int);

	psizesubtree = (int *) pLoc;
	pLoc += nsubtree * sizeof(int);

	pblkstree = (int *) pLoc;
	pLoc += (nblkstree+1) * sizeof(int);

	pblk2cputree = (int *) pLoc;
	pLoc += (nblkstree) * sizeof(int);

	pperf = (double *) pLoc;
	pLoc += nproc * sizeof(double);

	pmemo = (double *) pLoc;
	pLoc += nproc * sizeof(double);

	int j;

	for(j = 0; j < nproctot; j++) cpuidarr[j] = pcpuidarr[j];
	for(j = 0; j < nproc; j++) cpuidbeg[j] = pcpubeg[j];
	for(j = 0; j < nproc; j++) cpuidend[j] = pcpuend[j];
	for(j = 0; j < nproc; j++) perf2cpu[j] = pperf[j];
	for(j = 0; j < nproc; j++) memory2cpu[j] = pmemo[j];
	for(j = 0; j <= nblkstree; j++) blkstree[j] = pblkstree[j];
	for(j = 0; j <  nblkstree; j++) blk2cputree[j] = pblk2cputree[j];

	int lenloc;

	for (j=0;j<nnodes;j++) {

		lenloc = psize[j];
		nodes[j].UnPackNode (lenloc,pLoc);
		pLoc += lenloc;

	};

	for (j=0;j<nsubtree;j++) {

		lenloc = psizesubtree[j];
		subtreearr[j].UnPackTree (lenloc,pLoc);
		pLoc += lenloc;

	};

	lenloc = isizeablk;
	ablktree.UnPackMatrix (lenloc,pLoc);
	pLoc += lenloc;

};

// Author: Kharchenko S.A.
// CTree: Compute the list of level nodes for all cpu's
//========================================================================================
void CTree::LevelNodes (int *_lev2nodes) const { // Compute the list of level nodes for all cpu's

//	const char *funcname = "LevelNodes";

	int i, ilevbeg, ilevend;

	for (i=0;i<nproc;i++) {
		LevelNodes (i, ilevbeg, ilevend, _lev2nodes+i*nlev);
	};

};

// Author: Kharchenko S.A.
// CTree: Compute the list of level nodes for given cpu id
//========================================================================================
void CTree::LevelNodes (int _cpuid, int &_ilevbeg, int &_ilevend, int *_lev2node) const { // Compute the list of level nodes for given cpu id

//	const char *funcname = "LevelNodes";

	int i;

	for (i=0;i<nlev;i++) _lev2node[i] = -1;

	int nodeend = cpuidend[_cpuid];
	_ilevend = nodes[nodeend].nodelv;

	int nodebeg = cpuidbeg[_cpuid];
	_ilevbeg = nodes[nodebeg].nodelv;

	int inode = nodeend;
	int ilev, fathernode;
//	int childnode, ichild;

	while (inode >= 0) {

		ilev = nodes[inode].nodelv;
		_lev2node[ilev] = inode;

		fathernode = nodes[inode].fatherid;

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			inode = fathernode;
/*
			childnode = -1;
			for (ichild=0;ichild<nodes[inode].nchilds;ichild++) {
				int childnodeloc = nodes[inode].childs[ichild];
				if (nodes[childnodeloc].nodecpu == _cpuid) {
					childnode = childnodeloc;
				};
			};
			if (childnode != -1) {
				inode = childnode;
			} else {
				inode = nodes[inode].childs[0];
			};
*/
		};

	};

};

// Author: Kharchenko S.A.
// CTree: Compute Schur list for fct
//========================================================================================
void CTree::SchurList (int _cpuid, int &_nlistschur, int *_listschur) const { // Compute Schur list for fct

	const char *funcname = "SchurList";

// Compute lev2node

	int *lev2node;

	lev2node = new int [nlev];
	if (!lev2node) MemoryFail (funcname);

	int ilevbeg, ilevend;

	LevelNodes (_cpuid, ilevbeg, ilevend, lev2node);

// Compute the Schur list

	int ilev, inode, indbegloc, indendloc, j;

	_nlistschur = 0;

	for (ilev=0;ilev<=ilevend;ilev++) {
		inode = lev2node[ilev];
		indbegloc = nodes[inode].indbeg;
		indendloc = nodes[inode].indend;
		for (j=indbegloc;j<=indendloc;j++) {
			_listschur[_nlistschur] = j;
			_nlistschur++;
		};
	};

	qsort (_listschur, _nlistschur, sizeof(int), compint);

// Free work memory

	delete [] lev2node;

};

// Author: Kharchenko S.A.
// CTree: Reassign block numbers according to the partitioning
//========================================================================================
void CTree::ReassignBlockNumbers (int *_blko2n) { // Reassign block numbers according to the partitioning

//	const char *funcname = "ReassignBlockNumbers";

// Main cycle over the nodes

	int inode, ind;

	for (inode=0;inode<nnodes;inode++) {
		ind = nodes[inode].indbeg;
		if (ind >= 0) nodes[inode].indbeg = _blko2n[ind];
		ind = nodes[inode].indend;
		if (ind >= 0) nodes[inode].indend = _blko2n[ind+1]-1;
		ind = nodes[inode].indbegtot;
		if (ind >= 0) nodes[inode].indbegtot = _blko2n[ind];
		ind = nodes[inode].indendtot;
		if (ind >= 0) nodes[inode].indendtot = _blko2n[ind+1]-1;
	};

};

// Author: Kharchenko S.A.
// CTree: Output tree
//========================================================================================
ostream &operator<< (ostream &_stream, const CTree &_tree) { // Output tree

	_stream << " CTree:" << endl;

	_stream << " MyId = " << _tree.myid << endl;
	_stream << " NLev = " << _tree.nlev << " NNodes = " << _tree.nnodes << " NProc = " << _tree.nproc << " NProcTot = " << _tree.nproctot << endl;
	_stream << " RootId = " << _tree.rootid << " RootCpu = " << _tree.rootcpu << endl;
	_stream << " NblksTree = " << _tree.nblkstree << endl;

	OutArr (_stream, " CpuIdArr =", _tree.nproctot, _tree.cpuidarr);
	OutArr (_stream, " CpuIdBeg =", _tree.nproc, _tree.cpuidbeg);
	OutArr (_stream, " CpuIdEnd =", _tree.nproc, _tree.cpuidend);
	OutArr (_stream, " Perf2Cpu =", _tree.nproc, _tree.perf2cpu);
	OutArr (_stream, " Mem2Cpu =", _tree.nproc, _tree.memory2cpu);
	OutArr (_stream, " BlksTree =", _tree.nblkstree+1, _tree.blkstree);
	OutArr (_stream, " Blk2CpuTree =", _tree.nblkstree, _tree.blk2cputree);
	_stream << " ASpBlk = " << _tree.ablktree << endl;

	_stream << " List of nodes:" << endl;

	int i;

	for (i=0;i<_tree.nnodes;i++) {
		_stream << " Inode = " << i << endl;
		_stream << _tree.nodes[i] << endl;
	};

//	_stream << " Comm = " << _tree.comm << endl;

	_stream << " NSubTree = " << _tree.nsubtree << endl;
	_stream << " List of subtrees:" << endl;

	for (i=0;i<_tree.nsubtree;i++) {
		_stream << " Isubtree = " << i << endl;
		_stream << _tree.subtreearr[i] << endl;
	};

	return _stream;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
