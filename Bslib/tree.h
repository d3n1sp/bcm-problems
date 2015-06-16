//------------------------------------------------------------------------------------------------
// File: tree.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include <iostream>

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

// Node: Description of the node of the tree
//
//////////////////////////////////////////////////////////////////////////////

#include "globals.h"
#include "ExchangeMPI.h"
#include "Decomp.h"
#include "smatrix.h"

#ifndef __Node
#define __Node

// Preliminary declarations

enum nodetypes {undefnode=0, compnode, compnode2, exchangenode};

// nodetypes structure specifies the node types:
//
// undefnode     - the type of the node is undefined;
// compnode      - the node performes computations,   data flow from the fathernode of the current node to the childs and back
// compnode2     - the node performes computations,   data flow from the fathernode to the secondary fathernode and back
// exchangenode  - the node performes exchanges only, data flow from the childs to the secondary childs and back

class CTree;

class CNode
{
// Data
	int nodetype;    // nodetype   is the type of the node (compnode - computational node, exchangenode - exchange node)
	int nodeid;      // nodeid     is the ID    number of the current node
	int nodeidext;   // nodeidext  is the ID    number of the current node in the extended tree
	int nodeidsm;    // nodeidsm   is the ID    number of the current node in the included smaller (binary) tree (not collapsed one)
	int nodecpu;     // nodecpu    is the CPU   number of the current node
	int nodelv;      // nodelv     is the level number of the current node
	double nodeperf; // nodeperf   is the total performance of the computing node including the subtree of processors
	int indbeg;      // indbeg     is the starting index that describes data of the node
	int indend;      // indend     is the ending   index that describes data of the node
	int indbegtot;   // indbegtot  is the starting index that describes data of the node including the subtree
	int indendtot;   // indendtot  is the ending   index that describes data of the node including the subtree
	int fatherid;    // fatherid   is the ID  number of the father  node
	int fathercpu;   // fathercpu  is the CPU number of the father  node
	int fatherid2;   // fatherid2  is the ID  number of the exchange father  node
	int fathercpu2;  // fathercpu2 is the CPU number of the exchange father  node
	int nchilds;     // nchilds    is the total number of childs of the current node
	int nchilds2;    // nchilds2   is the total number of exchange childs of the current node
	int *childs;     // childs [nchilds]  array contains the list of child nodes
	int *childs2;    // childs2[nchilds2] array contains the list of exchange child nodes
	int indtree;     // index of the external tree structure (-1 if no external tree structure is available)
	int nprocgl;     // nprocgl is the total number of procs in the tree
	int nprocnode;   // nprocnode is the number of procs in the subtree starting from current node
	int *cpuidarrnode; // cpuidarrnode[nprocnode] is the array that specifies the list of cpu numbers involved in the subtree of current node
	int *cpuidarrnodeloc; // cpuidarrnodeloc[nprocnode] is the array that specifies the local list of cpu numbers involved in the subtree of current node
	int *cpuarrndgl2loc; // cpuidarrnode[nprocgl] is the array that specifies the list of cpu numbers involved in the subtree of current node
	CMPIComm commnode; // commnode structure contains current node communicator
public:
// Functions
// Constructors and destructor
	CNode (); // Zero data constructor
	CNode (int _nchilds); // Memory allocation constructor with zero data
	CNode (int _nchilds, int _nchilds2); // Memory allocation zero data constructor
	CNode (const CNode &_node2); // Copy constructor
	~CNode (); // Destructor
// Operator functions
	CNode &operator= (const CNode &_node2); // Equality operator
// Get/set functions
	int GetNodecpu () const {return nodecpu;};  // Get nodecpu
	int GetNodetype() const {return nodetype;}; // Get nodetype
	int GetNodelv () const {return nodelv;}; // Get nodelv
	int GetIndbeg  () const {return indbeg;};   // Get indbeg
	int IndBeg     () const {return indbeg;};   // Get indbeg
	int GetIndbegtot  () const {return indbegtot;}; // Get indbegtot
	int GetIndend  () const {return indend;};   // Get indend
	int IndEnd     () const {return indend;};   // Get indend
	int GetIndendtot  () const {return indendtot;}; // Get indendtot
	int FatherId   () const {return fatherid;}; // Get fatherid
	int GetNchilds () const {return nchilds;};  // Get nchilds
	int * GetChilds () const {return childs;};  // Get childs
	void SetIndbeg  (int _indbeg) {indbeg = _indbeg;}; // Set indbeg
	void SetIndend  (int _indend) {indend = _indend;}; // Set indend
	void SetIndbegtot (int _indbegtot) {indbegtot = _indbegtot;}; // Set indbegtot
	void SetIndendtot (int _indendtot) {indendtot = _indendtot;}; // Set indendtot
	CMPIComm & GetCommNode () {return commnode;}; // Get commnode
	void SetComm (CMPIComm &_comm) {commnode = _comm;}; // Set commnode
// Node functions
	void OrderIndices (int *_order); // Order node indices inside the node
	void PackNode (int &_length, char *&_obj); // Pack node data
	void UnPackNode (int _length, char *_obj); // UnPack node data
// Input and output
	friend std::ostream &operator<< (std::ostream &_stream, const CNode &_node); // Output node
// Friend classes
	friend class CDecomp;
	friend class CTree;
	friend class CFctDiagR;
	friend class CFctDiagC;
	friend class CMvm;
	friend class CMvmR;
	friend class CMvmC;
	friend class CQrd;
	friend class CQrdC;
	friend class CSMatrix;
	friend class CGSMatrix;
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
	friend class CGSMatrixCS;
};

#endif

// Tree: Description of CPU's tree
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __Tree
#define __Tree

class CSMatrix;

// Preliminary declarations

class CTree
{
// Data
	int myid;           // myid is the processor ID number
	int nlev;           // nlev   is the total number of levels    in the tree
	int nnodes;         // nnodes is the total number of nodes     in the tree
	int nproc;          // nproc  is the current number of processes involved in the tree
	int nproctot;       // nproctot is the total number of processes involved in the tree
	int rootid;         // rootid  is the ID  number of the root node
	int rootcpu;        // rootcpu is the CPU number of the root node
	CNode *nodes;       // nodes[nnodes] is the array of nodes
	int *cpuidarr;      // cpuidarr[nproctot] is the ID number of each local tree processor numbering in terms of up level communicator
	int *cpuidbeg;      // cpuidbeg[nproc] is the ID number of the first (highest level) node in the tree occupied by each CPU
	int *cpuidend;      // cpuidend[nproc] is the ID number of the last  (lowest  level) node in the tree occupied by each CPU
	double *perf2cpu;   // perf2cpu  [nproc] array contains the performances (in arbitrary relative form) of each cpu
	double *memory2cpu; // memory2cpu[nproc] array contains the available memory (in arbitrary relative form) available at each cpu
	int nsubtree;       // nsubtree is the number of elements in the subtree
	CTree *subtreearr;  // subtreearr[nsubtree] is the array of subtrees
	int nblkstree;      // nblkstree is the number of blocks to be stored in a tree
	int *blkstree;      // blkstree[nblkstree+1] is the blocks partitioning of the tree
	int *blk2cputree;   // blk2cputree[nblkstree] is the blocks distribution of the tree
	CSMatrix ablktree;  // ablktree stores the block sparsity to be defined in a tree
	CMPIComm comm;    // comm structure contains current communicator
public:
// Functions
// Constructors and destructor
	CTree (); // Zero length vector constructor
	CTree (int _nproc); // Memory allocation constructor with zero data
	CTree (const CTree &_tree); // Copy constructor
	CTree (int _nproc, int _nchilds); // Create the tree
	CTree (int _nproc, int _nchilds, double *_perf2cpu, double *_memory2cpu); // Create the tree
	~CTree (); // Destructor
// Get/set functions
	int GetMyid () const {return myid;};      // Get myid
	int GetNlev () const {return nlev;};    // Get nlev
	int GetNproc () const {return nproc;};    // Get nproc
	int GetNnodes () const {return nnodes;};  // Get nnodes
	int GetRootcpu () const {return rootcpu;};  // Get rootcpu
	CNode * GetNode () const {return nodes;};  // Get nodes
	int * GetCpuidarr () const {return cpuidarr;};  // Get cpuidarr
	int * GetCpuidend () const {return cpuidend;};  // Get cpuidend
	int GetNsubtree () const {return nsubtree;};  // Get nsubtree
	CTree * GetSubtreearr () const {return subtreearr;};  // Get subtreearr
	int GetNblkstree () const {return nblkstree;};  // Get nblkstree
	int * GetBlkstree () const {return blkstree;};  // Get blkstree
	int * GetBlk2cputree () const {return blk2cputree;};  // Get blk2cputree
	CSMatrix * GetAblktree () {return &ablktree;};  // Get ablktree
	void SetMyid (int _myid) {myid = _myid;}; // Set myid
	void SetNproc (int _nproc) {nproc = _nproc;}; // Set nproc
	void SetNproctot (int _nproctot) {nproctot = _nproctot;}; // Set nproctot
	void SetCpuidarr (int *_cpuidarr) {cpuidarr = _cpuidarr;};  // Set cpuidarr
	void SetNsubtree (int _nsubtree) {nsubtree = _nsubtree;}; // Set nsubtree
	void SetSubtreearr (CTree * _subtreearr) {subtreearr = _subtreearr;}; // Set subtreearr
	void SetNblkstree (int _nblkstree) {nblkstree = _nblkstree;}; // Set nblkstree
	void SetBlkstree (int * _blkstree) {blkstree = _blkstree;}; // Set blkstree
	void SetBlk2cputree (int * _blk2cputree) {blk2cputree = _blk2cputree;}; // Set blk2cputree
	const CMPIComm & GetComm () const {return comm;}; // Get comm
	CMPIComm & GetComm () {return comm;}; // Get comm
	void SetComm (CMPIComm &_comm) {comm = _comm; nproc = comm.GetNproc (); myid = comm.GetMyid ();}; // Set comm
	void SetCommNodes (CMPIComm &_comm); // Set comm for the nodes of the tree
	void SetCommTreeAndNodes (CMPIComm &_comm); // Set comm for the tree and for the nodes of the tree
	void SetCommSchur (CMPIComm &_comm); // Set comm
// Operator functions
	CTree &operator= (const CTree &_tree2); // Equality operator
// Partitioning
	bool IsTreeSorted () const; // Check if the tree is sorted
	CTree ExtendTree   () const; // Extend the tree
	CTree CollapseGeneralTree () const; // Collapse back general extended tree
	CTree CollapseTree () const; // Collapse back extended tree
	CTree TransformTree () const; // Transform the tree for easy tree scan
	void SortNodes (); // Sort nodes of the tree according to the block ordering assumed in the matrix
	void Partition (int &_nblks, int *&_blks) const; // Block partitioning according to the tree
	void PartitionTree (int &_nblks, int *&_blks, int *&_blk2cpu); // Block partitioning according to the tree
	void PartitionTreeSchur (int &_nblks, int *&_blks, int *&_blk2cpu); // Block partitioning according to the tree
	void StorePartitioning (int _nblks, int *_blks, int *_blk2cpu); // Store block partitioning in a tree
	void Block2Cpu (int *_bl2cpu) const; // For each block compute the processor that should compute it
	void Block2CpuSchur (int *_bl2cpu) const; // For each block compute the processor that should compute it (including the Schur complement)
	void TreeInBlocks (int _nsupmx, int &_nblks, int *&_blks); // Modify block partitioning, return modified block partitioning and return tree in terms of blocks
	void TreeInBlocks (int _nblks, int *_blks); // Return tree in terms of blocks
	void TreeInBlocksExt (int _nsupmx, int &_nblks, int *&_blks); // Modify block partitioning introduced on entry, return modified block partitioning and return tree in terms of blocks
	void IncludeBinaryTree (CTree &_treesm); // Include one extended binary tree into smaller one
	void IncludeBinaryTree (int _nproc, CTree &_treesm); // Include one extended binary tree into smaller one
	void IncludeBinaryTreeSchur (int _nproc, CTree &_treesm); // Include one extended binary tree into smaller one
	void IncludeBinaryTreeSchur (CTree &_treesm, int _nblks, int *_blks, CSMatrix &_amatr); // Include one extended binary tree into smaller one
	CSMatrix MaximalBlockSparsity (int _nblks); // Return maximal possible block sparsity
	CSMatrix MaximalBlockSparsityForSortedTree (int _nblks); // Return maximal possible block sparsity
	void LevelNodes (int _cpuid, int &_ilevbeg, int &_ilevend, int *_lev2node) const; // Compute the list of level nodes for given cpu id
	void LevelNodes (int *_lev2nodes) const; // Compute the list of level nodes for all cpu's
	void SchurList (int _cpuid, int &_nlistschur, int *_listschur) const; // Compute Schur list for fct
	void CpuListsForNodesSchur (); // Compute the list of processors for each node of the tree including Schur subtree's
	void CpuListsForNodes (); // Compute the list of processors for each node of the tree
	void ModifyCpuLists (int *_cpuarr); // Modify the proc numbers in the nodes lists
	void ReassignBlockNumbers (int *_blko2n); // Reassign block numbers according to the partitioning
// Input and output
	void PackTree (int &_length, char *&_obj); // Pack tree data
	void UnPackTree (int _length, char *_obj); // UnPack tree data
	friend std::ostream &operator<< (std::ostream &_stream, const CTree &_tree); // Output tree
// Friend classes
	friend class CDecomp;
	friend class CFctDiagR;
	friend class CFctDiagC;
	friend class CMvm;
	friend class CMvmR;
	friend class CMvmC;
	friend class CQrd;
	friend class CQrdC;
	friend class CSMatrix;
	friend class CGSMatrix;
	friend class CGSMatrixR;
	friend class CGSMatrixRS;
	friend class CGSMatrixCS;
	friend class CSVector;
	friend class CSVectorC;
};

#endif

