//------------------------------------------------------------------------------------------------
// File: CpuTree.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include <iostream>

#include "tree.h"

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

// CpuTree: Description of the CPU's tree
//
//////////////////////////////////////////////////////////////////////////////

#ifdef __FV_2004__
namespace flowvision {
#endif

#ifndef __CpuTree
#define __CpuTree

#ifdef __FV_2004__
class EqnSolver::CNode;
class EqnSolver::CTree;
#else
class CNode;
class CTree;
#endif

class CCpuTree
{
// Data
#ifdef __FV_2004__
	EqnSolver::CTree tree;    // tree structure describes the distribution of the data over processors
	EqnSolver::CTree treeext; // treeext structure describes the extended distribution of the data over processors
	EqnSolver::CTree treehc;  // treehc structure describes the hypercells tree according to the nested dessection ordering
#else
	CTree tree;    // tree structure describes the distribution of the data over processors
	CTree treeext; // treeext structure describes the extended distribution of the data over processors
	CTree treehc;  // treehc structure describes the hypercells tree according to the nested dessection ordering
#endif
// Functions
public:
// Set/get functions
	void SetNproc (int _nproc) {tree.SetNproc (_nproc); treeext.SetNproc (_nproc);} // Set nproc
	void SetComm (CMPIComm &_comm) {tree.SetComm (_comm); treeext.SetComm (_comm);} // Set communicator
	void SetCommSchur (CMPIComm &_comm) {tree.SetCommSchur (_comm); treeext.SetCommSchur (_comm);} // Set communicator
	void SetCommNodes (CMPIComm &_comm) {tree.SetCommNodes (_comm); treeext.SetCommNodes (_comm);} // Set communicator
	void SortNodes () {treeext.SortNodes (); tree = treeext.CollapseGeneralTree (); treehc.SortNodes ();}; // Sort tree nodes
#ifdef __FV_2004__
	void SetTree (EqnSolver::CTree &_tree) {tree = _tree;}; // Set tree
	void SetTreeExt (EqnSolver::CTree &_treeext) {treeext = _treeext;}; // Set treeext
	void SetTreeHc (EqnSolver::CTree &_treehc) {treehc = _treehc;}; // Set treehc
	EqnSolver::CNode *Nodes   () const {return tree.GetNode ();}; // Get nodes
	EqnSolver::CTree * GetTree () { return &tree;}; // Get tree
	EqnSolver::CTree * GetTreeext () { return &treeext;}; // Get treeext
	EqnSolver::CTree * GetTreeHc () { return &treehc;}; // Get treehc
#else
	void SetTree (CTree &_tree) {tree = _tree;}; // Set tree
	void SetTreeExt (CTree &_treeext) {treeext = _treeext;}; // Set treeext
	void SetTreeHc (CTree &_treehc) {treehc = _treehc;}; // Set treehc
	CNode * Nodes   () const {return tree.GetNode ();};  // Get nodes
	CTree * GetTree () { return &tree;}; // Get tree
	CTree * GetTreeext () { return &treeext;}; // Get treeext
	CTree * GetTreeHc () { return &treehc;}; // Get treehc
#endif
	int * CpuIdEnd () const {return tree.GetCpuidend ();};  // Get cpuidend
	friend std::ostream &operator<< (std::ostream &_stream, const CCpuTree &_tree) { // Output tree
		_stream << " CCpuTree:" << std::endl;
		_stream << " Tree = " << _tree.tree << std::endl;
		_stream << " TreeExt = " << _tree.treeext << std::endl;
		_stream << " TreeHc = " << _tree.treehc << std::endl;
		return _stream;
	};
// Friend classes
#ifdef __FV_2004__
	friend class EqnSolver::CDecomp;
#else
	friend class CDecomp;
#endif
};

#endif

#ifdef __FV_2004__
} // namespace flowvision
#endif
