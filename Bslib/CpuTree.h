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

#ifndef __CpuTree
#define __CpuTree

class CNode;
class CTree;

class CCpuTree
{
// Data
	CTree tree;    // tree structure describes the distribution of the data over processors
	CTree treeext; // treeext structure describes the extended distribution of the data over processors
	CTree treehc;  // treehc structure describes the hypercells tree according to the nested dessection ordering
// Functions
public:
// Set/get functions
	void SetNproc (int _nproc) {tree.SetNproc (_nproc); treeext.SetNproc (_nproc);} // Set nproc
	void SetComm (CMPIComm &_comm) {tree.SetComm (_comm); treeext.SetComm (_comm);} // Set communicator
	void SetCommSchur (CMPIComm &_comm) {tree.SetCommSchur (_comm); treeext.SetCommSchur (_comm);} // Set communicator
	void SetCommNodes (CMPIComm &_comm) {tree.SetCommNodes (_comm); treeext.SetCommNodes (_comm);} // Set communicator
	void SortNodes () {treeext.SortNodes (); tree = treeext.CollapseGeneralTree (); treehc.SortNodes ();}; // Sort tree nodes
	void SetTree (CTree &_tree) {tree = _tree;}; // Set tree
	void SetTreeExt (CTree &_treeext) {treeext = _treeext;}; // Set treeext
	void SetTreeHc (CTree &_treehc) {treehc = _treehc;}; // Set treehc
	CNode * Nodes   () const {return tree.GetNode ();};  // Get nodes
	CTree * GetTree () { return &tree;}; // Get tree
	CTree * GetTreeext () { return &treeext;}; // Get treeext
	CTree * GetTreeHc () { return &treehc;}; // Get treehc
	int * CpuIdEnd () const {return tree.GetCpuidend ();};  // Get cpuidend
	friend std::ostream &operator<< (std::ostream &_stream, const CCpuTree &_tree) { // Output tree
		_stream << " CCpuTree:" << std::endl;
		_stream << " Tree = " << _tree.tree << std::endl;
		_stream << " TreeExt = " << _tree.treeext << std::endl;
		_stream << " TreeHc = " << _tree.treehc << std::endl;
		return _stream;
	};
// Friend classes
	friend class CDecomp;
};

#endif
