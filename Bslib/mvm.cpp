//------------------------------------------------------------------------------------------------
// File: mvm.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#ifdef __PMPITRACE__
#include "mpi.h"
#include "mpe.h"
#endif

#include <iostream>
#include <fstream>
#include <cstring>

#include "tree.h"
#include "gsmatrix.h"
#include "smatrix.h"
#include "svector.h"
#include "globals.h"
#include "mvm.h"
#include "ExchangeMPI.h"

using namespace std;

// Author: Kharchenko S.A.
// CMvm: Constructor
//========================================================================================
CMvm::CMvm (int _nblks, int _nsups, int _nlev) { // Constructor

	const char *funcname = "CMvm_00";

// Allocate the data

	nlev  = _nlev;
	nsups = _nsups;
	nrhs  = 1;
	nloc  = 0;
	nlocext = 0;
	int nblks  = _nblks;

	int i;

	bsx = new int [nsups];
	if (!bsx) MemoryFail (funcname);
	for (i=0;i<nsups;i++) bsx[i] = -1;
	bspx = new int [nsups];
	if (!bspx) MemoryFail (funcname);
	for (i=0;i<nsups;i++) bspx[i] = -1;
	bsx2 = new int [nblks];
	if (!bsx2) MemoryFail (funcname);
	for (i=0;i<nblks;i++) bsx2[i] = -1;
	bspx2 = new int [nblks];
	if (!bspx2) MemoryFail (funcname);
	for (i=0;i<nblks;i++) bspx2[i] = -1;
	lev2node = new int [nlev];
	if (!lev2node) MemoryFail (funcname);
	for (i=0;i<nlev;i++) lev2node[i] = -1;
	blk2cpu = new int [nblks];
	if (!blk2cpu) MemoryFail (funcname);
	nchildsarr = new int [nlev];
	if (!nchildsarr) MemoryFail (funcname);
	nlistupsndarr = new int [nlev];
	if (!nlistupsndarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistupsndarr[i] = 0;
	nlistuprcvarr = new int [nlev];
	if (!nlistuprcvarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistuprcvarr[i] = 0;
	nlistdnsndarr = new int [nlev];
	if (!nlistdnsndarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistdnsndarr[i] = 0;
	nlistdnrcvarr = new int [nlev];
	if (!nlistdnrcvarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistdnrcvarr[i] = 0;
	nlistfctarr = new int [nlev];
	if (!nlistfctarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistfctarr[i] = 0;
	listupsndarr = new int * [nlev];
	if (!listupsndarr) MemoryFail (funcname);
	listuprcvarr = new int * [nlev];
	if (!listuprcvarr) MemoryFail (funcname);
	listdnsndarr = new int * [nlev];
	if (!listdnsndarr) MemoryFail (funcname);
	listdnrcvarr = new int * [nlev];
	if (!listdnrcvarr) MemoryFail (funcname);
	list2upsndarr = new int * [nlev];
	if (!list2upsndarr) MemoryFail (funcname);
	list2uprcvarr = new int * [nlev];
	if (!list2uprcvarr) MemoryFail (funcname);
	list2dnsndarr = new int * [nlev];
	if (!list2dnsndarr) MemoryFail (funcname);
	list2dnrcvarr = new int * [nlev];
	if (!list2dnrcvarr) MemoryFail (funcname);
	listfctarr = new int * [nlev];
	if (!listfctarr) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// CMvm: Constructor
//========================================================================================
CMvm::CMvm (bool _is2index, const CTree &_tree, const CGSMatrix &_gmtra) { // Constructor

	const char *funcname = "CMvm_01";

// Allocate the data

	int nblksloc = _gmtra.GetNblksr ();

	nlev  = _tree.nlev;
	nsups = _gmtra.nsupr;
	nrhs  = 0;

	int i;

	if (!_is2index) {
		bsx = new int [nsups];
		if (!bsx) MemoryFail (funcname);
		for (i=0;i<nsups;i++) bsx[i] = -1;
		bspx = new int [nsups];
		if (!bspx) MemoryFail (funcname);
		for (i=0;i<nsups;i++) bspx[i] = -1;
	} else {
		bsx = new int [nblksloc];
		if (!bsx) MemoryFail (funcname);
		for (i=0;i<nblksloc;i++) bsx[i] = -1;
		bspx = new int [nblksloc];
		if (!bspx) MemoryFail (funcname);
		for (i=0;i<nblksloc;i++) bspx[i] = -1;
	};
	bsx2 = new int [nblksloc];
	if (!bsx2) MemoryFail (funcname);
	for (i=0;i<nblksloc;i++) bsx2[i] = -1;
	bspx2 = new int [nblksloc];
	if (!bspx2) MemoryFail (funcname);
	for (i=0;i<nblksloc;i++) bspx2[i] = -1;
	lev2node = new int [nlev];
	if (!lev2node) MemoryFail (funcname);
	for (i=0;i<nlev;i++) lev2node[i] = -1;
	blk2cpu = new int [_gmtra.nblksr];
	if (!blk2cpu) MemoryFail (funcname);
	nchildsarr = new int [nlev];
	if (!nchildsarr) MemoryFail (funcname);
	nlistupsndarr = new int [nlev];
	if (!nlistupsndarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistupsndarr[i] = 0;
	nlistuprcvarr = new int [nlev];
	if (!nlistuprcvarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistuprcvarr[i] = 0;
	nlistdnsndarr = new int [nlev];
	if (!nlistdnsndarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistdnsndarr[i] = 0;
	nlistdnrcvarr = new int [nlev];
	if (!nlistdnrcvarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistdnrcvarr[i] = 0;
	nlistfctarr = new int [nlev];
	if (!nlistfctarr) MemoryFail (funcname);
	for (i=0;i<nlev;i++) nlistfctarr[i] = 0;
	listupsndarr = new int * [nlev];
	if (!listupsndarr) MemoryFail (funcname);
	listuprcvarr = new int * [nlev];
	if (!listuprcvarr) MemoryFail (funcname);
	listdnsndarr = new int * [nlev];
	if (!listdnsndarr) MemoryFail (funcname);
	listdnrcvarr = new int * [nlev];
	if (!listdnrcvarr) MemoryFail (funcname);
	list2upsndarr = new int * [nlev];
	if (!list2upsndarr) MemoryFail (funcname);
	list2uprcvarr = new int * [nlev];
	if (!list2uprcvarr) MemoryFail (funcname);
	list2dnsndarr = new int * [nlev];
	if (!list2dnsndarr) MemoryFail (funcname);
	list2dnrcvarr = new int * [nlev];
	if (!list2dnrcvarr) MemoryFail (funcname);
	listfctarr = new int * [nlev];
	if (!listfctarr) MemoryFail (funcname);

	for (i=0;i<nlev;i++) {
		list2upsndarr[i] = new int [0];
		if (!list2upsndarr[i]) MemoryFail (funcname);
		list2uprcvarr[i] = new int [0];
		if (!list2uprcvarr[i]) MemoryFail (funcname);
		list2dnsndarr[i] = new int [0];
		if (!list2dnsndarr[i]) MemoryFail (funcname);
		list2dnrcvarr[i] = new int [0];
		if (!list2dnrcvarr[i]) MemoryFail (funcname);
	};

// Compute blk2cpu array

	_tree.Block2Cpu (blk2cpu);

// Init nloc, bsx and bspx

	nloc = 0;

	int ilist;

	if (!_is2index) {
		for (ilist=0;ilist<_gmtra.nlist;ilist++) {
			int iblk = _gmtra.listb[ilist];
			for (int isup=_gmtra.blksr[iblk];isup<_gmtra.blksr[iblk+1];isup++) {
				int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
				bsx[isup] = nloc;
				bspx[isup] = nloc;
				nloc += blai;
			};
		};
	} else {
		for (ilist=0;ilist<_gmtra.nlist;ilist++) {
			int iblk = _gmtra.listb[ilist];
			bsx[iblk] = nloc;
			bspx[iblk] = nloc;
			nloc += _gmtra.bl2ndr[iblk+1]-_gmtra.bl2ndr[iblk];
		};
	};

// Create the list of nodes

	int *nodeslist;

	nodeslist = new int [nlev];
	if (!nodeslist) MemoryFail (funcname);

	for (i=0;i<nlev;i++) nodeslist[i] = -1;

// Search the tree

	int inode;

	int nodeend = _tree.cpuidend[_tree.myid];

	inode = nodeend;

	while (inode >= 0) {

		int ilev = _tree.nodes[inode].nodelv;
		int inodetype = _tree.nodes[inode].nodetype;

		lev2node[ilev] = inode;

		nodeslist[ilev] = inode;
		nchildsarr[ilev] = _tree.nodes[inode].nchilds;

		int fathernode = _tree.nodes[inode].fatherid;

// Go up level

		if (inode == fathernode) {
			inode = -1;
		} else {
			if (inodetype != exchangenode) {
				inode = fathernode;
			} else {
				int childnode = -1;
				for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {
					int childnodeloc = _tree.nodes[inode].childs2[ichild];
					if (_tree.nodes[childnodeloc].nodecpu == _tree.myid) {
						childnode = childnodeloc;
					};
				};
				if (childnode != -1) {
					inode = childnode;
				} else {
					inode = _tree.nodes[inode].childs2[0];
				};
			};
		};

	};

// Create the lists of send/receive blocks

	int *imask;
	int *listblk;

	imask = new int [_gmtra.nblksr];
	if (!imask) MemoryFail (funcname);
	listblk = new int [_gmtra.nblksr];
	if (!listblk) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<_gmtra.nblksr;i++) imask[i] = icycle;

	int nlistloc;
	int *listloc;

// Down from upper to the lower levels

// Send down lists

	int ilev;
	for (ilev=0;ilev<nlev;ilev++) {
		icycle++;
		int inode = nodeslist[ilev];
		nlistloc = 0;
		if (inode >= 0) {

			int inodetype = _tree.nodes[inode].nodetype;
			int nodecpuloc = _tree.nodes[inode].nodecpu;
			int iblkbeg = _tree.nodes[inode].indbeg;
			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
			int nchildsloc = _tree.nodes[inode].nchilds;

			if (nodecpuloc != _tree.myid) goto NextSendDn;

// Send comp to comp/comp2 

			if (inodetype == compnode) {

// Next node if this is the degenerate comp node

				if (nchildsloc <= 1) goto NextSendDn;

// Combine all previous data into the send list

				int fathernode = _tree.nodes[inode].fatherid;
				int fathertype = _tree.nodes[fathernode].nodetype;

				for (int jlev=0;jlev<=ilev;jlev++) {
					int jnode = nodeslist[jlev];
					int jnodetype  = _tree.nodes[jnode].nodetype;
					int jblkbeg;
					int jblkend;
					if (jlev < ilev && fathertype == exchangenode) {
						jblkbeg = _tree.nodes[jnode].indbegtot;
						jblkend = _tree.nodes[jnode].indendtot;
					} else {
						jblkbeg = _tree.nodes[jnode].indbeg;
						jblkend = _tree.nodes[jnode].indend;
					};
					if (jnodetype != exchangenode) {
						for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
							if (imask[jblk] != icycle) {
								listblk[nlistloc] = jblk;
								nlistloc++;
								imask[jblk] = icycle;
							};
						};
					};
				};
			};

// Send comp2 to exchange

			if (inodetype == compnode2) {

// Add data into the send list

				int fathernode2 = _tree.nodes[inode].fatherid2;

				if (_tree.nodes[fathernode2].nodecpu == _tree.myid) goto NextSendDn;

				for (int jblk=iblkbeg;jblk<=iblkend;jblk++) {
					if (imask[jblk] != icycle) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			};

// Send exchange to comp

			if (inodetype == exchangenode) {

				int fathernode = _tree.nodes[inode].fatherid;

				int jblkbeg = _tree.nodes[fathernode].indbeg;
				int jblkend = _tree.nodes[fathernode].indend;
				int jblk;
				for (jblk=jblkbeg;jblk<=jblkend;jblk++) {
					if (imask[jblk] != icycle) {
						imask[jblk] = icycle;
					};
				};

// Add data into the send list

				for (jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
					if (imask[jblk] != icycle) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			};

		};

NextSendDn:;

		qsort (listblk, nlistloc, sizeof(int), compint);

		listloc = new int [nlistloc];
		if (!listloc) MemoryFail (funcname);
		for (i=0;i<nlistloc;i++) {
			listloc[i] = listblk[i];
		};

		nlistdnsndarr[ilev] = nlistloc;
		listdnsndarr[ilev] = listloc;

	};

// Receive down lists

	for (ilev=0;ilev<nlev;ilev++) {
		icycle++;
		int inode = nodeslist[ilev];
		nlistloc = 0;
		if (inode >= 0) {

			int inodetype = _tree.nodes[inode].nodetype;
			int nodecpu = _tree.nodes[inode].nodecpu;
			int fathercpu = _tree.nodes[inode].fathercpu;
			int fathernode = _tree.nodes[inode].fatherid;
			int fathertype = _tree.nodes[fathernode].nodetype;
//			int iblkbeg = _tree.nodes[inode].indbeg;
//			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nchildsloc = _tree.nodes[inode].nchilds;

			if (nodecpu != _tree.myid) goto NextRecvDn;

			if ((inodetype == compnode && fathertype == compnode) || inodetype == compnode2) {

				if (fathercpu == _tree.myid) goto NextRecvDn;

// Combine all previous data into the receive list

				if (inodetype == compnode) {
					for (int jlev=0;jlev<=ilev-1;jlev++) {
						int jnode = nodeslist[jlev];
						int jnodetype  = _tree.nodes[jnode].nodetype;
						int jblkbeg;
						int jblkend;
						jblkbeg = _tree.nodes[jnode].indbeg;
						jblkend = _tree.nodes[jnode].indend;
						if (jnodetype != exchangenode && _tree.nodes[jnode].nodecpu != _tree.myid) {
							for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
								if (imask[jblk] != icycle) {
									listblk[nlistloc] = jblk;
									nlistloc++;
									imask[jblk] = icycle;
								};
							};
						};
					};
				} else {
					for (int jlev=0;jlev<=ilev-1;jlev++) {
						int jnode = nodeslist[jlev];
						int jnodetype  = _tree.nodes[jnode].nodetype;
						int jblkbeg;
						int jblkend;
						if (jlev != ilev-1) {
							jblkbeg = _tree.nodes[jnode].indbegtot;
							jblkend = _tree.nodes[jnode].indendtot;
						} else {
							jblkbeg = _tree.nodes[jnode].indbeg;
							jblkend = _tree.nodes[jnode].indend;
						};
						if (jnodetype != exchangenode && _tree.nodes[jnode].nodecpu != _tree.myid) {
							for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
								if (imask[jblk] != icycle) {
									listblk[nlistloc] = jblk;
									nlistloc++;
									imask[jblk] = icycle;
								};
							};
						};
					};
				};

			} else if (inodetype == compnode && fathertype == exchangenode) {

				if (fathercpu == _tree.myid) goto NextRecvDn;

				int fathernode3 = _tree.nodes[fathernode].fatherid;

				int jblkbegtot = _tree.nodes[fathernode3].indbeg;
				int jblkendtot = _tree.nodes[fathernode3].indend;

				int jblk;
				for (jblk=jblkbegtot;jblk<=jblkendtot;jblk++) {
					imask[jblk] = icycle;
				};

				jblkbegtot = _tree.nodes[fathernode].indbegtot;
				jblkendtot = _tree.nodes[fathernode].indendtot;

				for (jblk=jblkbegtot;jblk<=jblkendtot;jblk++) {
					if (imask[jblk] != icycle && blk2cpu[jblk] != _tree.myid) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			} else if (inodetype == exchangenode) {

				for (int jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
					if (imask[jblk] != icycle && blk2cpu[jblk] != _tree.myid) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			};

		};

NextRecvDn:;

		qsort (listblk, nlistloc, sizeof(int), compint);

		listloc = new int [nlistloc];
		if (!listloc) MemoryFail (funcname);
		for (i=0;i<nlistloc;i++) {
			listloc[i] = listblk[i];
		};

		nlistdnrcvarr[ilev] = nlistloc;
		listdnrcvarr[ilev] = listloc;

	};

// Up from the lower to an upper levels

// Send up lists

	for (ilev=0;ilev<nlev;ilev++) {
		icycle++;
		int inode = nodeslist[ilev];
		nlistloc = 0;
		if (inode >= 0) {

			int inodetype = _tree.nodes[inode].nodetype;
			int nodecpu = _tree.nodes[inode].nodecpu;
			int fathernode = _tree.nodes[inode].fatherid;
			int fathercpu = _tree.nodes[inode].fathercpu;
			int fathertype = _tree.nodes[fathernode].nodetype;
//			int iblkbeg = _tree.nodes[inode].indbeg;
//			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
//			int nchildsloc = _tree.nodes[inode].nchilds;

			if (nodecpu != _tree.myid) goto NextSendUp;

			if ((inodetype == compnode && fathertype == compnode) || inodetype == compnode2) {

// Next node if this is the degenerate comp node

				if (fathercpu == _tree.myid) goto NextSendUp;

// Combine all previous data into the send list

				for (int jlev=0;jlev<=ilev-1;jlev++) {
					int jnode = nodeslist[jlev];
					int jnodetype  = _tree.nodes[jnode].nodetype;
					int jblkbeg;
					int jblkend;
					if (jlev < ilev-1 && inodetype == compnode2) {
						jblkbeg = _tree.nodes[jnode].indbegtot;
						jblkend = _tree.nodes[jnode].indendtot;
					} else {
						jblkbeg = _tree.nodes[jnode].indbeg;
						jblkend = _tree.nodes[jnode].indend;
					};
					if (jnodetype != exchangenode) {
						for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
							if (imask[jblk] != icycle) {
								listblk[nlistloc] = jblk;
								nlistloc++;
								imask[jblk] = icycle;
							};
						};
					};
				};

			} else if (inodetype == compnode && fathertype == exchangenode) {

				if (fathercpu == _tree.myid) goto NextSendUp;

				int fathernode2 = _tree.nodes[fathernode].fatherid;

				int jblkbegtot = _tree.nodes[fathernode2].indbeg;
				int jblkendtot = _tree.nodes[fathernode2].indend;

				int jblk;

				for (jblk=jblkbegtot;jblk<=jblkendtot;jblk++) {
					if (imask[jblk] != icycle) {
						imask[jblk] = icycle;
					};
				};

				jblkbegtot = _tree.nodes[fathernode].indbegtot;
				jblkendtot = _tree.nodes[fathernode].indendtot;

				for (jblk=jblkbegtot;jblk<=jblkendtot;jblk++) {
					if (imask[jblk] != icycle) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			} else if (inodetype == exchangenode) {

				for (int jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
					if (imask[jblk] != icycle && blk2cpu[jblk] != _tree.myid) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			};

		};

NextSendUp:;

		qsort (listblk, nlistloc, sizeof(int), compint);

		listloc = new int [nlistloc];
		if (!listloc) MemoryFail (funcname);
		for (i=0;i<nlistloc;i++) {
			listloc[i] = listblk[i];
		};

		nlistupsndarr[ilev] = nlistloc;
		listupsndarr[ilev] = listloc;

	};

// Receive up lists

	for (ilev=0;ilev<nlev;ilev++) {
		icycle++;
		int inode = nodeslist[ilev];
		nlistloc = 0;
		if (inode >= 0) {

			int inodetype = _tree.nodes[inode].nodetype;
			int nodecpu = _tree.nodes[inode].nodecpu;
			int fathernode = _tree.nodes[inode].fatherid;
			int fathercpu = _tree.nodes[inode].fathercpu;
			int fathertype = _tree.nodes[fathernode].nodetype;
//			int iblkbeg = _tree.nodes[inode].indbeg;
//			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
			int nchildsloc = _tree.nodes[inode].nchilds;

			if (nodecpu != _tree.myid) goto NextRecvUp;

			if (inodetype == compnode) {

// Next node if this is the degenerate node

				if (nchildsloc <= 1) goto NextRecvUp;

// Combine all previous data into the receive list

				for (int jlev=0;jlev<=ilev;jlev++) {
					int jnode = nodeslist[jlev];
					int jnodetype  = _tree.nodes[jnode].nodetype;
					int jblkbeg;
					int jblkend;
					if (jlev < ilev && fathertype == exchangenode) {
						jblkbeg = _tree.nodes[jnode].indbegtot;
						jblkend = _tree.nodes[jnode].indendtot;
					} else {
						jblkbeg = _tree.nodes[jnode].indbeg;
						jblkend = _tree.nodes[jnode].indend;
					};
					if (jnodetype != exchangenode) {
						for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
							if (imask[jblk] != icycle) {
								listblk[nlistloc] = jblk;
								nlistloc++;
								imask[jblk] = icycle;
							};
						};
					};
				};

			} else if (inodetype == compnode2) {

				if (fathercpu == _tree.myid) goto NextRecvUp;

				for (int jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
					if (imask[jblk] != icycle && blk2cpu[jblk] == _tree.myid) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			} else if (inodetype == exchangenode) {

				int jblkbegtot = _tree.nodes[fathernode].indbeg;
				int jblkendtot = _tree.nodes[fathernode].indend;

				int jblk;

				for (jblk=jblkbegtot;jblk<=jblkendtot;jblk++) {
					if (imask[jblk] != icycle) {
						imask[jblk] = icycle;
					};
				};

				for (jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
					if (imask[jblk] != icycle) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};

			};
		};

NextRecvUp:;

		qsort (listblk, nlistloc, sizeof(int), compint);

		listloc = new int [nlistloc];
		if (!listloc) MemoryFail (funcname);
		for (i=0;i<nlistloc;i++) {
			listloc[i] = listblk[i];
		};

		nlistuprcvarr[ilev] = nlistloc;
		listuprcvarr[ilev] = listloc;

	};

// Determine the tree type

	int treetype = 0;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = nodeslist[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
			if (inodetype != compnode) treetype = 1;
		};
	};

// Fct lists

	for (ilev=0;ilev<nlev;ilev++) {
		icycle++;
		int inode = nodeslist[ilev];
		nlistloc = 0;
		if (inode >= 0) {

			int inodetype = _tree.nodes[inode].nodetype;
//			int nodecpu = _tree.nodes[inode].nodecpu;
//			int fathernode = _tree.nodes[inode].fatherid;
//			int fathercpu = _tree.nodes[inode].fathercpu;
//			int fathertype = _tree.nodes[fathernode].nodetype;
//			int iblkbeg = _tree.nodes[inode].indbeg;
//			int iblkend = _tree.nodes[inode].indend;
//			int iblkbegtot = _tree.nodes[inode].indbegtot;
//			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
//			int nchildsloc = _tree.nodes[inode].nchilds;

			if (true) {

// Combine all previous data into the list

				for (int jlev=0;jlev<=ilev-1;jlev++) {
					int jnode = nodeslist[jlev];
					int jnodetype  = _tree.nodes[jnode].nodetype;
					int jblkbeg;
					int jblkend;
					if (treetype == 1) {
						if (inodetype == compnode2 && jlev == ilev-1) {
							jblkbeg = _tree.nodes[jnode].indbeg;
							jblkend = _tree.nodes[jnode].indend;
						} else {
							jblkbeg = _tree.nodes[jnode].indbegtot;
							jblkend = _tree.nodes[jnode].indendtot;
						};
					} else {
						jblkbeg = _tree.nodes[jnode].indbeg;
						jblkend = _tree.nodes[jnode].indend;
					};
					if (jnodetype != exchangenode) {
						for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
							if (imask[jblk] != icycle) {
								listblk[nlistloc] = jblk;
								nlistloc++;
								imask[jblk] = icycle;
							};
						};
					};
				};
			};
		};

		qsort (listblk, nlistloc, sizeof(int), compint);

		listloc = new int [nlistloc];
		if (!listloc) MemoryFail (funcname);
		for (i=0;i<nlistloc;i++) {
			listloc[i] = listblk[i];
		};

		nlistfctarr[ilev] = nlistloc;
		listfctarr[ilev] = listloc;

	};

// Create the list of blocks to be added

	icycle++;

	for (ilist=0;ilist<_gmtra.nlist;ilist++) {
		int iblk = _gmtra.listb[ilist];
		imask[iblk] = icycle;
	};

	nlistloc = 0;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = nodeslist[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
//			int nodecpu = _tree.nodes[inode].nodecpu;
			int iblkbeg = _tree.nodes[inode].indbeg;
			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
			int nchildsloc = _tree.nodes[inode].nchilds;

			if (inodetype != exchangenode) {

// Next node if this is the degenerate node

				if (nchildsloc <= 1) goto NextNode;

// Combine all previous data into the receive list

				int jblkbeg;
				int jblkend;
				if (treetype == 1) {
					jblkbeg = iblkbegtot;
					jblkend = iblkendtot;
				} else {
					jblkbeg = iblkbeg;
					jblkend = iblkend;
				};

				for (int jblk=jblkbeg;jblk<=jblkend;jblk++) {
					if (imask[jblk] != icycle) {
						listblk[nlistloc] = jblk;
						nlistloc++;
						imask[jblk] = icycle;
					};
				};
			};
		};

NextNode:;

	};

	qsort (listblk, nlistloc, sizeof(int), compint);

// Compute nlocext and bspx

	nlocext = nloc;

	if (!_is2index) {

		for (ilist=0;ilist<nlistloc;ilist++) {
			int iblk = listblk[ilist];
			for (int isup=_gmtra.blksr[iblk];isup<_gmtra.blksr[iblk+1];isup++) {
				int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
				bspx[isup] = nlocext;
				nlocext += blai;
			};
		};

	} else {

		for (ilist=0;ilist<nlistloc;ilist++) {
			int iblk = listblk[ilist];
			bspx[iblk] = nlocext;
			nlocext += _gmtra.bl2ndr[iblk+1]-_gmtra.bl2ndr[iblk];
		};

	};

// Free work arrays

	delete [] nodeslist;
	delete [] imask;
	delete [] listblk;

};

//========================================================================
// CMvm: Destructor
//========================================================================================
CMvm::~CMvm () { // Destructor
//	std::cout << " On entry to CMvm destructor " << std::endl;
	nsups   = 0;
	nloc    = 0;
	nlocext = 0;
	nrhs    = 0;
	delete [] bsx;
	delete [] bspx;
	delete [] bsx2;
	delete [] bspx2;
	delete [] lev2node;
	delete [] blk2cpu;
	delete [] nchildsarr;
	delete [] nlistupsndarr;
	delete [] nlistuprcvarr;
	delete [] nlistdnsndarr;
	delete [] nlistdnrcvarr;
	delete [] nlistfctarr;
	for (int ilev=0;ilev<nlev;ilev++) {
		delete [] listupsndarr [ilev];
		delete [] listuprcvarr [ilev];
		delete [] listdnsndarr [ilev];
		delete [] listdnrcvarr [ilev];
		delete [] list2upsndarr[ilev];
		delete [] list2uprcvarr[ilev];
		delete [] list2dnsndarr[ilev];
		delete [] list2dnrcvarr[ilev];
		delete [] listfctarr   [ilev];
	};
	delete [] listupsndarr;
	delete [] listuprcvarr;
	delete [] listdnsndarr;
	delete [] listdnrcvarr;
	delete [] list2upsndarr;
	delete [] list2uprcvarr;
	delete [] list2dnsndarr;
	delete [] list2dnrcvarr;
	delete [] listfctarr;
	nlev = 0;
//	std::cout << " On return from CMvm destructor " << std::endl;
};

// Author: Kharchenko S.A.
// CMvmR: Constructor
//========================================================================================
CMvmR::CMvmR (const CTree &_tree, const CGSMatrixRS &_gmtra, int _nrhs): CMvm (false, _tree, _gmtra) { // Constructor

	const char *funcname = "CMvmR_00";

// Allocate the data

	CSVector temp(nlocext,1,_nrhs);

	px = temp;
	pax = temp;

	CSVector vectdummy;

	temp = vectdummy;

// Allocate and init work arrays

	int *imask;
	int *listblk;

	imask = new int [_gmtra.nblksr];
	if (!imask) MemoryFail (funcname);
	listblk = new int [_gmtra.nblksr];
	if (!listblk) MemoryFail (funcname);

	int icycle = -1;

	for (int i=0;i<_gmtra.nblksr;i++) imask[i] = icycle;

// Determine the tree type

	int treetype = 0;

	int ilev;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = lev2node[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
			if (inodetype != compnode) treetype = 1;
		};
	};

// Create the list of blocks to be added

	icycle++;

	int nlistloc = 0;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = lev2node[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
//			int nodecpu = _tree.nodes[inode].nodecpu;
			int iblkbeg = _tree.nodes[inode].indbeg;
			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
			int nchildsloc = _tree.nodes[inode].nchilds;

			if (inodetype != exchangenode) {

// Next node if this is the degenerate node

				if (nchildsloc <= 1) goto NextNode;

// Combine all previous data into the receive list

				if (treetype == 0) {
					for (int jblk=iblkbeg;jblk<=iblkend;jblk++) {
						if (imask[jblk] != icycle) {
							listblk[nlistloc] = jblk;
							nlistloc++;
							imask[jblk] = icycle;
						};
					};
				} else {
					for (int jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
						if (imask[jblk] != icycle) {
							listblk[nlistloc] = jblk;
							nlistloc++;
							imask[jblk] = icycle;
						};
					};
				};
			};
		};

NextNode:;

	};

	qsort (listblk, nlistloc, sizeof(int), compint);

// Compute nlocext and bspx

	int nadd = 0;

	for (int ilist=0;ilist<nlistloc;ilist++) {
		int iblk = listblk[ilist];
		for (int isup=_gmtra.blksr[iblk];isup<_gmtra.blksr[iblk+1];isup++) {
			int blai = _gmtra.sprndr[isup+1]-_gmtra.sprndr[isup];
			nadd += blai;
		};
	};

	vectsndrcv = new double [_nrhs*nadd];
	if (!vectsndrcv) MemoryFail (funcname);

	nrhs = _nrhs;

	delete [] imask;
	delete [] listblk;

};

// Author: Kharchenko S.A.
// CMvmR: Constructor
//========================================================================================
CMvmR::CMvmR (const CTree &_tree, const CGSMatrixR &_gmtra, int _nrhs): CMvm (true, _tree, _gmtra) { // Constructor

	const char *funcname = "CMvmR_00";

// Allocate the data

	CSVector temp(nlocext,1,_nrhs);

	px = temp;
	pax = temp;

	CSVector vectdummy;

	temp = vectdummy;

// Allocate and init work arrays

	int *imask;
	int *listblk;

	imask = new int [_gmtra.nblksr];
	if (!imask) MemoryFail (funcname);
	listblk = new int [_gmtra.nblksr];
	if (!listblk) MemoryFail (funcname);

	int icycle = -1;

	for (int i=0;i<_gmtra.nblksr;i++) imask[i] = icycle;

// Determine the tree type

	int treetype = 0;

	int ilev;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = lev2node[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
			if (inodetype != compnode) treetype = 1;
		};
	};

// Create the list of blocks to be added

	icycle++;

	int nlistloc = 0;

	for (ilev=0;ilev<nlev;ilev++) {
		int inode = lev2node[ilev];
		if (inode >= 0) {
			int inodetype = _tree.nodes[inode].nodetype;
//			int nodecpu = _tree.nodes[inode].nodecpu;
			int iblkbeg = _tree.nodes[inode].indbeg;
			int iblkend = _tree.nodes[inode].indend;
			int iblkbegtot = _tree.nodes[inode].indbegtot;
			int iblkendtot = _tree.nodes[inode].indendtot;
//			int nblksloc = iblkend-iblkbeg+1;
			int nchildsloc = _tree.nodes[inode].nchilds;

			if (inodetype != exchangenode) {

// Next node if this is the degenerate node

				if (nchildsloc <= 1) goto NextNode;

// Combine all previous data into the receive list

				if (treetype == 0) {
					for (int jblk=iblkbeg;jblk<=iblkend;jblk++) {
						if (imask[jblk] != icycle) {
							listblk[nlistloc] = jblk;
							nlistloc++;
							imask[jblk] = icycle;
						};
					};
				} else {
					for (int jblk=iblkbegtot;jblk<=iblkendtot;jblk++) {
						if (imask[jblk] != icycle) {
							listblk[nlistloc] = jblk;
							nlistloc++;
							imask[jblk] = icycle;
						};
					};
				};
			};
		};

NextNode:;

	};

	qsort (listblk, nlistloc, sizeof(int), compint);

// Compute nlocext and bspx

	int nadd = 0;

	for (int ilist=0;ilist<nlistloc;ilist++) {
		int iblk = listblk[ilist];
		nadd += (_gmtra.blksr[iblk+1]-_gmtra.blksr[iblk]);
	};

	vectsndrcv = new double [_nrhs*nadd];
	if (!vectsndrcv) MemoryFail (funcname);

	nrhs = _nrhs;

	delete [] imask;
	delete [] listblk;

};

// Author: Kharchenko S.A.
// CMvmR: Constructor
//========================================================================================
CMvmR::CMvmR (bool _bsecond, const CTree &_tree, const CGSMatrixRS &_gmtra, int _nrhs): CMvm (_gmtra.nblksr, _gmtra.nsupr, _tree.nproc) { // Constructor

	const char *funcname = "CMvmR_01";

// Compute blk2cpu array

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

//	_tree.Block2Cpu (blk2cpu);
	int nblkstree = _tree.GetNblkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();

	int i;

	for (i=0;i<nblkstree;i++) blk2cpu[i] = pblk2cputree[i];

// Init the mask and sp2blk arrays

	int nsuploc = _gmtra.nsupr;
	int nblksloc = _gmtra.nblksr;

	int *sp2blk;
	int *listr;
	int *listc;
	int *imaskr;
	int *imaskc;

	sp2blk = new int [nsuploc];
	if (!sp2blk) MemoryFail (funcname);
	listr = new int [nsuploc];
	if (!listr) MemoryFail (funcname);
	listc = new int [nsuploc];
	if (!listc) MemoryFail (funcname);
	imaskr = new int [nsuploc];
	if (!imaskr) MemoryFail (funcname);
	imaskc = new int [nsuploc];
	if (!imaskc) MemoryFail (funcname);

	int j;

	int *pblksr, *pblksc;

	_gmtra.GetBlks (pblksr, pblksc);

	for (i=0;i<nblksloc;i++) {
		for (j=pblksr[i];j<pblksr[i+1];j++) {
			sp2blk[j] = i;
		};
	};

	int icycle = -1;
	for (i=0;i<nsuploc;i++) imaskr[i] = icycle;
	for (i=0;i<nsuploc;i++) imaskc[i] = icycle;

// Create the local lists of rows and columns by scanning local block rows

	CSMatrixRS *pmtrarr = _gmtra.GetMtrarr ();

	int nlistr = 0;
	int nlistc = 0;

	icycle++;

	int iblk, ilist, jj;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			for (ilist=0;ilist<pmtrarr[iblk].nlist;ilist++) {
				listr[nlistr] = pmtrarr[iblk].list[ilist];
				nlistr++;
				for (j=pmtrarr[iblk].ia[ilist];j<pmtrarr[iblk].ia[ilist+1];j++) {
					jj = pmtrarr[iblk].ja[j];
					if (imaskc[jj] != icycle) {
						listc[nlistc] = jj;
						nlistc++;
						imaskc[jj] = icycle;
					};
				};
			};
		};
	};

	qsort (listr, nlistr, sizeof(int), compint);
	qsort (listc, nlistc, sizeof(int), compint);

//	cout << " Nlistr = " << nlistr << " Nlistc = " << nlistc << endl;

// Create local data to be exchanged

	int *listsnd;
	int *iasnd;
	int *iptr;

	listsnd = new int [nsuploc];
	if (!listsnd) MemoryFail (funcname);
	iasnd = new int [nproc+1];
	if (!iasnd) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iasnd[i] = 0;

	int iproc;

	for (i=0;i<nlistc;i++) {
		jj = listc[i];
		iblk = sp2blk[jj];
		iproc = blk2cpu[iblk];
		iasnd[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iasnd[i+1] = iasnd[i] + iasnd[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iasnd[i];

	int k;

	for (i=0;i<nlistc;i++) {
		jj = listc[i];
		iblk = sp2blk[jj];
		iproc = blk2cpu[iblk];
		k = iptr[iproc];
		listsnd[k] = jj;
		iptr[iproc]++;
	};

// Exchange the data

	int NObjSend = nproc-1;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ObjTypeSend[iobj] = 1;
			ObjIDSend[iobj] = myid;
			CpuIDSend[iobj] = iproc;
			ObjSizeSend[iobj] = (iasnd[iproc+1]-iasnd[iproc]) * sizeof (int);
			ObjSend[iobj] = (char *) (listsnd+iasnd[iproc]);
			iobj++;
		};
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info = CMPIExchange::DataExchangeMPI (pcomm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
	if(info) throw " CMvmR::CMvmR: Error in DataExchangeMPI";

// Free send arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Create the array of local numbering of global indices

	int *gl2loc;
	int *listgl;

	gl2loc = new int [nsuploc];
	if(!gl2loc) MemoryFail(funcname);
	listgl = new int [nsuploc];
	if(!listgl) MemoryFail(funcname);

	int ind=0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			for (i=pblksr[iblk];i<pblksr[iblk+1];i++) {
				gl2loc[i] = ind;
				listgl[ind] = i;
				ind++;
			};
		} else {
			for (i=pblksr[iblk];i<pblksr[iblk+1];i++) {
				gl2loc[i] = -1;
			};
		};
	};

	int nlistgl = ind;

//	cout << " Nlistgl = " << nlistgl << endl;

// Create combined list of indices to be used and to be sent to other cpu's

	int *imask;
	int *listtot;

	imask = new int [nsuploc];
	if(!imask) MemoryFail(funcname);
	listtot = new int [nsuploc];
	if(!listtot) MemoryFail(funcname);

	for (i=0;i<nsuploc;i++) imask[i] = -1;

	int nlist = 0;

	for (i=0;i<nlistr;i++) {
		jj = listr[i];
		listtot[nlist] = jj;
		nlist++;
		imask[jj] = 1;
	};

	for (i=0;i<nlistc;i++) {
		jj = listc[i];
		iblk = sp2blk[jj];
		iproc = blk2cpu[iblk];
		if (iproc == myid) {
			if (imask[jj] < 0) {
				listtot[nlist] = jj;
				nlist++;
				imask[jj] = 1;
			};
		};
	};

//	cout << " Nlist = " << nlist << endl;

	int *iarr;

	int ni, kii;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		ni = ObjSizeRecv[iobj] / sizeof (int);
		iarr = (int *) ObjRecv[iobj];
		for (j=0;j<ni;j++) {
			jj = iarr[j];
			if (imask[jj] < 0) {
				listtot[nlist] = jj;
				nlist++;
				imask[jj] = 1;
			};
		};
	};

	qsort (listtot, nlist, sizeof(int), compint);

	for (i=0;i<nlist;i++) {
		jj = listtot[i];
		imask[jj] = i;
	};

	nloc = nlist;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			for (j=iasnd[iproc];j<iasnd[iproc+1];j++) {
				jj = listsnd[j];
				if (imask[jj] == -1) {
					listtot[nlist] = jj;
					imask[jj] = nlist;
					nlist++;
				};
			};
		};
	};

	nlocext = nlist;

//	cout << " Nloc = " << nloc << " nlocext = " << endl;

// Create send and receive lists (up and down)

	iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistupsndarr[myid] = nloc;
			listupsndarr[myid] = new int [nloc];
			if(!listupsndarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nloc;kii++) {
				jj = listtot[kii];
				listupsndarr[myid][kii] = jj;
			};
		} else {
			ni = ObjSizeRecv[iobj] / sizeof (int);
			iarr = (int *) ObjRecv[iobj];
			nlistupsndarr[iproc] = ni;
			listupsndarr[iproc] = new int [ni];
			if(!listupsndarr[iproc]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii];
				listupsndarr[iproc][kii] = jj;
			};
			iobj++;
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistuprcvarr[myid] = nlistr;
			listuprcvarr[myid] = new int [nlistr];
			if(!listuprcvarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nlistr;kii++) {
				jj = listr[kii];
				listuprcvarr[myid][kii] = jj;
			};
		} else {
			ni = iasnd[iproc+1]-iasnd[iproc];
			iarr = listsnd+iasnd[iproc];
			nlistuprcvarr[iproc] = ni;
			listuprcvarr[iproc] = new int [ni];
			if(!listuprcvarr[iproc]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii];
				listuprcvarr[iproc][kii] = jj;
			};
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistdnsndarr[myid] = nlistr;
			listdnsndarr[myid] = new int [nlistr];
			if(!listdnsndarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nlistr;kii++) {
				jj = listr[kii];
				listdnsndarr[myid][kii] = jj;
			};
		} else {
			ni = iasnd[iproc+1]-iasnd[iproc];
			iarr = listsnd+iasnd[iproc];
			nlistdnsndarr[iproc] = ni;
			listdnsndarr[iproc] = new int [ni];
			if(!listdnsndarr[iproc]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii];
				listdnsndarr[iproc][kii] = jj;
			};
		};
	};

	iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistdnrcvarr[myid] = nloc;
			listdnrcvarr[myid] = new int [nloc];
			if(!listdnrcvarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nloc;kii++) {
				jj = listtot[kii];
				listdnrcvarr[myid][kii] = jj;
			};
		} else {
			ni = ObjSizeRecv[iobj] / sizeof (int);
			iarr = (int *) ObjRecv[iobj];
			nlistdnrcvarr[iproc] = ni;
			listdnrcvarr[iproc] = new int [ni];
			if(!listdnrcvarr[iproc]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii];
				listdnrcvarr[iproc][kii] = jj;
			};
			iobj++;
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		ni = 0;
		nlistfctarr[iproc] = ni;
		listfctarr[iproc] = new int [ni];
		if(!listfctarr[iproc]) MemoryFail(funcname);
	};

// Compute nlocext, bsx and bspx

	int *psprndr, *psprndc;

	_gmtra.GetSups (psprndr, psprndc);

	int nlocpt, nlocextpt = 0;

	int blai;

	for (ilist=0;ilist<nloc;ilist++) {
		jj = listtot[ilist];
		blai = psprndr[jj+1]-psprndr[jj];
		bspx[jj] = nlocextpt;
		nlocextpt += blai;
	};

	nlocpt = nlocextpt;

	for (ilist=nloc;ilist<nlocext;ilist++) {
		jj = listtot[ilist];
		blai = psprndr[jj+1]-psprndr[jj];
		bspx[jj] = nlocextpt;
		nlocextpt += blai;
	};

	nloc = nlocpt;
	nlocext = nlocextpt;

	nlocpt = 0;

	for (ilist=0;ilist<nlistgl;ilist++) {
		jj = listgl[ilist];
		blai = psprndr[jj+1]-psprndr[jj];
		bsx[jj] = nlocpt;
		nlocpt += blai;
	};

// Allocate the vector data

	CSVector temp(nlocext,1,_nrhs);

	px = temp;
	pax = temp;

	CSVector vectdummy;

	temp = vectdummy;

	int nadd = nlocext;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ni = iasnd[iproc+1]-iasnd[iproc];
			iarr = listsnd+iasnd[iproc];
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii];
				blai = psprndr[jj+1]-psprndr[jj];
				nadd += blai;
			};
		};
	};

	int nadd1 = nlocext;

	iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ni = ObjSizeRecv[iobj] / sizeof (int);
			iarr = (int *) ObjRecv[iobj];
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii];
				blai = psprndr[jj+1]-psprndr[jj];
				nadd1 += blai;
			};
			iobj++;
		};
	};

	if (nadd1 > nadd) nadd = nadd1;

	vectsndrcv = new double [nadd];
	if (!vectsndrcv) MemoryFail (funcname);

//	cout << " Nadd = " << nadd << endl;

// Free receive arrays

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Free work arrays

	delete [] sp2blk;
	delete [] listr;
	delete [] listc;
	delete [] imaskr;
	delete [] imaskc;
	delete [] listsnd;
	delete [] iasnd;
	delete [] iptr;
	delete [] gl2loc;
	delete [] listgl;
	delete [] imask;
	delete [] listtot;

};

// Author: Kharchenko S.A.
// CMvmR: Constructor
//========================================================================================
CMvmR::CMvmR (bool _bsecond2ind, const CTree &_tree, const CGSMatrixR &_gmtra, int _nrhs): CMvm (_gmtra.nblksr, 0, _tree.nproc) { // Constructor

	const char *funcname = "CMvmR_02";

// Compute blk2cpu array

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

//	_tree.Block2Cpu (blk2cpu);
	int nblkstree = _tree.GetNblkstree ();
	int *pblk2cputree = _tree.GetBlk2cputree ();

	int i;

	for (i=0;i<nblkstree;i++) blk2cpu[i] = pblk2cputree[i];

// Init the mask and sp2blk arrays

	int nblksloc = _gmtra.nblksr;
	int *pblks = _gmtra.blksr;

	int *listrblk;
	int *listcblk;
	int *imaskrblk;
	int *imaskcblk;

	listrblk = new int [nblksloc];
	if (!listrblk) MemoryFail (funcname);
	listcblk = new int [nblksloc];
	if (!listcblk) MemoryFail (funcname);
	imaskrblk = new int [nblksloc];
	if (!imaskrblk) MemoryFail (funcname);
	imaskcblk = new int [nblksloc];
	if (!imaskcblk) MemoryFail (funcname);

	int j;

	int icycle = -1;

	for (i=0;i<nblksloc;i++) imaskrblk[i] = icycle;
	for (i=0;i<nblksloc;i++) imaskcblk[i] = icycle;

// Create the local lists of rows and columns by scanning local block rows

	CSMatrixR *pmtrarr = _gmtra.GetMtrarr ();

	int nlistrblk = 0;
	int nlistcblk = 0;

	icycle++;

	int iblk, ilist, jj, jblk;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			for (ilist=0;ilist<pmtrarr[iblk].nlist;ilist++) {
				jblk = pmtrarr[iblk].list2[ilist];
				if (imaskrblk[jblk] != icycle) {
					listrblk[nlistrblk] = jblk;
					nlistrblk++;
					imaskrblk[jblk] = icycle;
				};
				for (j=pmtrarr[iblk].ia[ilist];j<pmtrarr[iblk].ia[ilist+1];j++) {
					jblk = pmtrarr[iblk].ja2[j];
					if (imaskcblk[jblk] != icycle) {
						listcblk[nlistcblk] = jblk;
						nlistcblk++;
						imaskcblk[jblk] = icycle;
					};
				};
			};
		};
	};

	qsort (listrblk, nlistrblk, sizeof(int), compint);
	qsort (listcblk, nlistcblk, sizeof(int), compint);

// Create segmented mask arrays and create complete two indices rows and columns lists

	int *ibsr, *ibsc;

	ibsr = new int [nblksloc];
	if (!ibsr) MemoryFail (funcname);
	ibsc = new int [nblksloc];
	if (!ibsc) MemoryFail (funcname);

	for (i=0;i<nblksloc;i++) ibsr[i] = -1;
	for (i=0;i<nblksloc;i++) ibsc[i] = -1;

	int nzrmax = 0;

	for (i=0;i<nlistrblk;i++) {
		iblk = listrblk[i];
		ibsr[iblk] = nzrmax;
		nzrmax += pblks[iblk+1]-pblks[iblk];
	};

	int *imaskr, *imaskc;

	imaskr = new int [nzrmax];
	if (!imaskr) MemoryFail (funcname);

	for (i=0;i<nzrmax;i++) imaskr[i] = -1;

	int nzcmax = 0;

	for (i=0;i<nlistcblk;i++) {
		iblk = listcblk[i];
		ibsc[iblk] = nzcmax;
		nzcmax += pblks[iblk+1]-pblks[iblk];
	};

	imaskc = new int [nzcmax];
	if (!imaskc) MemoryFail (funcname);

	for (i=0;i<nzcmax;i++) imaskc[i] = -1;

	int *listr, *listc, *listr2, *listc2;

	listr = new int [nzrmax];
	if (!listr) MemoryFail (funcname);
	listc = new int [nzcmax];
	if (!listc) MemoryFail (funcname);
	listr2 = new int [nzrmax];
	if (!listr2) MemoryFail (funcname);
	listc2 = new int [nzcmax];
	if (!listc2) MemoryFail (funcname);

	int nlistr = 0;
	int nlistc = 0;

	int ibs;

	icycle++;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			for (ilist=0;ilist<pmtrarr[iblk].nlist;ilist++) {
				jj = pmtrarr[iblk].list[ilist];
				jblk = pmtrarr[iblk].list2[ilist];
				ibs = ibsr[jblk]+jj;
				if (imaskr[ibs] != icycle) {
					listr[nlistr] = jj;
					listr2[nlistr] = jblk;
					nlistr++;
					imaskr[ibs] = icycle;
				};
				for (j=pmtrarr[iblk].ia[ilist];j<pmtrarr[iblk].ia[ilist+1];j++) {
					jj = pmtrarr[iblk].ja[j];
					jblk = pmtrarr[iblk].ja2[j];
					ibs = ibsc[jblk]+jj;
					if (imaskc[ibs] != icycle) {
						listc[nlistc] = jj;
						listc2[nlistc] = jblk;
						nlistc++;
						imaskc[ibs] = icycle;
					};
				};
			};
		};
	};

// Sort two index lists

	int nlistmax = nlistr;
	if (nlistc > nlistmax) nlistmax = nlistc;

	CInd2Int *i2arr;

	i2arr = new CInd2Int [nlistmax];
	if (!i2arr) MemoryFail (funcname);

	for (i=0;i<nlistr;i++) {
		i2arr[i].indx = listr2[i];
		i2arr[i].indy = listr[i];
		i2arr[i].intvalue = i;
	};

	qsort (i2arr, nlistr, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

	for (i=0;i<nlistr;i++) {
		listr2[i] = i2arr[i].indx;
		listr[i] = i2arr[i].indy;
	};

	for (i=0;i<nlistc;i++) {
		i2arr[i].indx = listc2[i];
		i2arr[i].indy = listc[i];
		i2arr[i].intvalue = i;
	};

	qsort (i2arr, nlistc, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

	for (i=0;i<nlistc;i++) {
		listc2[i] = i2arr[i].indx;
		listc[i] = i2arr[i].indy;
	};

	delete [] i2arr;

//	cout << " Nlistr = " << nlistr << " Nlistc = " << nlistc << endl;

// Create local data to be exchanged

	int *listsnd;
	int *iasnd;
	int *iptr;

	listsnd = new int [nlistmax*2];
	if (!listsnd) MemoryFail (funcname);
	iasnd = new int [nproc+1];
	if (!iasnd) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) iasnd[i] = 0;

	int iproc;

	for (i=0;i<nlistc;i++) {
		jj = listc[i];
		iblk = listc2[i];
		iproc = blk2cpu[iblk];
		iasnd[iproc+1]++;
	};

	for (i=0;i<nproc;i++) iasnd[i+1] = iasnd[i] + iasnd[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iasnd[i];

	int k;

	for (i=0;i<nlistc;i++) {
		jj = listc[i];
		iblk = listc2[i];
		iproc = blk2cpu[iblk];
		k = iptr[iproc];
		listsnd[k*2] = jj;
		listsnd[k*2+1] = iblk;
		iptr[iproc]++;
	};

// Exchange the data

	int NObjSend = nproc-1;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ObjTypeSend[iobj] = 1;
			ObjIDSend[iobj] = myid;
			CpuIDSend[iobj] = iproc;
			ObjSizeSend[iobj] = 2 * (iasnd[iproc+1]-iasnd[iproc]) * sizeof (int);
			ObjSend[iobj] = (char *) (listsnd+2*iasnd[iproc]);
			iobj++;
		};
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info = CMPIExchange::DataExchangeMPI (pcomm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
	if(info) throw " CMvmR::CMvmR: Error in DataExchangeMPI";

// Free send arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Create the array of local numbering of global indices

	int *imaskblk;
	int *listblk;

	imaskblk = new int [nblksloc];
	if (!imaskblk) MemoryFail(funcname);
	listblk = new int [nblksloc];
	if (!listblk) MemoryFail(funcname);

	for (iblk=0;iblk<nblksloc;iblk++) imaskblk[iblk] = -1;

	int nlistblk = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (blk2cpu[iblk] == myid) {
			listblk[nlistblk] = iblk;
			nlistblk++;
			imaskblk[iblk] = 1;
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			for (j=iasnd[iproc];j<iasnd[iproc+1];j++) {
				jj = listsnd[j*2];
				jblk = listsnd[j*2+1];
				if (imaskblk[jblk] < 0) {
					listblk[nlistblk] = jblk;
					nlistblk++;
					imaskblk[jblk] = 1;
				};
			};
		};
	};

	qsort (listblk, nlistblk, sizeof(int), compint);

	int *ibsgl;

	ibsgl = new int [nblksloc];
	if(!ibsgl) MemoryFail(funcname);

	for (iblk=0;iblk<nblksloc;iblk++) ibsgl[iblk] = -1;

	int nz = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		ibsgl[iblk] = nz;
		if (imaskblk[iblk] >= 0) {
			nz += pblks[iblk+1]-pblks[iblk];
		};
	};

//	cout << " Nlistgl = " << nlistgl << endl;

// Create combined list of indices to be used and to be sent to other cpu's

	int *imask;
	int *listtot;
	int *list2tot;

	imask = new int [nz];
	if(!imask) MemoryFail(funcname);
	listtot = new int [nz];
	if(!listtot) MemoryFail(funcname);
	list2tot = new int [nz];
	if(!list2tot) MemoryFail(funcname);

	for (i=0;i<nz;i++) imask[i] = -1;

	int nlist = 0;

	for (i=0;i<nlistr;i++) {
		jj = listr[i];
		jblk = listr2[i];
		iproc = blk2cpu[jblk];
		ibs = ibsgl[jblk]+jj;
		if (iproc == myid) {
			if (imask[ibs] < 0) {
				listtot[nlist] = jj;
				list2tot[nlist] = jblk;
				nlist++;
				imask[ibs] = 1;
			};
		};
	};

	for (i=0;i<nlistc;i++) {
		jj = listc[i];
		jblk = listc2[i];
		iproc = blk2cpu[jblk];
		ibs = ibsgl[jblk]+jj;
		if (iproc == myid) {
			if (imask[ibs] < 0) {
				listtot[nlist] = jj;
				list2tot[nlist] = jblk;
				nlist++;
				imask[ibs] = 1;
			};
		};
	};

//	cout << " Nlist = " << nlist << endl;

	int *iarr;

	int ni, kii;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		ni = ObjSizeRecv[iobj] / (sizeof (int)*2);
		iarr = (int *) ObjRecv[iobj];
		for (j=0;j<ni;j++) {
			jj = iarr[j*2];
			jblk = iarr[j*2+1];
			iproc = blk2cpu[jblk];
			ibs = ibsgl[jblk]+jj;
			if (iproc == myid) {
				if (imask[ibs] < 0) {
					listtot[nlist] = jj;
					list2tot[nlist] = jblk;
					nlist++;
					imask[ibs] = 1;
				};
			};
		};
	};

	i2arr = new CInd2Int [nlist];
	if (!i2arr) MemoryFail (funcname);

	for (i=0;i<nlist;i++) {
		i2arr[i].indx = list2tot[i];
		i2arr[i].indy = listtot[i];
		i2arr[i].intvalue = i;
	};

	qsort (i2arr, nlist, sizeof(CInd2Int), CInd2Int::CompareInd2Int);

	for (i=0;i<nlist;i++) {
		list2tot[i] = i2arr[i].indx;
		listtot[i] = i2arr[i].indy;
	};

	delete [] i2arr;

	for (i=0;i<nlist;i++) {
		jj = listtot[i];
		jblk = list2tot[i];
		ibs = ibsgl[jblk]+jj;
		imask[ibs] = i;
	};

	nloc = nlist;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			for (j=iasnd[iproc];j<iasnd[iproc+1];j++) {
				jj = listsnd[j*2];
				jblk = listsnd[j*2+1];
				ibs = ibsgl[jblk]+jj;
				if (imask[ibs] < 0) {
					listtot[nlist] = jj;
					list2tot[nlist] = jblk;
					imask[ibs] = nlist;
					nlist++;
				};
			};
		};
	};

	nlocext = nlist;

//	cout << " Nloc = " << nloc << " nlocext = " << endl;

// Create send and receive lists (up and down)

	int iproc1;

	iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistupsndarr[myid] = nloc;
			listupsndarr[myid] = new int [nloc];
			if(!listupsndarr[myid]) MemoryFail(funcname);
			list2upsndarr[myid] = new int [nloc];
			if(!list2upsndarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nloc;kii++) {
				jj = listtot[kii];
				jblk = list2tot[kii];
				listupsndarr[myid][kii] = jj;
				list2upsndarr[myid][kii] = jblk;
			};
		} else {
			ni = ObjSizeRecv[iobj] / (2*sizeof (int));
			iarr = (int *) ObjRecv[iobj];
			iproc1 = CpuIDRecv[iobj];
			nlistupsndarr[iproc1] = ni;
			listupsndarr[iproc1] = new int [ni];
			if(!listupsndarr[iproc1]) MemoryFail(funcname);
			list2upsndarr[iproc1] = new int [ni];
			if(!list2upsndarr[iproc1]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii*2];
				jblk = iarr[kii*2+1];
				listupsndarr[iproc1][kii] = jj;
				list2upsndarr[iproc1][kii] = jblk;
			};
			iobj++;
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistuprcvarr[myid] = nlistr;
			listuprcvarr[myid] = new int [nlistr];
			if(!listuprcvarr[myid]) MemoryFail(funcname);
			list2uprcvarr[myid] = new int [nlistr];
			if(!list2uprcvarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nlistr;kii++) {
				jj = listr[kii];
				jblk = listr2[kii];
				listuprcvarr[myid][kii] = jj;
				list2uprcvarr[myid][kii] = jblk;
			};
		} else {
			ni = iasnd[iproc+1]-iasnd[iproc];
			iarr = listsnd+iasnd[iproc]*2;
			nlistuprcvarr[iproc] = ni;
			listuprcvarr[iproc] = new int [ni];
			if(!listuprcvarr[iproc]) MemoryFail(funcname);
			list2uprcvarr[iproc] = new int [ni];
			if(!list2uprcvarr[iproc]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii*2];
				jblk = iarr[kii*2+1];
				listuprcvarr[iproc][kii] = jj;
				list2uprcvarr[iproc][kii] = jblk;
			};
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistdnsndarr[myid] = nlistr;
			listdnsndarr[myid] = new int [nlistr];
			if(!listdnsndarr[myid]) MemoryFail(funcname);
			list2dnsndarr[myid] = new int [nlistr];
			if(!list2dnsndarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nlistr;kii++) {
				jj = listr[kii];
				jblk = listr2[kii];
				listdnsndarr[myid][kii] = jj;
				list2dnsndarr[myid][kii] = jblk;
			};
		} else {
			ni = iasnd[iproc+1]-iasnd[iproc];
			iarr = listsnd+iasnd[iproc]*2;
			nlistdnsndarr[iproc] = ni;
			listdnsndarr[iproc] = new int [ni];
			if(!listdnsndarr[iproc]) MemoryFail(funcname);
			list2dnsndarr[iproc] = new int [ni];
			if(!list2dnsndarr[iproc]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii*2];
				jblk = iarr[kii*2+1];
				listdnsndarr[iproc][kii] = jj;
				list2dnsndarr[iproc][kii] = jblk;
			};
		};
	};

	iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc == myid) {
			nlistdnrcvarr[myid] = nloc;
			listdnrcvarr[myid] = new int [nloc];
			if(!listdnrcvarr[myid]) MemoryFail(funcname);
			list2dnrcvarr[myid] = new int [nloc];
			if(!list2dnrcvarr[myid]) MemoryFail(funcname);
			for (kii=0;kii<nloc;kii++) {
				jj = listtot[kii];
				jblk = list2tot[kii];
				listdnrcvarr[myid][kii] = jj;
				list2dnrcvarr[myid][kii] = jblk;
			};
		} else {
			ni = ObjSizeRecv[iobj] / (2*sizeof (int));
			iarr = (int *) ObjRecv[iobj];
			iproc1 = CpuIDRecv[iobj];
			nlistdnrcvarr[iproc1] = ni;
			listdnrcvarr[iproc1] = new int [ni];
			if(!listdnrcvarr[iproc1]) MemoryFail(funcname);
			list2dnrcvarr[iproc1] = new int [ni];
			if(!list2dnrcvarr[iproc1]) MemoryFail(funcname);
			for (kii=0;kii<ni;kii++) {
				jj = iarr[kii*2];
				jblk = iarr[kii*2+1];
				listdnrcvarr[iproc1][kii] = jj;
				list2dnrcvarr[iproc1][kii] = jblk;
			};
			iobj++;
		};
	};

	for (iproc=0;iproc<nproc;iproc++) {
		ni = 0;
		nlistfctarr[iproc] = ni;
		listfctarr[iproc] = new int [ni];
		if(!listfctarr[iproc]) MemoryFail(funcname);
	};

// Compute bsx, bsx2, bspx and bspx2

	delete [] bsx;
	delete [] bsx2;
	delete [] bspx;
	delete [] bspx2;

	bsx = new int [nblksloc];
	if (!bsx) MemoryFail(funcname);
	bsx2 = new int [nblksloc];
	if (!bsx2) MemoryFail(funcname);
	bspx2 = new int [nblksloc];
	if (!bspx2) MemoryFail(funcname);

	nz = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		bsx[iblk] = nz;
		if (blk2cpu[iblk] == myid) {
			nz += pblks[iblk+1]-pblks[iblk];
		};
	};

	nlistblk = 0;

	for (i=0;i<nblksloc;i++) imaskblk[i] = -1;

	for (i=0;i<nlist;i++) {
		jj = listtot[i];
		jblk = list2tot[i];
		if (imaskblk[jblk] == -1) {
			listblk[nlistblk] = jblk;
			nlistblk++;
			imaskblk[jblk] = 1;
		};
	};

	qsort (listblk, nlistblk, sizeof(int), compint);

	for (i=0;i<nblksloc;i++) bspx2[i] = -1;

	nz = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		bspx2[iblk] = nz;
		nz += pblks[iblk+1]-pblks[iblk];
	};

	bspx = new int [nz];
	if (!bspx) MemoryFail(funcname);

	for (i=0;i<nz;i++) bspx[i] = -1;

	for (i=0;i<nlist;i++) {
		jj = listtot[i];
		jblk = list2tot[i];
		ibs = bspx2[jblk]+jj;
		bspx[ibs] = i;
	};

// Allocate the vector data

	CSVector temp(nlocext,1,_nrhs);

	px = temp;
	pax = temp;

	CSVector vectdummy;

	temp = vectdummy;

	int nadd = nlocext;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ni = iasnd[iproc+1]-iasnd[iproc];
			nadd += ni;
		};
	};

	int nadd1 = nlocext;

	iobj = 0;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ni = ObjSizeRecv[iobj] / (2*sizeof (int));
			nadd1 += ni;
			iobj++;
		};
	};

	if (nadd1 > nadd) nadd = nadd1;

	vectsndrcv = new double [nadd];
	if (!vectsndrcv) MemoryFail (funcname);

//	cout << " Nadd = " << nadd << endl;

// Free receive arrays

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Free work arrays

	delete [] listrblk;
	delete [] listcblk;
	delete [] imaskrblk;
	delete [] imaskcblk;
	delete [] ibsr;
	delete [] ibsc;
	delete [] imaskr;
	delete [] imaskc;
	delete [] listr;
	delete [] listc;
	delete [] listr2;
	delete [] listc2;
	delete [] listsnd;
	delete [] iasnd;
	delete [] iptr;
	delete [] ibsgl;
	delete [] imask;
	delete [] listtot;
	delete [] list2tot;
	delete [] imaskblk;
	delete [] listblk;

};

// Author: Kharchenko S.A.
// CMvmR: Destructor
//========================================================================================
CMvmR::~CMvmR () { // Destructor
//	std::cout << " On entry to CMvmR destructor " << std::endl;

	delete [] vectsndrcv;

//	std::cout << " On return from CMvmR destructor " << std::endl;
};

// Author: Kharchenko S.A.
// CSMatrixR: General matrix by vector multiplication
//========================================================================================
void CSMatrixR::MvmA (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication

//	const char *funcname = "MvmA_00";

	int i,j,jj;

	_ax.SetSVect (0.0e0);

	for (i=0;i<n;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			_ax.vect[i] += a[j] * _x.vect[jj];
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: General matrix by vector multiplication for chosen nodes
//========================================================================================
void CSMatrixR::MvmA (const CSVector &_x, CSVector &_ax, const bool *_barr) const { // General matrix by vector multiplication for chosen nodes

//	const char *funcname = "MvmA_00";

	int i,j,jj;

	_ax.SetSVect (0.0e0);

	for (i=0;i<n;i++) {
		if (_barr[i]) {
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				_ax.vect[i] += a[j] * _x.vect[jj];
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: General matrix by vector multiplication
//========================================================================================
void CSMatrixRS::MvmA (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication

//	const char *funcname = "MvmA_02";

	int i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv, irhs;

	double *px, *pax, *pa, *pablk;;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	_ax.SetSVect (0.0e0);

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0;i<nsupr;i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = sprndc[jj+1]-sprndc[jj];
				jbsv = sprndc[jj];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
						*pax += *pa * *px++;
						pa += blai;
					};
					pax++;
				};
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: General transposed matrix by vector multiplication
//========================================================================================
void CSMatrixR::MvmAt (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication

//	const char *funcname = "MvmAt_00";

	int i,j,jj;

	_ax.SetSVect (0.0e0);

	for (i=0;i<n;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			_ax.vect[jj] += a[j] * _x.vect[i];
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixRS: General transposed matrix by vector multiplication for matrix stored by columns
//========================================================================================
void CSMatrixRS::MvmAt (const CSVector &_x, CSVector &_ax) const { // General transposed matrix by vector multiplication

//	const char *funcname = "MvmAt_01";

	int irhs,i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	double *px, *pax, *pa, *pablk;

	_ax.SetSVect (0.0e0);

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0;i<nsupr;i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = sprndc[jj+1]-sprndc[jj];
				jbsv = sprndc[jj];
				px = _x.vect+irhs*nvx+ibsv;
				for (kii=0;kii<blai;kii++) {
					pax = _ax.vect+irhs*nvax+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[jbsv+kjj] += a[ibs+kjj*blai+kii] * _x.vect[ibsv+kii];
						*pax++ += *pa * *px;
						pa += blai;
					};
					px++;
				};
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Symmetric matrix by vector multiplication
//========================================================================================
void CSMatrixR::MvmASymm (const CSVector &_x, CSVector &_ax) const { // Symmetric matrix by vector multiplication

	int i,j,jj;

	_ax.SetSVect (0.0e0);

	for (i=0;i<n;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			_ax.vect[i] += a[j] * _x.vect[jj];
			if (i != jj) {
				_ax.vect[jj] += a[j] * _x.vect[i];
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixRB: Symmetric matrix by vector multiplication
//========================================================================================
void CSMatrixRB::MvmASymm (const CSVector &_x, CSVector &_ax) const { // Symmetric matrix by vector multiplication

	int i,j,jj,kii,kjj;

	int blal = bla;
	int bla_2l = bla_2;

	_ax.SetSVect (0.0e0);

	for (i=0;i<nsupr;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			for (kii=0;kii<blal;kii++) {
				for (kjj=0;kjj<blal;kjj++) {
					_ax.vect[i*blal+kii] += a[j*bla_2l+kjj*blal+kii] * _x.vect[jj*blal+kjj];
				};
			};
			if (i != jj) {
				for (kii=0;kii<blal;kii++) {
					for (kjj=0;kjj<blal;kjj++) {
						_ax.vect[jj*blal+kii] += a[j*bla_2l+kii*blal+kjj] * _x.vect[i*blal+kjj];
					};
				};
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixRS: General matrix by vector multiplication for matrix stored by columns
//========================================================================================
void CSMatrixRS::MvmACol (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns

//	const char *funcname = "MvmACol_00";

	int irhs,i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	double *px, *pax, *pa, *pablk;

	_ax.SetSVect (0.0e0);

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0;i<nsupc;i++) {
			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = sprndr[jj+1]-sprndr[jj];
				jbsv = sprndr[jj];
				px = _x.vect+irhs*nvx+ibsv;
				for (kii=0;kii<blai;kii++) {
					pax = _ax.vect+irhs*nvax+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[jbsv+kjj] += a[ibs+kjj*blai+kii] * _x.vect[ibsv+kii];
						*pax++ += *pa * *px;
						pa += blai;
					};
					px++;
				};
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixRS: General matrix by vector multiplication for matrix stored by columns in out-of-core mode
//========================================================================================
void CSMatrixRS::MvmACol (int _nblks, int *_blks, FILE **_mtrafiles, // General matrix by vector multiplication for matrix stored by columns in out-of-core mode
							const CSVector &_x, CSVector &_ax) const {

	const char *funcname = "MvmACol_01";

	int irhs,i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibs0, ibsv, jbsv;
	int iblk, nzloc, nzaloc;
	int i1, i2, ibeg, iend, jloc, jjloc;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

// Allocate work arrays for reading 

	const int nsupalloc = 20;

	int nsupmx = 0;

	for (i=0;i<nsupc;i++) {
		int nsupi = ia[i+1]-ia[i];
		if (nsupi>nsupmx) nsupmx = nsupi;
	};

	if (nsupmx > nsupalloc) nsupmx = nsupalloc;

	int sizealoc = nsupmx * blamx * blamx;

	double *aloc;

	aloc = new double [sizealoc];
	if (!aloc) MemoryFail (funcname);

	double *px, *pax, *pa, *pablk;

	_ax.SetSVect (0.0e0);

// Main cycle over the block columns

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {

			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];

// Reading cycle

			i1 = ia[i];
			i2 = ia[i+1]-1;

			iend = i1-1;

			while (iend < i2) {

				ibeg = iend + 1;
				nzaloc = 0;

				while (nzaloc <= sizealoc && iend <= i2) {

					jloc = iend + 1;
					jjloc = ja[jloc];
					blaj = sprndr[jjloc+1]-sprndr[jjloc];
					nzloc = blai*blaj;

					if (nzaloc+nzloc > sizealoc) break;

					nzaloc += nzloc;
					iend = jloc;

					if (iend == i2) break;

				};

// Read data of the supernode

				ibs0 = bsa[ibeg];

				FGet (_mtrafiles[iblk],nzaloc,aloc,ibs0);

// Multiply

				for (irhs=0;irhs<nrhsx;irhs++) {
					for (j=ibeg;j<=iend;j++) {
						jj = ja[j];
						ibs = bsa[j];
						pablk = aloc+ibs-ibs0;
						blaj = sprndr[jj+1]-sprndr[jj];
						jbsv = sprndr[jj];
						px = _x.vect+irhs*nvx+ibsv;
						for (kii=0;kii<blai;kii++) {
							pax = _ax.vect+irhs*nvax+jbsv;
							pa = pablk+kii;
							for (kjj=0;kjj<blaj;kjj++) {
//								_ax.vect[jbsv+kjj] += a[ibs+kjj*blai+kii] * _x.vect[ibsv+kii];
								*pax++ += *pa * *px;
								pa += blai;
							};
							px++;
						};
					};
				};

			};

		};
	};

// Delete work arrays

	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: General matrix by vector multiplication for matrix stored by columns
//========================================================================================
void CSMatrixRS::MvmAtCol (const CSVector &_x, CSVector &_ax) const { // General matrix by vector multiplication for matrix stored by columns

//	const char *funcname = "MvmAtCol_00";

	int i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv, irhs;

	double *px, *pax, *pa, *pablk;;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	_ax.SetSVect (0.0e0);

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0;i<nsupc;i++) {
			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];
			for (j=ia[i];j<ia[i+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = sprndr[jj+1]-sprndr[jj];
				jbsv = sprndr[jj];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
						*pax += *pa * *px++;
						pa += blai;
					};
					pax++;
				};
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixRS: General transposed matrix by vector multiplication for matrix stored by columns in out-of-core mode
//========================================================================================
void CSMatrixRS::MvmAtCol (int _nblks, int *_blks, FILE **_mtrafiles, // General transposed matrix by vector multiplication for matrix stored by columns in out-of-core mode
							const CSVector &_x, CSVector &_ax) const {

	const char *funcname = "MvmAtCol_01";

	int i,j,jj,kii,kjj, irhs;
	int blai, blaj, ibs, ibs0, ibsv, jbsv;
	int iblk, nzaloc, nzloc;
	int i1, i2, ibeg, iend, jloc, jjloc;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

// Allocate work arrays for reading 

	const int nsupalloc = 20;

	int nsupmx = 0;

	for (i=0;i<nsupc;i++) {
		int nsupi = ia[i+1]-ia[i];
		if (nsupi>nsupmx) nsupmx = nsupi;
	};

	if (nsupmx > nsupalloc) nsupmx = nsupalloc;

	int sizealoc = nsupmx * blamx * blamx;

	double *aloc;

	aloc = new double [sizealoc];
	if (!aloc) MemoryFail (funcname);

	double *px, *pax, *pa, *pablk;;

	_ax.SetSVect (0.0e0);

// Main cycle over the block columns

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {

			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];

// Reading cycle

			i1 = ia[i];
			i2 = ia[i+1]-1;

			iend = i1-1;

			while (iend < i2) {

				ibeg = iend + 1;
				nzaloc = 0;

				while (nzaloc <= sizealoc && iend <= i2) {

					jloc = iend + 1;
					jjloc = ja[jloc];
					blaj = sprndr[jjloc+1]-sprndr[jjloc];
					nzloc = blai*blaj;

					if (nzaloc+nzloc > sizealoc) break;

					nzaloc += nzloc;
					iend = jloc;

					if (iend == i2) break;

				};

// Read data of the supernode

				ibs0 = bsa[ibeg];

				FGet (_mtrafiles[iblk],nzaloc,aloc,ibs0);

// Multiply

				for (irhs=0;irhs<nrhsx;irhs++) {
					for (j=ibeg;j<=iend;j++) {
						jj = ja[j];
						ibs = bsa[j];
						pablk = aloc+ibs-ibs0;
						blaj = sprndr[jj+1]-sprndr[jj];
						jbsv = sprndr[jj];
						pax = _ax.vect+irhs*nvax+ibsv;
						for (kii=0;kii<blai;kii++) {
							px = _x.vect+irhs*nvx+jbsv;
							pa = pablk+kii;
							for (kjj=0;kjj<blaj;kjj++) {
//								_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
								*pax += *pa * *px++;
								pa += blai;
							};
							pax++;
						};
					};
				};
			};
		};
	};

// Delete work arrays

	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Symmetric matrix by vector multiplication
//========================================================================================
void CSMatrixRS::MvmASymm (const CSVector &_x, CSVector &_ax) const { // Symmetric matrix by vector multiplication

	int i,j,jj,kii,kjj;
	int blai, blaj, ibs, ibsv, jbsv;

	_ax.SetSVect (0.0e0);

	for (i=0;i<nsupr;i++) {
		blai = sprndr[i+1]-sprndr[i];
		ibsv = sprndr[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			ibs = bsa[j];
			blaj = sprndr[jj+1]-sprndr[jj];
			jbsv = sprndr[jj];
			for (kii=0;kii<blai;kii++) {
				for (kjj=0;kjj<blaj;kjj++) {
					_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
				};
			};
			if (i != jj) {
				for (kii=0;kii<blaj;kii++) {
					for (kjj=0;kjj<blai;kjj++) {
						_ax.vect[jbsv+kii] += a[ibs+kii*blai+kjj] * _x.vect[ibsv+kjj];
					};
				};
			};
		};
	};
};

// Author: Kharchenko S.A.
// CSMatrixR: Solve triangular system with L stored by columns
//========================================================================================
void CSMatrixR::SolveL (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns

	int i;

// Direct solve

	for (i=0;i<n;i++) _lx.vect[i] = _x.vect[i];

	for (i=0; i<n; i++) {
		int i1   = ia[i]+1;
		int i2   = ia[i+1];
		_lx.vect[i] *= a[ia[i]];
		for (int j=i1; j<i2; j++) {
			int k = ja[j];
			_lx.vect[k] -= a[j] * _lx.vect[i];
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Solve for chosen components triangular system with L stored by rows in the whole matrix stored by rows
//========================================================================================
void CSMatrixR::SolveLA (const CSVector &_diagl, const CSVector &_x, CSVector &_lx, const bool *_barr) const { // Solve for chosen components triangular system with L stored by rows in the whole matrix stored by rows

	int i;

// Direct solve

	for (i=0;i<n;i++) _lx.vect[i] = _x.vect[i];

	for (i=0; i<n; i++) {
		if (_barr[i]) {
			int i1   = ia[i];
			int i2   = ia[i+1];
			for (int j=i1; j<i2; j++) {
				int k = ja[j];
				if (k >= i) break;
				_lx.vect[i] -= a[j] * _lx.vect[k];
			};
			_lx.vect[i] *= _diagl.vect[i];
		} else {
			_lx.vect[i] = 0.0e0;
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixRB: Solve triangular system with L stored by columns
//========================================================================================
void CSMatrixRB::SolveL (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns

	const char *funcname = "SolveL_01";

	int i,kii,kjj,j;
	int nloc=n;

	int blal = bla;
	int bla_2l = bla_2;

	double *vloc;

	vloc = new double [blal];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	for (i=0;i<nloc;i++) _lx.vect[i] = _x.vect[i];

	for (i=0; i<nsupr; i++) {
		int i1 = ia[i]+1;
		int i2 = ia[i+1];
		j = ia[i];
		for (kii=0;kii<blal;kii++) vloc[kii] = 0.0e0;
		for (kii=0;kii<blal;kii++) {
			for (kjj=0;kjj<blal;kjj++) {
				vloc[kii] += a[j*bla_2l+kii*blal+kjj]*_lx.vect[i*blal+kjj];
			};
		};
		for (kii=0;kii<blal;kii++) _lx.vect[i*blal+kii] = vloc[kii];
		for (int j=i1; j<i2; j++) {
			int k = ja[j];
			for (kii=0;kii<blal;kii++) {
				for (kjj=0;kjj<blal;kjj++) {
					_lx.vect[k*blal+kii] -= a[j*bla_2l+kii*blal+kjj] * _lx.vect[i*blal+kjj];
				};
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with L stored by columns
//========================================================================================
void CSMatrixRS::SolveL (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by columns

	const char *funcname = "SolveL_02";

	int i,kii,kjj,j;
	int nloc=n;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	int blai, blaj, ibs, ibsv, jbsv, kiib, irhs;

	int nrhsx = _x.nrhs;
//	int nvx = _x.nv;
	int nvlx = _lx.nv;

	double *px, *plx, *pl, *pablk;

	px = _x.vect;
	plx = _lx.vect;

	for (i=0;i<nloc*nrhsx;i++) *plx++ = *px++;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=0; i<nsupr; i++) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			int i1 = ia[i]+1;
			int i2 = ia[i+1];
			j = ia[i];
			ibs = bsa[j];
			pablk = a+ibs;
			plx = vloc;
			for (kii=0;kii<blai;kii++) *plx++ = 0.0e0;
			plx = vloc;
			kiib = 0;
			for (kii=0;kii<blai;kii++) {
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk+kiib;
				for (kjj=0;kjj<blai;kjj++) {
//					vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
					*plx += *pl++ * *px++;
				};
				plx++;
				kiib += blai;
			};
			px = vloc;
			plx = _lx.vect+irhs*nvlx+ibsv;
			for (kii=0;kii<blai;kii++) *plx++ = *px++;
			for (int j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = sprndr[k+1]-sprndr[k];
				jbsv = sprndr[k];
				ibs = bsa[j];
				pablk = a+ibs;
				plx = _lx.vect+irhs*nvlx+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
						*plx -= *pl++ * *px++;
					};
					plx++;
					kiib += blai;
				};
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with L stored by columns in out-of-core mode
//========================================================================================
void CSMatrixRS::SolveL (int _nblks, int *_blks, FILE **_mtrlfiles, // Solve triangular system with L stored by columns in out-of-core mode
							const CSVector &_x, CSVector &_lx) const {

	const char *funcname = "SolveL_03";

	int i,kii,kjj,j;
	int nloc=n;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Allocate work arrays for reading 

	int nsupmx = 0;

	for (i=0;i<nsupc;i++) {
		int nsupi = ia[i+1]-ia[i];
		if (nsupi>nsupmx) nsupmx = nsupi;
	};

	double *aloc;

	aloc = new double [nsupmx*blamx*blamx];
	if (!aloc) MemoryFail (funcname);

// Direct solve

	int blai, blaj, ibs, ibs0, ibsv, jbsv, kiib;
	int iblk, irhs;

	double *px, *plx, *pl, *pablk;

	int nrhsx = _x.nrhs;
//	int nvx = _x.nv;
	int nvlx = _lx.nv;

	px = _x.vect;
	plx = _lx.vect;

	for (i=0;i<nloc*nrhsx;i++) *plx++ = *px++;

// Main cycle over the block columns

	for (iblk=0;iblk<_nblks;iblk++) {
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {

			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];

// Read data of the supernode

			ibs0 = bsa[ia[i]];

			int nzloc = 0;

			for (j=ia[i];j<ia[i+1];j++) {
				int jj = ja[j];
				blaj = sprndr[jj+1]-sprndr[jj];
				nzloc += blai*blaj;
			};

			FGet (_mtrlfiles[iblk],nzloc,aloc,ibs0);

// Solve

			for (irhs=0;irhs<nrhsx;irhs++) {
				int i1 = ia[i]+1;
				int i2 = ia[i+1];
				j = ia[i];
				ibs = bsa[j];
				pablk = aloc+ibs-ibs0;
				plx = vloc;
				for (kii=0;kii<blai;kii++) *plx++ = 0.0e0;
				plx = vloc;
				kiib = 0;
				for (kii=0;kii<blai;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
						*plx += *pl++ * *px++;
					};
					plx++;
					kiib += blai;
				};
				px = vloc;
				plx = _lx.vect+irhs*nvlx+ibsv;
				for (kii=0;kii<blai;kii++) *plx++ = *px++;
				for (int j=i1; j<i2; j++) {
					int k = ja[j];
					blaj = sprndr[k+1]-sprndr[k];
					jbsv = sprndr[k];
					ibs = bsa[j];
					pablk = aloc+ibs-ibs0;
					plx = _lx.vect+irhs*nvlx+jbsv;
					kiib = 0;
					for (kii=0;kii<blaj;kii++) {
						px = _lx.vect+irhs*nvlx+ibsv;
						pl = pablk+kiib;
						for (kjj=0;kjj<blai;kjj++) {
//							_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
							*plx -= *pl++ * *px++;
						};
						plx++;
						kiib += blai;
					};
				};
			};
		};
	};

	delete [] vloc;
	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Solve triangular system with U stored by rows
//========================================================================================
void CSMatrixR::SolveU (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows

	int i;

// Backward solve

	for (i=n-1; i>=0; i--) {
		int i1   = ia[i]+1;
		int i2   = ia[i+1];
		double s = _x.vect[i];
		for (int j=i1; j<i2; j++) {
			int k = ja[j];
			s    -= a[j] * _ux.vect[k];
		};
		_ux.vect[i] = s * a[ia[i]];
	};
};

// Author: Kharchenko S.A.
// CSMatrixRB: Solve triangular system with U stored by rows
//========================================================================================
void CSMatrixRB::SolveU (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows

	const char *funcname = "SolveU_01";

	int i, kii, kjj;
//	int nloc=n;

	int blal = bla;
	int bla_2l = bla_2;

	double *vloc;

	vloc = new double [blal];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	for (i=nsupr-1;i>=0;i--) {
		int i1 = ia[i]+1;
		int i2 = ia[i+1];
		for (kii=0;kii<blal;kii++) vloc[kii] = _x.vect[i*blal+kii];
		int j;
		for (j=i1; j<i2; j++) {
			int k = ja[j];
			for (kii=0;kii<blal;kii++) {
				for (kjj=0;kjj<blal;kjj++) {
					vloc[kii] -= a[j*bla_2l+kjj*blal+kii] * _ux.vect[k*blal+kjj];
				};
			};
		};
		j = ia[i];
		for (kii=0;kii<blal;kii++) _ux.vect[i*blal+kii] = 0.0e0;
		for (kii=0;kii<blal;kii++) {
			for (kjj=0;kjj<blal;kjj++) {
				_ux.vect[i*blal+kii] += a[j*bla_2l+kjj*blal+kii] * vloc[kjj];
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with U stored by rows
//========================================================================================
void CSMatrixRS::SolveU (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by rows

	const char *funcname = "SolveU_02";

	int i, kii, kjj;
//	int nloc=n;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	int blai, blaj, ibs, ibsv, jbsv, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	double *px, *pux, *pu, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (i=nsupr-1;i>=0;i--) {
			blai = sprndr[i+1]-sprndr[i];
			ibsv = sprndr[i];
			int i1 = ia[i]+1;
			int i2 = ia[i+1];
			pux = vloc;
			px = _x.vect+irhs*nvx+ibsv;
			for (kii=0;kii<blai;kii++) {
//				vloc[kii] = _x.vect[ibsv+kii];
				*pux++ = *px++;
			};
			int j;
			for (j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = sprndr[k+1]-sprndr[k];
				jbsv = sprndr[k];
				ibs = bsa[j];
				pablk = a+ibs;
				pux = vloc;
				for (kii=0;kii<blai;kii++) {
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
						*pux -= *pu * *px++;
						pu += blai;
					};
					pux++;
				};
			};
			j = ia[i];
			ibs = bsa[j];
			pablk = a+ibs;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) *pux++ = 0.0e0;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) {
				px = vloc;
				pu = pablk+kii;
				for (kjj=0;kjj<blai;kjj++) {
//					_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
					*pux += *pu * *px++;
					pu += blai;
				};
				pux++;
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with U stored by rows in out-of-core mode
//========================================================================================
void CSMatrixRS::SolveU (int _nblks, int *_blks, FILE **_mtrufiles, // Solve triangular system with U stored by rows in out-of-core mode
							const CSVector &_x, CSVector &_ux) const {

	const char *funcname = "SolveU_03";

	int i, kii, kjj;
//	int nloc=n;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Allocate work arrays for reading 

	int nsupmx = 0;

	for (i=0;i<nsupc;i++) {
		int nsupi = ia[i+1]-ia[i];
		if (nsupi>nsupmx) nsupmx = nsupi;
	};

	double *aloc;

	aloc = new double [nsupmx*blamx*blamx];
	if (!aloc) MemoryFail (funcname);

// Backward solve

	int blai, blaj, ibs, ibs0, ibsv, jbsv;
	int iblk, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	double *px, *pux, *pu, *pablk;

	for (iblk=_nblks-1;iblk>=0;iblk--) {
		for (i=_blks[iblk+1]-1;i>=_blks[iblk];i--) {

			blai = sprndc[i+1]-sprndc[i];
			ibsv = sprndc[i];

// Read data of the supernode

			ibs0 = bsa[ia[i]];

			int nzloc = 0;

			int j;
			for (j=ia[i];j<ia[i+1];j++) {
				int jj = ja[j];
				blaj = sprndr[jj+1]-sprndr[jj];
				nzloc += blai*blaj;
			};

			FGet (_mtrufiles[iblk],nzloc,aloc,ibs0);

// Solve

			for (irhs=0;irhs<nrhsx;irhs++) {
				int i1 = ia[i]+1;
				int i2 = ia[i+1];
				pux = vloc;
				px = _x.vect+irhs*nvx+ibsv;
				for (kii=0;kii<blai;kii++) {
//					vloc[kii] = _x.vect[ibsv+kii];
					*pux++ = *px++;
				};
				for (j=i1; j<i2; j++) {
					int k = ja[j];
					blaj = sprndr[k+1]-sprndr[k];
					jbsv = sprndr[k];
					ibs = bsa[j];
					pablk = aloc+ibs-ibs0;
					pux = vloc;
					for (kii=0;kii<blai;kii++) {
						px = _ux.vect+irhs*nvux+jbsv;
						pu = pablk+kii;
						for (kjj=0;kjj<blaj;kjj++) {
//							vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
							*pux -= *pu * *px++;
							pu += blai;
						};
						pux++;
					};
				};
				j = ia[i];
				ibs = bsa[j];
				pablk = aloc+ibs-ibs0;
				pux = _ux.vect+irhs*nvux+ibsv;
				for (kii=0;kii<blai;kii++) *pux++ = 0.0e0;
				pux = _ux.vect+irhs*nvux+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = vloc;
					pu = pablk+kii;
					for (kjj=0;kjj<blai;kjj++) {
//						_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
						*pux += *pu * *px++;
						pu += blai;
					};
					pux++;
				};
			};
		};
	};

	delete [] vloc;
	delete [] aloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Multiply by the current block row of Au
//========================================================================================
void CSMatrixRS::MvmBlockRowAu (int *_sprnds, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const {

//	const char *funcname = "MvmBlockRowAu";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsax[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprnds[jj+1]-_sprnds[jj];
				jbsv = _bsx[jj];
				pax = _ax.vect+irhs*nvax+ibsv;
				for (kii=0;kii<blai;kii++) {
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						_ax.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * _x.vect[jbsv+kjj];
//						*pax += *pa * *px++;
						*pax += *pa * *px++;
						pa += blai;
					};
					pax++;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Multiply by the current block row of Au
//========================================================================================
void CSMatrixRS::MvmPointBlockRowAu (int *_sprnds, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const {

//	const char *funcname = "MvmPointBlockRowAu";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj;
	int ilist, ibs, ibsv, jbsv, irhs;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				ibsv = _bsax[i];
				for (j=ia[ilist];j<ia[ilist+1];j++) {
					jj = ja[j];
					ibs = bsa[j];
					pablk = a+ibs;
					jbsv = _bsx[jj];
					pax = _ax.vect+irhs*nvax+ibsv;
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk;
					*pax += *pa * *px;
				};
			};
		};
	} else {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			ibsv = _bsax[i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jbsv = _bsx[jj];
				_ax.vect[ibsv] += a[j] * _x.vect[jbsv];
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Multiply by the current block row of Au
//========================================================================================
void CSMatrixR::MvmBlockRowAu2Index (// Multiply by the current block row of Au
									int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const {

//	const char *funcname = "MvmBlockRowAu2Index";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i, j, jj;
	int ilist, ibsv, jbsv, irhs, iblk, jblk;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsax[iblk]+i;
				for (j=ia[ilist];j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2[j];
					jbsv = _bsx[jblk]+jj;
					pablk = a+j;
					pax = _ax.vect+irhs*nvax+ibsv;
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk;
					*pax += *pa * *px;
				};
			};
		};
	} else {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			iblk = list2[ilist];
			ibsv = _bsax[iblk]+i;
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jblk = ja2[j];
				jbsv = _bsx[jblk]+jj;
				_ax.vect[ibsv] += a[j] * _x.vect[jbsv];
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Multiply by the current block row of Au
//========================================================================================
void CSMatrixR::MvmBlockRowAu2Index2 (// Multiply by the current block row of Au
									int *_bsx, int *_bsx2, const CSVector &_x, int *_bsax, int *_bsax2, CSVector &_ax) const {

//	const char *funcname = "MvmBlockRowAu2Index2";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i, j, jj;
	int ilist, ibsv, jbsv, irhs, iblk, jblk;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsax[_bsax2[iblk]+i];
				for (j=ia[ilist];j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2[j];
					jbsv = _bsx[_bsx2[jblk]+jj];
					pablk = a+j;
					pax = _ax.vect+irhs*nvax+ibsv;
					px = _x.vect+irhs*nvx+jbsv;
					pa = pablk;
					*pax += *pa * *px;
				};
			};
		};
	} else {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			iblk = list2[ilist];
			ibsv = _bsax[_bsax2[iblk]+i];
			for (j=ia[ilist];j<ia[ilist+1];j++) {
				jj = ja[j];
				jblk = ja2[j];
				jbsv = _bsx[_bsx2[jblk]+jj];
				_ax.vect[ibsv] += a[j] * _x.vect[jbsv];
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Multiply by the current block column of Al with/without diagonal supernode
//========================================================================================
void CSMatrixRS::MvmBlockColumnAl (bool _bdiag, // Multiply by the current block column of Al with/without diagonal supernode
									int *_sprnds,
									int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const {

//	const char *funcname = "MvmBlockColumnAl";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,kii,kjj,ibeg;
	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs, kiib;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0;ilist<nlist;ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			if (_bdiag) {
				ibeg = ia[ilist];
			} else {
				ibeg = ia[ilist]+1;
			};
			for (j=ibeg;j<ia[ilist+1];j++) {
				jj = ja[j];
				ibs = bsa[j];
				pablk = a+ibs;
				blaj = _sprnds[jj+1]-_sprnds[jj];
				jbsv = _bsax[jj];
				pax = _ax.vect+irhs*nvax+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_ax.vect[jbsv+kii] += a[ibs+kii*blai+kjj] * _x.vect[ibsv+kjj];
//						*pax += *pa++ * *px++;
						*pax += *pa++ * *px++;
					};
					pax++;
					kiib += blai;
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixRS: Multiply by the current block column of Al with/without diagonal supernode
//========================================================================================
void CSMatrixRS::MvmPointBlockColumnAl (bool _bdiag, // Multiply by the current block column of Al with/without diagonal supernode
														int *_sprnds,
														int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const {

//	const char *funcname = "MvmPointBlockColumnAl";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,ibeg;
	int ilist, ibs, ibsv, jbsv, irhs;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				ibsv = _bsx[i];
				if (_bdiag) {
					ibeg = ia[ilist];
				} else {
					ibeg = ia[ilist]+1;
				};
				for (j=ibeg;j<ia[ilist+1];j++) {
					jj = ja[j];
					ibs = bsa[j];
					pablk = a+ibs;
					jbsv = _bsax[jj];
					pax = _ax.vect+irhs*nvax+jbsv;
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk;
					*pax += *pa * *px;
				};
			};
		};
	} else {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				ibsv = _bsx[i];
				if (_bdiag) {
					ibeg = ia[ilist];
				} else {
					ibeg = ia[ilist]+1;
				};
				for (j=ibeg;j<ia[ilist+1];j++) {
					jj = ja[j];
					jbsv = _bsax[jj];
					_ax.vect[jbsv] += a[j] * _x.vect[ibsv];
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Multiply by the current block column of Al with/without diagonal supernode
//========================================================================================
void CSMatrixR::MvmBlockColumnAl2Index (bool _bdiag, // Multiply by the current block column of Al with/without diagonal supernode
										int *_bsx, const CSVector &_x, int *_bsax, CSVector &_ax) const {

//	const char *funcname = "MvmBlockColumnAl2Index";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,ibeg;
	int ilist, ibsv, jbsv, irhs, iblk, jblk;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsx[iblk]+i;
				if (_bdiag) {
					ibeg = ia[ilist];
				} else {
					ibeg = ia[ilist]+1;
				};
				for (j=ibeg;j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2[j];
					pablk = a+j;
					jbsv = _bsax[jblk]+jj;
					pax = _ax.vect+irhs*nvax+jbsv;
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk;
					*pax += *pa * *px;
				};
			};
		};
	} else {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsx[iblk]+i;
				if (_bdiag) {
					ibeg = ia[ilist];
				} else {
					ibeg = ia[ilist]+1;
				};
				for (j=ibeg;j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2[j];
					jbsv = _bsax[jblk]+jj;
					_ax.vect[jbsv] += a[j] * _x.vect[ibsv];
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Multiply by the current block column of Al with/without diagonal supernode
//========================================================================================
void CSMatrixR::MvmBlockColumnAl2Index2 (bool _bdiag, // Multiply by the current block column of Al with/without diagonal supernode
										int *_bsx, int *_bsx2, const CSVector &_x, int *_bsax, int *_bsax2, CSVector &_ax) const {

//	const char *funcname = "MvmBlockColumnAl2Index2";

// Attention: ax on entry should be initialized by zeroes if necessary

	int i,j,jj,ibeg;
	int ilist, ibsv, jbsv, irhs, iblk, jblk;

	double *px, *pax, *pa, *pablk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvax = _ax.nv;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsx[_bsx2[iblk]+i];
				if (_bdiag) {
					ibeg = ia[ilist];
				} else {
					ibeg = ia[ilist]+1;
				};
				for (j=ibeg;j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2[j];
					pablk = a+j;
					jbsv = _bsax[_bsax2[jblk]+jj];
					pax = _ax.vect+irhs*nvax+jbsv;
					px = _x.vect+irhs*nvx+ibsv;
					pa = pablk;
					*pax += *pa * *px;
				};
			};
		};
	} else {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0;ilist<nlist;ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsx[_bsx2[iblk]+i];
				if (_bdiag) {
					ibeg = ia[ilist];
				} else {
					ibeg = ia[ilist]+1;
				};
				for (j=ibeg;j<ia[ilist+1];j++) {
					jj = ja[j];
					jblk = ja2[j];
					jbsv = _bsax[_bsax2[jblk]+jj];
					_ax.vect[jbsv] += a[j] * _x.vect[ibsv];
				};
			};
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: General matrix by vector multiplication
//========================================================================================
void CGSMatrixRS::MvmA (const CGSMatrixRS &_gmtral, const CGSMatrixRS &_gmtrau, // General matrix by vector multiplication
						const CSVector &_x, CSVector &_ax) const {

//	const char *funcname = "MvmA_00";

	double dzero = 0.0e0;

// Initialize ax by zeroes

	int nrhsx = _x.nrhs;
	int nvax = _ax.nv;

	for (int kii=0;kii<nrhsx*nvax;kii++) _ax.vect[kii] = dzero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<_gmtral.nblksr;iblk++) {

		_gmtrau.mtrarr[iblk].MvmBlockRowAu    (_gmtral.sprndr, _gmtral.sprndr, _x, _gmtral.sprndr, _ax);

		_gmtral.mtrarr[iblk].MvmBlockColumnAl (false, _gmtral.sprndr, _gmtral.sprndr, _x, _gmtral.sprndr, _ax);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: General matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixRS::MvmA (const CTree &_tree, CMvmR &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixRS &_gmtral, CGSMatrixRS &_gmtrau, 
							const CSVector &_x, CSVector &_ax) const {

//	const char *funcname = "MvmA_01";

	double dzero = 0.0e0;

// Determine the case

	int icase = 1;
	if (_gmtral.n == _gmtral.nsupr) icase = 0;

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = dzero;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.pax.vect[kii] = dzero;

// Initialize px by x data

	int irhs;
	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_mvm.px.vect[irhs*_mvm.nlocext+kii] = _x.vect[irhs*_mvm.nloc+kii];
		};
	};

// Search the tree

	int inode, inodetype;
	int fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;
	int ilev;

	for (ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode1;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (false, _tree.comm, 'A', '=', _tree.myid, fathercpuid, _tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (false, _tree.comm, 'R', '=', _tree.myid, childcpuid, _tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (false, _tree.comm, 'r', _tree.myid, childcpuid, _tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (false, _tree.comm, 'A', _tree.myid, fathercpuid2, _tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

NextNode1:;

	};

// Perform all local computations

	for (int ilist=0;ilist<_gmtral.nlist;ilist++) {
		int iblk = _gmtral.listb[ilist];

		_gmtrau.ReadBlock (iblk);
		if (icase == 1) {
			_gmtrau.mtrarr[iblk].MvmBlockRowAu    (_gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		} else {
			_gmtrau.mtrarr[iblk].MvmPointBlockRowAu    (_gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		};
		_gmtrau.FreeBlock (iblk);

		_gmtral.ReadBlock (iblk);
		if (icase == 1) {
			_gmtral.mtrarr[iblk].MvmBlockColumnAl (false, _gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		} else {
			_gmtral.mtrarr[iblk].MvmPointBlockColumnAl (false, _gmtral.sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		};
		_gmtral.FreeBlock (iblk);

	};

// Search the tree one more time

	nodebeg = _tree.cpuidbeg[_tree.myid];
	nodeend = _tree.cpuidend[_tree.myid];

	ilevbeg = _tree.nodes[nodebeg].nodelv;
	ilevend = _tree.nodes[nodeend].nodelv;

	for (ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode2;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (false, _tree.comm, 'A', '+', _tree.myid, childcpuid, 2*_tree.nnodes+inode,
										_mvm.pax,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (false, _tree.comm, 'A', '=', _tree.myid, fathercpuid2, 2*_tree.nnodes+inode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (false, _tree.comm, 'A', _tree.myid, fathercpuid, 2*_tree.nnodes+fathernode,
								_mvm.pax,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (false, _tree.comm, 'R', _tree.myid, childcpuid, 2*_tree.nnodes+childnode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

NextNode2:;

	};

// Store the result of computations in ax

	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ax.vect[irhs*_mvm.nloc+kii] = _mvm.pax.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: General matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixR::MvmA2Index (const CTree &_tree, CMvmR &_mvm, // General matrix by vector multiplication in parallel mode
							CGSMatrixR &_gmtral, CGSMatrixR &_gmtrau, 
							const CSVector &_x, CSVector &_ax) const {

//	const char *funcname = "MvmA2Index";

#ifdef __PMPITRACE__
	static int i_set = 0;
	static int mvma_beg;
	static int mvma_end;
	static int mvmaloc_beg;
	static int mvmaloc_end;

	if (i_set == 0) {

		mvma_beg = MPE_Log_get_event_number ();
		mvma_end = MPE_Log_get_event_number ();
		MPE_Describe_state (mvma_beg, mvma_end, "MvmA", "gold");

		mvmaloc_beg = MPE_Log_get_event_number ();
		mvmaloc_end = MPE_Log_get_event_number ();
		MPE_Describe_state (mvmaloc_beg, mvmaloc_end, "MvmALoc", "brown");

		i_set = 1;

	};

	MPE_Log_event (mvma_beg, 0, NULL);
#endif

	double dzero = 0.0e0;

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = dzero;
	for (kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.pax.vect[kii] = dzero;

// Initialize px by x data

	int irhs;
	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_mvm.px.vect[irhs*_mvm.nlocext+kii] = _x.vect[irhs*_mvm.nloc+kii];
		};
	};

// Search the tree

	int inode, inodetype;
	int fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;
	int ilev;

	for (ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode1;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (true, _tree.comm, 'A', '=', _tree.myid, fathercpuid, _tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (true, _tree.comm, 'R', '=', _tree.myid, childcpuid, _tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (true, _tree.comm, 'r', _tree.myid, childcpuid, _tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (true, _tree.comm, 'A', _tree.myid, fathercpuid2, _tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

NextNode1:;

	};

// Perform all local computations

#ifdef __PMPITRACE__
	MPE_Log_event (mvmaloc_beg, 0, NULL);
#endif

	for (int ilist=0;ilist<_gmtral.nlist;ilist++) {
		int iblk = _gmtral.listb[ilist];

		_gmtrau.ReadBlock (iblk);
		_gmtrau.mtrarr[iblk].MvmBlockRowAu2Index (_mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtrau.FreeBlock (iblk);

		_gmtral.ReadBlock (iblk);
		_gmtral.mtrarr[iblk].MvmBlockColumnAl2Index (false, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		_gmtral.FreeBlock (iblk);

	};

#ifdef __PMPITRACE__
	MPE_Log_event (mvmaloc_end, 0, NULL);
#endif

// Search the tree one more time

	nodebeg = _tree.cpuidbeg[_tree.myid];
	nodeend = _tree.cpuidend[_tree.myid];

	ilevbeg = _tree.nodes[nodebeg].nodelv;
	ilevend = _tree.nodes[nodeend].nodelv;

	for (ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode2;

		inodetype = _tree.nodes[inode].nodetype;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (true, _tree.comm, 'A', '+', _tree.myid, childcpuid, 2*_tree.nnodes+inode,
										_mvm.pax,
										nlistloc, listloc,
										_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (true, _tree.comm, 'A', '=', _tree.myid, fathercpuid2, 2*_tree.nnodes+inode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (true, _tree.comm, 'A', _tree.myid, fathercpuid, 2*_tree.nnodes+fathernode,
								_mvm.pax,
								nlistloc, listloc,
								_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (true, _tree.comm, 'R', _tree.myid, childcpuid, 2*_tree.nnodes+childnode,
									_mvm.pax,
									nlistloc, listloc,
									_gmtral.blksr, _gmtral.sprndr, _gmtral.bl2ndr);

				};

			};

		};

NextNode2:;

	};

// Store the result of computations in ax

	for (irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ax.vect[irhs*_mvm.nloc+kii] = _mvm.pax.vect[irhs*_mvm.nlocext+kii];
		};
	};

#ifdef __PMPITRACE__
	MPE_Log_event (mvma_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CGSMatrixRS: General secondary matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixRS::MvmA2 (const CTree &_tree, CMvmR &_mvm, // General secondary matrix by vector multiplication in parallel mode
							const CSVector &_x, CSVector &_ax) {

	const char *funcname = "MvmA2";

// Remark: on entry to the routine vector _ax is assumed to be initialized

	double dzero = 0.0e0;

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Init the result

//	_ax.SetSVect (dzero);

// Determine the case

	int icase = 1;
	if (n == nsupr) icase = 0;

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.px.vect[kii] = dzero;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.pax.vect[kii] = dzero;

// Initialize px by local x data

	int kjj, jj, ibsg, ibsl, blai;

	for (kii=0;kii<_mvm.nlistupsndarr[myid];kii++) {
		jj = _mvm.listupsndarr[myid][kii];
		ibsg = _mvm.bsx[jj];
		ibsl = _mvm.bspx[jj];
		blai = sprndr[jj+1]-sprndr[jj];
		for (kjj=0;kjj<blai;kjj++) {
			_mvm.px.vect[ibsl+kjj] = _x.vect[ibsg+kjj];
		};
	};

// Exchange the data

	int NObjSend = nproc-1;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	int nz = 0;

	int iproc, nzloc;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ObjTypeSend[iobj] = 1;
			ObjIDSend[iobj] = myid;
			CpuIDSend[iobj] = iproc;
			ObjSend[iobj] = (char *) (_mvm.vectsndrcv+nz);
			nzloc = 0;
			for (kii=0;kii<_mvm.nlistupsndarr[iproc];kii++) {
				jj = _mvm.listupsndarr[iproc][kii];
				ibsg = _mvm.bsx[jj];
				ibsl = _mvm.bspx[jj];
				blai = sprndr[jj+1]-sprndr[jj];
				for (kjj=0;kjj<blai;kjj++) {
					_mvm.vectsndrcv[nz+nzloc+kjj] = _x.vect[ibsg+kjj];
				};
				nzloc += blai;
			};
			ObjSizeSend[iobj] = nzloc * sizeof (double);
			nz += nzloc;
			iobj++;
		};
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info = CMPIExchange::DataExchangeMPI (pcomm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
	if(info) throw " CGSMatrixRS::MvmA2: Error in DataExchangeMPI";

	double *darr;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		iproc = ObjIDRecv[iobj];
		darr = (double *) ObjRecv[iobj];
		nzloc = 0;
		for (kii=0;kii<_mvm.nlistuprcvarr[iproc];kii++) {
			jj = _mvm.listuprcvarr[iproc][kii];
			ibsg = _mvm.bsx[jj];
			ibsl = _mvm.bspx[jj];
			blai = sprndr[jj+1]-sprndr[jj];
			for (kjj=0;kjj<blai;kjj++) {
				_mvm.px.vect[ibsl+kjj] = darr[nzloc+kjj];
			};
			nzloc += blai;
		};
	};

// Free work arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Perform all local computations

	for (int ilist=0;ilist<nlist;ilist++) {
		int iblk = listb[ilist];

		ReadBlock (iblk);
		if (icase == 1) {
			mtrarr[iblk].MvmBlockRowAu      (sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		} else {
			mtrarr[iblk].MvmPointBlockRowAu (sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		};
		FreeBlock (iblk);

	};

// Store the result of computations in ax

	for (kii=0;kii<_mvm.nlistuprcvarr[myid];kii++) {
		jj = _mvm.listuprcvarr[myid][kii];
		ibsg = _mvm.bsx[jj];
		ibsl = _mvm.bspx[jj];
		blai = sprndr[jj+1]-sprndr[jj];
		for (kjj=0;kjj<blai;kjj++) {
			_ax.vect[ibsg+kjj] += _mvm.pax.vect[ibsl+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: General transposed secondary matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixRS::MvmAh2 (const CTree &_tree, CMvmR &_mvm, // General transposed secondary matrix by vector multiplication in parallel mode
									const CSVector &_x, CSVector &_ax) {

	const char *funcname = "MvmAh2";

// Remark: on entry to the routine vector _ax is assumed to be initialized

	double dzero = 0.0e0;

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Init the result

//	_ax.SetSVect (dzero);

// Determine the case

	int icase = 1;
	if (n == nsupr) icase = 0;

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.px.vect[kii] = dzero;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.pax.vect[kii] = dzero;

// Initialize px by local x data

	int kjj, jj, ibsg, ibsl, blai;

	for (kii=0;kii<_mvm.nlistdnsndarr[myid];kii++) {
		jj = _mvm.listdnsndarr[myid][kii];
		ibsg = _mvm.bsx[jj];
		ibsl = _mvm.bspx[jj];
		blai = sprndr[jj+1]-sprndr[jj];
		for (kjj=0;kjj<blai;kjj++) {
			_mvm.px.vect[ibsl+kjj] = _x.vect[ibsg+kjj];
		};
	};

// Perform all local computations

	for (int ilist=0;ilist<nlist;ilist++) {
		int iblk = listb[ilist];

		ReadBlock (iblk);
		if (icase == 1) {
			mtrarr[iblk].MvmBlockColumnAl      (true, sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		} else {
			mtrarr[iblk].MvmPointBlockColumnAl (true, sprndr, _mvm.bspx, _mvm.px, _mvm.bspx, _mvm.pax);
		};
		FreeBlock (iblk);

	};

// Exchange the resulting data

	int NObjSend = nproc-1;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	int nz = 0;
	int nzloc;

	int iproc;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ObjTypeSend[iobj] = 1;
			ObjIDSend[iobj] = myid;
			CpuIDSend[iobj] = iproc;
			ObjSend[iobj] = (char *) (_mvm.vectsndrcv+nz);
			nzloc = 0;
			for (kii=0;kii<_mvm.nlistdnsndarr[iproc];kii++) {
				jj = _mvm.listdnsndarr[iproc][kii];
				ibsg = _mvm.bsx[jj];
				ibsl = _mvm.bspx[jj];
				blai = sprndr[jj+1]-sprndr[jj];
				for (kjj=0;kjj<blai;kjj++) {
					_mvm.vectsndrcv[nz+nzloc+kjj] = _mvm.pax.vect[ibsl+kjj];
				};
				nzloc += blai;
			};
			ObjSizeSend[iobj] = nzloc * sizeof (double);
			nz += nzloc;
			iobj++;
		};
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info = CMPIExchange::DataExchangeMPI (pcomm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
	if(info) throw " CGSMatrixRS::MvmA2: Error in DataExchangeMPI";

	double *darr;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		iproc = ObjIDRecv[iobj];
		darr = (double *) ObjRecv[iobj];
		nzloc = 0;
		for (kii=0;kii<_mvm.nlistdnrcvarr[iproc];kii++) {
			jj = _mvm.listdnrcvarr[iproc][kii];
			ibsg = _mvm.bsx[jj];
			ibsl = _mvm.bspx[jj];
			blai = sprndr[jj+1]-sprndr[jj];
			for (kjj=0;kjj<blai;kjj++) {
				_mvm.pax.vect[ibsl+kjj] += darr[nzloc+kjj];
			};
			nzloc += blai;
		};
	};

// Free work arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Store the result of computations in ax

	for (kii=0;kii<_mvm.nlistdnrcvarr[myid];kii++) {
		jj = _mvm.listdnrcvarr[myid][kii];
		ibsg = _mvm.bsx[jj];
		ibsl = _mvm.bspx[jj];
		blai = sprndr[jj+1]-sprndr[jj];
		for (kjj=0;kjj<blai;kjj++) {
			_ax.vect[ibsg+kjj] += _mvm.pax.vect[ibsl+kjj];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: General secondary matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixR::MvmA2Index2 (const CTree &_tree, CMvmR &_mvm, // General secondary matrix by vector multiplication in parallel mode
							const CSVector &_x, CSVector &_ax) {

	const char *funcname = "MvmA2Index2";

#ifdef __PMPITRACE__
	static int i_set = 0;
	static int mvma2_beg;
	static int mvma2_end;

	if (i_set == 0) {

		mvma2_beg = MPE_Log_get_event_number ();
		mvma2_end = MPE_Log_get_event_number ();
		MPE_Describe_state (mvma2_beg, mvma2_end, "MvmA2", "green");

		i_set = 1;

	};

	MPE_Log_event (mvma2_beg, 0, NULL);
#endif

// Remark: on entry to the routine vector _ax is assumed to be initialized

	double dzero = 0.0e0;

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Init the result

//	_ax.SetSVect (dzero);

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.px.vect[kii] = dzero;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.pax.vect[kii] = dzero;

// Initialize px by local x data

	int jj, ibsg, ibsl, jblk;

	for (kii=0;kii<_mvm.nlistupsndarr[myid];kii++) {
		jblk = _mvm.list2upsndarr[myid][kii];
		jj = _mvm.listupsndarr[myid][kii];
		ibsg = _mvm.bsx[jblk]+jj;
		ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
		_mvm.px.vect[ibsl] = _x.vect[ibsg];
	};

// Exchange the data

	int NObjSend = nproc-1;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	int nz = 0;

	int iproc, nzloc;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ObjTypeSend[iobj] = 1;
			ObjIDSend[iobj] = myid;
			CpuIDSend[iobj] = iproc;
			ObjSend[iobj] = (char *) (_mvm.vectsndrcv+nz);
			nzloc = 0;
			for (kii=0;kii<_mvm.nlistupsndarr[iproc];kii++) {
				jblk = _mvm.list2upsndarr[iproc][kii];
				jj = _mvm.listupsndarr[iproc][kii];
				ibsg = _mvm.bsx[jblk]+jj;
				_mvm.vectsndrcv[nz+nzloc] = _x.vect[ibsg];
				nzloc++;
			};
			ObjSizeSend[iobj] = nzloc * sizeof (double);
			nz += nzloc;
			iobj++;
		};
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info = CMPIExchange::DataExchangeMPI (pcomm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
	if(info) throw " CGSMatrixRS::MvmA2: Error in DataExchangeMPI";

	double *darr;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		iproc = ObjIDRecv[iobj];
		darr = (double *) ObjRecv[iobj];
		nzloc = 0;
		for (kii=0;kii<_mvm.nlistuprcvarr[iproc];kii++) {
			jblk = _mvm.list2uprcvarr[iproc][kii];
			jj = _mvm.listuprcvarr[iproc][kii];
			ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
			_mvm.px.vect[ibsl] = darr[nzloc];
			nzloc ++;
		};
	};

// Free work arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Perform all local computations

	for (int ilist=0;ilist<nlist;ilist++) {
		int iblk = listb[ilist];

		ReadBlock (iblk);
		mtrarr[iblk].MvmBlockRowAu2Index2 (_mvm.bspx, _mvm.bspx2, _mvm.px, _mvm.bspx, _mvm.bspx2, _mvm.pax);
		FreeBlock (iblk);

	};

// Store the result of computations in ax

	for (kii=0;kii<_mvm.nlistuprcvarr[myid];kii++) {
		jblk = _mvm.list2uprcvarr[myid][kii];
		jj = _mvm.listuprcvarr[myid][kii];
		ibsg = _mvm.bsx[jblk]+jj;
		ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
		_ax.vect[ibsg] += _mvm.pax.vect[ibsl];
	};

#ifdef __PMPITRACE__
	MPE_Log_event (mvma2_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CGSMatrixR: General transposed secondary matrix by vector multiplication in parallel mode
//========================================================================================
void CGSMatrixR::MvmAh2Index2 (const CTree &_tree, CMvmR &_mvm, // General transposed secondary matrix by vector multiplication in parallel mode
									const CSVector &_x, CSVector &_ax) {

	const char *funcname = "MvmAh2Index2";

#ifdef __PMPITRACE__
	static int i_set = 0;
	static int mvma2h_beg;
	static int mvma2h_end;

	if (i_set == 0) {

		mvma2h_beg = MPE_Log_get_event_number ();
		mvma2h_end = MPE_Log_get_event_number ();
		MPE_Describe_state (mvma2h_beg, mvma2h_end, "MvmA2H", "red");

		i_set = 1;

	};

	MPE_Log_event (mvma2h_beg, 0, NULL);
#endif

// Remark: on entry to the routine vector _ax is assumed to be initialized

	double dzero = 0.0e0;

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Init the result

//	_ax.SetSVect (dzero);

// Initialize px and pax by zeroes

	int kii;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.px.vect[kii] = dzero;
	for (kii=0;kii<_mvm.nlocext;kii++) _mvm.pax.vect[kii] = dzero;

// Initialize px by local x data

	int jj, ibsg, ibsl, jblk;

	for (kii=0;kii<_mvm.nlistdnsndarr[myid];kii++) {
		jblk = _mvm.list2dnsndarr[myid][kii];
		jj = _mvm.listdnsndarr[myid][kii];
		ibsg = _mvm.bsx[jblk]+jj;
		ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
		_mvm.px.vect[ibsl] = _x.vect[ibsg];
	};

// Perform all local computations

	for (int ilist=0;ilist<nlist;ilist++) {
		int iblk = listb[ilist];

		ReadBlock (iblk);
		mtrarr[iblk].MvmBlockColumnAl2Index2 (true, _mvm.bspx, _mvm.bspx2, _mvm.px, _mvm.bspx, _mvm.bspx2, _mvm.pax);
		FreeBlock (iblk);

	};

// Exchange the resulting data

	int NObjSend = nproc-1;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	ObjTypeSend = new int [NObjSend];
	if(!ObjTypeSend) MemoryFail(funcname);
	ObjIDSend = new int [NObjSend];
	if(!ObjIDSend) MemoryFail(funcname);
	CpuIDSend = new int [NObjSend];
	if(!CpuIDSend) MemoryFail(funcname);
	ObjSizeSend = new int [NObjSend];
	if(!ObjSizeSend) MemoryFail(funcname);
	ObjSend = new char * [NObjSend];
	if(!ObjSend) MemoryFail(funcname);

	int iobj = 0;

	int nz = 0;
	int nzloc;

	int iproc;

	for (iproc=0;iproc<nproc;iproc++) {
		if (iproc != myid) {
			ObjTypeSend[iobj] = 1;
			ObjIDSend[iobj] = myid;
			CpuIDSend[iobj] = iproc;
			ObjSend[iobj] = (char *) (_mvm.vectsndrcv+nz);
			nzloc = 0;
			for (kii=0;kii<_mvm.nlistdnsndarr[iproc];kii++) {
				jblk = _mvm.list2dnsndarr[iproc][kii];
				jj = _mvm.listdnsndarr[iproc][kii];
				ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
				_mvm.vectsndrcv[nz+nzloc] = _mvm.pax.vect[ibsl];
				nzloc++;
			};
			ObjSizeSend[iobj] = nzloc * sizeof (double);
			nz += nzloc;
			iobj++;
		};
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	CMPIComm pcomm = _tree.GetComm ();

	int info = CMPIExchange::DataExchangeMPI (pcomm,
															NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
															ObjSizeSend, ObjSend,
															NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
															ObjSizeRecv, ObjRecv);
	if(info) throw " CGSMatrixRS::MvmA2: Error in DataExchangeMPI";

	double *darr;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		iproc = ObjIDRecv[iobj];
		darr = (double *) ObjRecv[iobj];
		nzloc = 0;
		for (kii=0;kii<_mvm.nlistdnrcvarr[iproc];kii++) {
			jblk = _mvm.list2dnrcvarr[iproc][kii];
			jj = _mvm.listdnrcvarr[iproc][kii];
			ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
			_mvm.pax.vect[ibsl] += darr[nzloc];
			nzloc++;
		};
	};

// Free work arrays

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (iobj=0;iobj<NObjRecv;iobj++) {
		delete [] ObjRecv[iobj];
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	delete [] ObjRecv;

// Store the result of computations in ax

	for (kii=0;kii<_mvm.nlistdnrcvarr[myid];kii++) {
		jblk = _mvm.list2dnrcvarr[myid][kii];
		jj = _mvm.listdnrcvarr[myid][kii];
		ibsg = _mvm.bsx[jblk]+jj;
		ibsl = _mvm.bspx[_mvm.bspx2[jblk]+jj];
		_ax.vect[ibsg] += _mvm.pax.vect[ibsl];
	};

#ifdef __PMPITRACE__
	MPE_Log_event (mvma2h_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with block column of L
//========================================================================================
void CSMatrixRS::SolveBlockColumnL (int *_sprnds, // Solve triangular system with block column of L
									int *_bsx, const CSVector &_x, 
									int *_bslx, CSVector &_lx) const { 

	const char *funcname = "SolveBlockColumnL";

// Attention: lx on entry should be initialized (by zeroes or update values, not by x data) !!!

	int i,kii,kjj,j;

	double dzero = 0.0e0;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	int ilist, blai, blaj, ibs, ibsv, jbsv, kiib, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvlx = _lx.nv;

	double *px, *plx, *pl, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			jbsv = _bslx[i];
			px = _x.vect+irhs*nvx+ibsv;
			plx = _lx.vect+irhs*nvlx+jbsv;
			for (kii=0;kii<blai;kii++) *plx++ += *px++;
		};
	};

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bslx[i];
			int i1 = ia[ilist]+1;
			int i2 = ia[ilist+1];
			j = ia[ilist];
			ibs = bsa[j];
			pablk = a+ibs;
			plx = vloc;
			for (kii=0;kii<blai;kii++) *plx++ = dzero;
			plx = vloc;
			kiib = 0;
			for (kii=0;kii<blai;kii++) {
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk+kiib;
				for (kjj=0;kjj<blai;kjj++) {
//					vloc[kii] += a[ibs+kii*blai+kjj]*_lx.vect[ibsv+kjj];
//					*plx += *pl++ * *px++;
					*plx += *pl++ * *px++;
				};
				plx++;
				kiib += blai;
			};
			ibsv = _bslx[i];
			px = vloc;
			plx = _lx.vect+irhs*nvlx+ibsv;
			for (kii=0;kii<blai;kii++) *plx++ = *px++;
			for (int j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = _sprnds[k+1]-_sprnds[k];
				jbsv = _bslx[k];
				ibs = bsa[j];
				pablk = a+ibs;
				plx = _lx.vect+irhs*nvlx+jbsv;
				kiib = 0;
				for (kii=0;kii<blaj;kii++) {
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk+kiib;
					for (kjj=0;kjj<blai;kjj++) {
//						_lx.vect[jbsv+kii] -= a[ibs+kii*blai+kjj] * _lx.vect[ibsv+kjj];
//						*plx -= *pl++ * *px++;
						*plx -= *pl++ * *px++;
					};
					plx++;
					kiib += blai;
				};
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with block column of L
//========================================================================================
void CSMatrixRS::SolvePointBlockColumnL (int *_sprnds, // Solve triangular system with block column of L
									int *_bsx, const CSVector &_x, 
									int *_bslx, CSVector &_lx) const { 

	const char *funcname = "SolvePointBlockColumnL";

// Attention: lx on entry should be initialized (by zeroes or update values, not by x data) !!!

	int i,j;

	double dzero = 0.0e0;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Direct solve

	int ilist, ibs, ibsv, jbsv, irhs;
	int i1, i2, k;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvlx = _lx.nv;

	double *px, *plx, *pl, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			ibsv = _bsx[i];
			jbsv = _bslx[i];
			px = _x.vect+irhs*nvx+ibsv;
			plx = _lx.vect+irhs*nvlx+jbsv;
			*plx += *px;
		};
	};

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0; ilist<nlist; ilist++) {
				i = list[ilist];
				ibsv = _bslx[i];
				i1 = ia[ilist]+1;
				i2 = ia[ilist+1];
				j = ia[ilist];
				ibs = bsa[j];
				pablk = a+ibs;
				plx = vloc;
				*plx = dzero;
				plx = vloc;
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk;
				*plx += *pl * *px;
				ibsv = _bslx[i];
				px = vloc;
				plx = _lx.vect+irhs*nvlx+ibsv;
				*plx = *px;
				for (j=i1; j<i2; j++) {
					k = ja[j];
					jbsv = _bslx[k];
					ibs = bsa[j];
					pablk = a+ibs;
					plx = _lx.vect+irhs*nvlx+jbsv;
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk;
					*plx -= *pl * *px;
				};
			};
		};
	} else {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			ibsv = _bslx[i];
			i1 = ia[ilist]+1;
			i2 = ia[ilist+1];
			j = ia[ilist];
			_lx.vect[ibsv] = a[j] * _lx.vect[ibsv];
			for (j=i1; j<i2; j++) {
				k = ja[j];
				jbsv = _bslx[k];
				_lx.vect[jbsv] -= a[j] * _lx.vect[ibsv];
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Solve triangular system with block column of L
//========================================================================================
void CSMatrixR::SolveBlockColumnL2Index (// Solve triangular system with block column of L
											int *_bsx, const CSVector &_x, 
											int *_bslx, CSVector &_lx) const { 

//	const char *funcname = "SolveBlockColumnL2Index";

// Attention: lx on entry should be initialized (by zeroes or update values, not by x data) !!!

	int i,j;

	double dzero = 0.0e0;

	double vloc;
	double *pvloc = &vloc;

// Direct solve

	int ilist, ibsv, jbsv, irhs;
	int i1, i2, k, iblk, jblk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvlx = _lx.nv;

	double *px, *plx, *pl, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			iblk = list2[ilist];
			ibsv = _bsx[iblk]+i;
			jbsv = _bslx[iblk]+i;
			px = _x.vect+irhs*nvx+ibsv;
			plx = _lx.vect+irhs*nvlx+jbsv;
			*plx += *px;
		};
	};

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0; ilist<nlist; ilist++) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bslx[iblk]+i;
				i1 = ia[ilist]+1;
				i2 = ia[ilist+1];
				j = ia[ilist];
				pablk = a+j;
				plx = pvloc;
				*plx = dzero;
				plx = pvloc;
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk;
				*plx += *pl * *px;
				ibsv = _bslx[iblk]+i;
				px = pvloc;
				plx = _lx.vect+irhs*nvlx+ibsv;
				*plx = *px;
				for (j=i1; j<i2; j++) {
					k = ja[j];
					jblk = ja2[j];
					jbsv = _bslx[jblk]+k;
					pablk = a+j;
					plx = _lx.vect+irhs*nvlx+jbsv;
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk;
					*plx -= *pl * *px;
				};
			};
		};
	} else {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			iblk = list2[ilist];
			ibsv = _bslx[iblk]+i;
			i1 = ia[ilist]+1;
			i2 = ia[ilist+1];
			j = ia[ilist];
			_lx.vect[ibsv] = a[j] * _lx.vect[ibsv];
			for (j=i1; j<i2; j++) {
				k = ja[j];
				jblk = ja2[j];
				jbsv = _bslx[jblk]+k;
				_lx.vect[jbsv] -= a[j] * _lx.vect[ibsv];
			};
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Solve triangular system with block column of L
//========================================================================================
void CSMatrixR::SolveBlockColumnL (// Solve triangular system with block column of L
												const CSVector &_x, CSVector &_lx) const { 

//	const char *funcname = "SolvePointBlockColumnL";

// Attention: lx on entry should be initialized (by zeroes or update values, not by x data) !!!

	int i,j;

	double dzero = 0.0e0;

	double vloc;
	double *pvloc = &vloc;

// Direct solve

	int ilist, ibs, ibsv, jbsv, irhs;
	int i1, i2, k;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvlx = _lx.nv;

	double *px, *plx, *pl, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			ibsv = i;
			jbsv = i;
			px = _x.vect+irhs*nvx+ibsv;
			plx = _lx.vect+irhs*nvlx+jbsv;
			*plx += *px;
		};
	};

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=0; ilist<nlist; ilist++) {
				i = list[ilist];
				ibsv = i;
				i1 = ia[ilist]+1;
				i2 = ia[ilist+1];
				j = ia[ilist];
				ibs = j;
				pablk = a+ibs;
				plx = pvloc;
				*plx = dzero;
				plx = pvloc;
				px = _lx.vect+irhs*nvlx+ibsv;
				pl = pablk;
				*plx += *pl * *px;
				ibsv = i;
				px = pvloc;
				plx = _lx.vect+irhs*nvlx+ibsv;
				*plx = *px;
				for (j=i1; j<i2; j++) {
					k = ja[j];
					jbsv = k;
					ibs = j;
					pablk = a+ibs;
					plx = _lx.vect+irhs*nvlx+jbsv;
					px = _lx.vect+irhs*nvlx+ibsv;
					pl = pablk;
					*plx -= *pl * *px;
				};
			};
		};
	} else {
		for (ilist=0; ilist<nlist; ilist++) {
			i = list[ilist];
			i1 = ia[ilist]+1;
			i2 = ia[ilist+1];
			j = ia[ilist];
			_lx.vect[i] = a[j] * _lx.vect[i];
			for (j=i1; j<i2; j++) {
				k = ja[j];
				_lx.vect[k] -= a[j] * _lx.vect[i];
			};
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Solve triangular system with L stored by block columns
//========================================================================================
void CGSMatrixR::SolveL (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by block columns

//	const char *funcname = "SolveL_00";

	double dzero = 0.0e0;

// Initialize lx by zeroes

	int nrhsx = _x.GetNrhs ();
	int nvlx = _lx.GetNv ();

	double *pvect = _lx.GetVect ();

	for (int kii=0;kii<nrhsx*nvlx;kii++) pvect[kii] = dzero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		mtrarr[iblk].SolveBlockColumnL (_x, _lx);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Solve triangular system with L stored by block columns
//========================================================================================
void CGSMatrixRS::SolveL (const CSVector &_x, CSVector &_lx) const { // Solve triangular system with L stored by block columns

//	const char *funcname = "SolveL_00";

	double dzero = 0.0e0;

// Initialize lx by zeroes

	int nrhsx = _x.nrhs;
	int nvlx = _lx.nv;

	for (int kii=0;kii<nrhsx*nvlx;kii++) _lx.vect[kii] = dzero;

// Cycle over the block columns 

	int iblk;

	for (iblk=0;iblk<nblksr;iblk++) {

		mtrarr[iblk].SolveBlockColumnL (sprndr, sprndr, _x, sprndr, _lx);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Solve triangular system with L stored by block columns in parallel mode
//========================================================================================
void CGSMatrixRS::SolveL (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with L stored by block columns in parallel mode
							const CSVector &_x, CSVector &_lx) {

//	const char *funcname = "SolveL_01";

	double dzero = 0.0e0;

// Determine the case

	int icase = 1;
	if (n == nsupr) icase = 0;

// Initialize px by zeroes

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = dzero;

// Search the tree

	int inode, inodetype;
	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode;

		inodetype = _tree.nodes[inode].nodetype;

		iblkbeg = _tree.nodes[inode].indbeg;
		iblkend = _tree.nodes[inode].indend;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary and store data

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (false, _tree.comm, 'A', '+', _tree.myid, childcpuid, 3*_tree.nnodes+inode,
											_mvm.px,
											nlistloc, listloc,
											blksr, sprndr, bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (false, _tree.comm, 'A', '=', _tree.myid, fathercpuid2, 3*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

				ReadBlock (iblk);
				if (icase == 1) {
					mtrarr[iblk].SolveBlockColumnL (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				} else {
					mtrarr[iblk].SolvePointBlockColumnL (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				};
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (false, _tree.comm, 'A', _tree.myid, fathercpuid, 3*_tree.nnodes+fathernode,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (false, _tree.comm, 'R', _tree.myid, childcpuid, 3*_tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

				};

			};

		};

NextNode:;

	};

// Store the result of computations in lx

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_lx.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Solve triangular system with L stored by block columns in parallel mode
//========================================================================================
void CGSMatrixR::SolveL2Index (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with L stored by block columns in parallel mode
							const CSVector &_x, CSVector &_lx) {

//	const char *funcname = "SolveL2Index";

//	ofstream fout ("ChkSlvL.dat",ios::app);

#ifdef __PMPITRACE__
	static int i_set = 0;
	static int slvl_beg;
	static int slvl_end;

	if (i_set == 0) {

		slvl_beg = MPE_Log_get_event_number ();
		slvl_end = MPE_Log_get_event_number ();
		MPE_Describe_state (slvl_beg, slvl_end, "SolveL", "white");

		i_set = 1;

	};

	MPE_Log_event (slvl_beg, 0, NULL);
#endif

	double dzero = 0.0e0;

// Initialize px by zeroes

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = dzero;

// Search the tree

	int inode, inodetype;
	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevend;ilev>=ilevbeg;ilev--) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode;

		inodetype = _tree.nodes[inode].nodetype;

		iblkbeg = _tree.nodes[inode].indbeg;
		iblkend = _tree.nodes[inode].indend;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != compnode2) {

// Receive data from childs if necessary and store data

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

				int childnode = _tree.nodes[inode].childs[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (true, _tree.comm, 'A', '+', _tree.myid, childcpuid, 3*_tree.nnodes+inode,
											_mvm.px,
											nlistloc, listloc,
											blksr, sprndr, bl2ndr);

				};

			};

		} else {

// Receive data from secondary father node and store data

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistuprcvarr[ilev];
				int *listloc = _mvm.listuprcvarr[ilev];

				_mvm.ReceiveVectData (true, _tree.comm, 'A', '=', _tree.myid, fathercpuid2, 3*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkbeg;iblk<=iblkend;iblk++) {

				ReadBlock (iblk);
				mtrarr[iblk].SolveBlockColumnL2Index (_mvm.bsx, _x, _mvm.bspx, _mvm.px);
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != exchangenode) {

// Prepare and send data to the father node if necessary

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				_mvm.SendVectData (true, _tree.comm, 'A', _tree.myid, fathercpuid, 3*_tree.nnodes+fathernode,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		} else {

// Prepare and send data to the secondary child nodes

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				int nlistloc = _mvm.nlistupsndarr[ilev];
				int *listloc = _mvm.listupsndarr[ilev];

				if (childcpuid != _tree.myid) {

					_mvm.SendVectData (true, _tree.comm, 'R', _tree.myid, childcpuid, 3*_tree.nnodes+childnode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

				};

			};

		};

NextNode:;

	};

// Store the result of computations in lx

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_lx.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

#ifdef __PMPITRACE__
	MPE_Log_event (slvl_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with the block row of U
//========================================================================================
void CSMatrixRS::SolveBlockRowU (int *_sprnds, // Solve triangular system with the block row of U
									int *_bsx, const CSVector &_x, 
									int *_bsux, CSVector &_ux) const {

	const char *funcname = "SolveBlockRowU";

	int i, kii, kjj;
//	int nloc=n;

	double dzero = 0.0e0;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	int ilist, blai, blaj, ibs, ibsv, jbsv, irhs;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	double *px, *pux, *pu, *pablk;

	for (irhs=0;irhs<nrhsx;irhs++) {
		for (ilist=nlist-1;ilist>=0;ilist--) {
			i = list[ilist];
			blai = _sprnds[i+1]-_sprnds[i];
			ibsv = _bsx[i];
			int i1 = ia[ilist]+1;
			int i2 = ia[ilist+1];
			pux = vloc;
			px = _x.vect+irhs*nvx+ibsv;
			for (kii=0;kii<blai;kii++) {
//				vloc[kii] = _x.vect[ibsv+kii];
				*pux++ = *px++;
			};
			int j;
			for (j=i1; j<i2; j++) {
				int k = ja[j];
				blaj = _sprnds[k+1]-_sprnds[k];
				jbsv = _bsux[k];
				ibs = bsa[j];
				pablk = a+ibs;
				pux = vloc;
				for (kii=0;kii<blai;kii++) {
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk+kii;
					for (kjj=0;kjj<blaj;kjj++) {
//						vloc[kii] -= a[ibs+kjj*blai+kii] * _ux.vect[jbsv+kjj];
//						*pux -= *pu * *px++;
						*pux -= *pu * *px++;
						pu += blai;
					};
					pux++;
				};
			};
			j = ia[ilist];
			ibs = bsa[j];
			pablk = a+ibs;
			ibsv = _bsux[i];
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) *pux++ = dzero;
			pux = _ux.vect+irhs*nvux+ibsv;
			for (kii=0;kii<blai;kii++) {
				px = vloc;
				pu = pablk+kii;
				for (kjj=0;kjj<blai;kjj++) {
//					_ux.vect[ibsv+kii] += a[ibs+kjj*blai+kii] * vloc[kjj];
//					*pux += *pu * *px++;
					*pux += *pu * *px++;
					pu += blai;
				};
				pux++;
			};
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixRS: Solve triangular system with the block row of U
//========================================================================================
void CSMatrixRS::SolvePointBlockRowU (int *_sprnds, // Solve triangular system with the block row of U
									int *_bsx, const CSVector &_x, 
									int *_bsux, CSVector &_ux) const {

	const char *funcname = "SolvePointBlockRowU";

	int i;
//	int nloc=n;

	double dzero = 0.0e0;

	double *vloc;

	vloc = new double [blamx];
	if (!vloc) MemoryFail (funcname);

// Backward solve

	int ilist, ibs, ibsv, jbsv, irhs;
	int i1, i2, j, k;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	double *px, *pux, *pu, *pablk;
	double aux;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=nlist-1;ilist>=0;ilist--) {
				i = list[ilist];
				ibsv = _bsx[i];
				int i1 = ia[ilist]+1;
				int i2 = ia[ilist+1];
				pux = vloc;
				px = _x.vect+irhs*nvx+ibsv;
				*pux = *px;
				for (j=i1; j<i2; j++) {
					k = ja[j];
					jbsv = _bsux[k];
					ibs = bsa[j];
					pablk = a+ibs;
					pux = vloc;
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk;
					*pux -= *pu * *px;
				};
				j = ia[ilist];
				ibs = bsa[j];
				pablk = a+ibs;
				ibsv = _bsux[i];
				pux = _ux.vect+irhs*nvux+ibsv;
				*pux = dzero;
				pux = _ux.vect+irhs*nvux+ibsv;
				px = vloc;
				pu = pablk;
				*pux += *pu * *px;
			};
		};
	} else {
		for (ilist=nlist-1;ilist>=0;ilist--) {
			i = list[ilist];
			ibsv = _bsx[i];
			i1 = ia[ilist]+1;
			i2 = ia[ilist+1];
			aux = _x.vect[ibsv];
			for (j=i1; j<i2; j++) {
				k = ja[j];
				jbsv = _bsux[k];
				aux -= a[j] * _ux.vect[jbsv];
			};
			ibsv = _bsux[i];
			j = ia[ilist];
			_ux.vect[ibsv] = aux * a[j];
		};
	};

	delete [] vloc;

};

// Author: Kharchenko S.A.
// CSMatrixR: Solve triangular system with the block row of U
//========================================================================================
void CSMatrixR::SolveBlockRowU2Index (// Solve triangular system with the block row of U
										int *_bsx, const CSVector &_x, 
										int *_bsux, CSVector &_ux) const {

//	const char *funcname = "SolveBlockRowU2Index";

	int i;
//	int nloc=n;

	double dzero = 0.0e0;

	double vloc;
	double *pvloc = &vloc;

// Backward solve

	int ilist, ibsv, jbsv, irhs;
	int i1, i2, j, k, iblk, jblk;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	double *px, *pux, *pu, *pablk;
	double aux;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=nlist-1;ilist>=0;ilist--) {
				i = list[ilist];
				iblk = list2[ilist];
				ibsv = _bsx[iblk]+i;
				int i1 = ia[ilist]+1;
				int i2 = ia[ilist+1];
				pux = pvloc;
				px = _x.vect+irhs*nvx+ibsv;
				*pux = *px;
				for (j=i1; j<i2; j++) {
					k = ja[j];
					jblk = ja2[j];
					jbsv = _bsux[jblk]+k;
					pablk = a+j;
					pux = pvloc;
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk;
					*pux -= *pu * *px;
				};
				j = ia[ilist];
				pablk = a+j;
				ibsv = _bsux[iblk]+j;
				pux = _ux.vect+irhs*nvux+ibsv;
				*pux = dzero;
				pux = _ux.vect+irhs*nvux+ibsv;
				px = pvloc;
				pu = pablk;
				*pux += *pu * *px;
			};
		};
	} else {
		for (ilist=nlist-1;ilist>=0;ilist--) {
			i = list[ilist];
			iblk = list2[ilist];
			ibsv = _bsx[iblk]+i;
			i1 = ia[ilist]+1;
			i2 = ia[ilist+1];
			aux = _x.vect[ibsv];
			for (j=i1; j<i2; j++) {
				k = ja[j];
				jblk = ja2[j];
				jbsv = _bsux[jblk]+k;
				aux -= a[j] * _ux.vect[jbsv];
			};
			j = ia[ilist];
			ibsv = _bsux[iblk]+i;
			_ux.vect[ibsv] = aux * a[j];
		};
	};

};

// Author: Kharchenko S.A.
// CSMatrixR: Solve triangular system with the block row of U
//========================================================================================
void CSMatrixR::SolveBlockRowU (// Solve triangular system with the block row of U
													const CSVector &_x, 
													CSVector &_ux) const {

//	const char *funcname = "SolveBlockRowU";

	int i;
//	int nloc=n;

	double dzero = 0.0e0;

	double vloc;
	double *pvloc = &vloc;

// Backward solve

	int ilist, ibs, ibsv, jbsv, irhs;
	int i1, i2, j, k;

	int nrhsx = _x.nrhs;
	int nvx = _x.nv;
	int nvux = _ux.nv;

	double *px, *pux, *pu, *pablk;
	double aux;

	if (nrhsx != 1) {
		for (irhs=0;irhs<nrhsx;irhs++) {
			for (ilist=nlist-1;ilist>=0;ilist--) {
				i = list[ilist];
				ibsv = i;
				int i1 = ia[ilist]+1;
				int i2 = ia[ilist+1];
				pux = pvloc;
				px = _x.vect+irhs*nvx+ibsv;
				*pux = *px;
				for (j=i1; j<i2; j++) {
					k = ja[j];
					jbsv = k;
					ibs = j;
					pablk = a+ibs;
					pux = pvloc;
					px = _ux.vect+irhs*nvux+jbsv;
					pu = pablk;
					*pux -= *pu * *px;
				};
				j = ia[ilist];
				ibs = j;
				pablk = a+ibs;
				ibsv = i;
				pux = _ux.vect+irhs*nvux+ibsv;
				*pux = dzero;
				pux = _ux.vect+irhs*nvux+ibsv;
				px = pvloc;
				pu = pablk;
				*pux += *pu * *px;
			};
		};
	} else {
		for (ilist=nlist-1;ilist>=0;ilist--) {
			i = list[ilist];
			i1 = ia[ilist]+1;
			i2 = ia[ilist+1];
			aux = _x.vect[i];
			for (j=i1; j<i2; j++) {
				k = ja[j];
				aux -= a[j] * _ux.vect[k];
			};
			j = ia[ilist];
			_ux.vect[i] = aux * a[j];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Solve triangular system with U stored by block rows
//========================================================================================
void CGSMatrixR::SolveU (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by block rows

//	const char *funcname = "SolveU_00";

// Cycle over the block rows

	int iblk;

	for (iblk=nblksr-1;iblk>=0;iblk--) {

		mtrarr[iblk].SolveBlockRowU (_x, _ux);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Solve triangular system with U stored by block rows
//========================================================================================
void CGSMatrixRS::SolveU (const CSVector &_x, CSVector &_ux) const { // Solve triangular system with U stored by block rows

//	const char *funcname = "SolveU_00";

// Cycle over the block rows

	int iblk;

	for (iblk=nblksr-1;iblk>=0;iblk--) {

		mtrarr[iblk].SolveBlockRowU (sprndr, sprndr, _x, sprndr, _ux);

	};

};

// Author: Kharchenko S.A.
// CGSMatrixRS: Solve triangular system with U stored by block rows in parallel mode
//========================================================================================
void CGSMatrixRS::SolveU (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with U stored by block rows in parallel mode
							const CSVector &_x, CSVector &_ux) {

//	const char *funcname = "SolveU_01";

	double dzero = 0.0e0;

// Determine the case

	int icase = 1;
	if (n == nsupr) icase = 0;

// Initialize px

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = dzero;

// Search the tree

	int inode, inodetype;

	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode;

		inodetype = _tree.nodes[inode].nodetype;

		iblkbeg = _tree.nodes[inode].indbeg;
		iblkend = _tree.nodes[inode].indend;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (false, _tree.comm, 'A', '=', _tree.myid, fathercpuid, 4*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (false, _tree.comm, 'R', '=', _tree.myid, childcpuid, 4*_tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

				};

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkend;iblk>=iblkbeg;iblk--) {

				ReadBlock (iblk);
				if (icase == 1) {
					mtrarr[iblk].SolveBlockRowU (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				} else {
					mtrarr[iblk].SolvePointBlockRowU (sprndr, _mvm.bsx, _x, _mvm.bspx, _mvm.px);
				};
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			if (_tree.nodes[inode].nchilds > 1) {

				for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

					int childnode = _tree.nodes[inode].childs[ichild];
					int childcpuid = _tree.nodes[childnode].nodecpu;

					int nlistloc = _mvm.nlistdnsndarr[ilev];
					int *listloc = _mvm.listdnsndarr[ilev];

					if (childcpuid != _tree.myid) {

						_mvm.SendVectData (false, _tree.comm, 'r', _tree.myid, childcpuid, 4*_tree.nnodes+childnode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

					};

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (false, _tree.comm, 'A', _tree.myid, fathercpuid2, 4*_tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		};

NextNode:;

	};

// Store the result of computations in ux

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ux.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

};

// Author: Kharchenko S.A.
// CGSMatrixR: Solve triangular system with U stored by block rows in parallel mode
//========================================================================================
void CGSMatrixR::SolveU2Index (const CTree &_tree, CMvmR &_mvm, // Solve triangular system with U stored by block rows in parallel mode
							const CSVector &_x, CSVector &_ux) {

//	const char *funcname = "SolveU2Index";

#ifdef __PMPITRACE__
	static int i_set = 0;
	static int slvu_beg;
	static int slvu_end;

	if (i_set == 0) {

		slvu_beg = MPE_Log_get_event_number ();
		slvu_end = MPE_Log_get_event_number ();
		MPE_Describe_state (slvu_beg, slvu_end, "SolveU", "violet");

		i_set = 1;

	};

	MPE_Log_event (slvu_beg, 0, NULL);
#endif

	double dzero = 0.0e0;

// Determine the case

	int icase = 1;
	if (n == nsupr) icase = 0;

// Initialize px

	for (int kii=0;kii<_mvm.nlocext*_mvm.nrhs;kii++) _mvm.px.vect[kii] = dzero;

// Search the tree

	int inode, inodetype;

	int iblk, iblkbeg, iblkend, fathernode, fathercpuid, fathernode2, fathercpuid2;

	int nodebeg = _tree.cpuidbeg[_tree.myid];
	int nodeend = _tree.cpuidend[_tree.myid];

	int ilevbeg = _tree.nodes[nodebeg].nodelv;
	int ilevend = _tree.nodes[nodeend].nodelv;

	for (int ilev=ilevbeg;ilev<=ilevend;ilev++) {

		inode = _mvm.lev2node[ilev];

		int iproc = _tree.nodes[inode].nodecpu;

		if (iproc != _tree.myid) goto NextNode;

		inodetype = _tree.nodes[inode].nodetype;

		iblkbeg = _tree.nodes[inode].indbeg;
		iblkend = _tree.nodes[inode].indend;

		fathernode = _tree.nodes[inode].fatherid;
		fathercpuid = _tree.nodes[inode].fathercpu;
		fathernode2 = _tree.nodes[inode].fatherid2;
		fathercpuid2 = _tree.nodes[inode].fathercpu2;

// Switch the cases

		if (inodetype != exchangenode) {

// Receive data from father node if necessary and store them

			if (fathercpuid != _tree.myid) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				_mvm.ReceiveVectData (true, _tree.comm, 'A', '=', _tree.myid, fathercpuid, 4*_tree.nnodes+inode,
									_mvm.px,
									nlistloc, listloc,
									blksr, sprndr, bl2ndr);

			};

		} else {

// Receive data from child nodes if necessary and store them

			for (int ichild=0;ichild<_tree.nodes[inode].nchilds2;ichild++) {

				int nlistloc = _mvm.nlistdnrcvarr[ilev];
				int *listloc = _mvm.listdnrcvarr[ilev];

				int childnode = _tree.nodes[inode].childs2[ichild];
				int childcpuid = _tree.nodes[childnode].nodecpu;

				if (childcpuid != _tree.myid) {

					_mvm.ReceiveVectData (true, _tree.comm, 'R', '=', _tree.myid, childcpuid, 4*_tree.nnodes+inode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

				};

			};

		};

// Perform local computations

		if (inodetype != exchangenode) {

			for (iblk=iblkend;iblk>=iblkbeg;iblk--) {

				ReadBlock (iblk);
				mtrarr[iblk].SolveBlockRowU2Index (_mvm.bsx, _x, _mvm.bspx, _mvm.px);
				FreeBlock (iblk);

			};

		};

// Switch the cases

		if (inodetype != compnode2) {

// Prepare and send data to the child nodes if necessary

			if (_tree.nodes[inode].nchilds > 1) {

				for (int ichild=0;ichild<_tree.nodes[inode].nchilds;ichild++) {

					int childnode = _tree.nodes[inode].childs[ichild];
					int childcpuid = _tree.nodes[childnode].nodecpu;

					int nlistloc = _mvm.nlistdnsndarr[ilev];
					int *listloc = _mvm.listdnsndarr[ilev];

					if (childcpuid != _tree.myid) {

						_mvm.SendVectData (true, _tree.comm, 'r', _tree.myid, childcpuid, 4*_tree.nnodes+childnode,
										_mvm.px,
										nlistloc, listloc,
										blksr, sprndr, bl2ndr);

					};

				};

			};

		} else {

// Prepare and send data to the secondary father node

			if (fathercpuid2 != _tree.myid) {

				int nlistloc = _mvm.nlistdnsndarr[ilev];
				int *listloc = _mvm.listdnsndarr[ilev];

				_mvm.SendVectData (true, _tree.comm, 'A', _tree.myid, fathercpuid2, 4*_tree.nnodes+fathernode2,
								_mvm.px,
								nlistloc, listloc,
								blksr, sprndr, bl2ndr);

			};

		};

NextNode:;

	};

// Store the result of computations in ux

	for (int irhs=0;irhs<_mvm.nrhs;irhs++) {
		for (int kii=0;kii<_mvm.nloc;kii++) {
			_ux.vect[irhs*_mvm.nloc+kii] = _mvm.px.vect[irhs*_mvm.nlocext+kii];
		};
	};

#ifdef __PMPITRACE__
	MPE_Log_event (slvu_end, 0, NULL);
#endif

};

// Author: Kharchenko S.A.
// CMvmSchur: Implement abstract block Mvm with Schur complement parallelization
//========================================================================================
void CMvmSchur::AbstractBlockMvmA (char _transp, // Implement abstract block Mvm with Schur complement parallelization
												CSVector &_x, CSVector &_ax) {

//	const char *funcname = "AbstractBlockMvmA";

// Exchange necessary x data

	ExchangeX (_x);

// Init necessary Ax by zeroes

	InitAxByZeroes ();

// Multiply

	PerformLocalMultiplications (_transp);

// Exchange ax and store

	ExchangeAX (_ax);

};

// Author: Kharchenko S.A.
// CMvmSchur: Implement abstract block L solve with Schur complement parallelization
//========================================================================================
void CMvmSchur::AbstractBlockSolveL (char _lutype, // Implement abstract block L solve with Schur complement parallelization
													CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
													CSVector &_x, CSVector &_px) {

	const char *funcname = "AbstractBlockSolveL";

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Open the tree structure

	int nlevtot = _tree.GetNlev ();
//	int nnodestot = _tree.GetNnodes ();
	int *pcpuidend = _tree.GetCpuidend ();
	CNode *pnodes = _tree.GetNode ();

// Open the maximal factorization block sparsity structure

	int *piablksp = _ablksp.GetIa ();
	int *pjablksp = _ablksp.GetJa ();

// Get description of the tree for the current cpu

	int *i_lev_arr;

	i_lev_arr = new int [nlevtot];
	if (!i_lev_arr) MemoryFail (funcname);

	int ilev, ifather, inodeloc;

	inodeloc = pcpuidend[myid];

	while (true) {
		ilev = pnodes[inodeloc].GetNodelv ();
		i_lev_arr[ilev] = inodeloc;
		ifather = pnodes[inodeloc].FatherId ();
		if (ifather == inodeloc) break;
		inodeloc = ifather;
	};

// Allocate the work arrays

	int *iagroup;
	int *grp2cpu;
	int *grptype;
	int *iagrp2cpu;
	int *jagrp2cpu;
	int *imaskcpu;
	int *listcpu;
	int *iaset2grp;
	int *listfct;
	int *listschur;
	int *imaskschur;

	iagroup = new int [_nblks+1];
	if (!iagroup) MemoryFail (funcname);
	grp2cpu = new int [_nblks];
	if (!grp2cpu) MemoryFail (funcname);
	grptype = new int [_nblks];
	if (!grptype) MemoryFail (funcname);
	iagrp2cpu = new int [_nblks+1];
	if (!iagrp2cpu) MemoryFail (funcname);
	jagrp2cpu = new int [_nblks*nproc];
	if (!jagrp2cpu) MemoryFail (funcname);
	imaskcpu = new int [nproc];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);
	iaset2grp = new int [_nblks+1];
	if (!iaset2grp) MemoryFail (funcname);
	listfct = new int [_nblks];
	if (!listfct) MemoryFail (funcname);
	listschur = new int [_nblks];
	if (!listschur) MemoryFail (funcname);
	imaskschur = new int [_nblks];
	if (!imaskschur) MemoryFail (funcname);

	int icyclecpu = -1;

	int i;

	for (i=0;i<nproc;i++) imaskcpu[i] = icyclecpu;

	int icycleschur = -1;

	for (i=0;i<_nblks;i++) imaskschur[i] = icycleschur;

// Init Schur blocks

	int j, iblk, jjblk;

	icycleschur++;

	int nlistschur = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid) {
			if (imaskschur[iblk] != icycleschur) {
				listschur[nlistschur] = iblk;
				nlistschur++;
				imaskschur[iblk] = icycleschur;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskschur[jjblk] != icycleschur) {
						listschur[nlistschur] = jjblk;
						nlistschur++;
						imaskschur[jjblk] = icycleschur;
					};
				};
			};
		};
	};

// Init Schur blocks

	InitPxBlocks (_x);

// Main cycle over the nodes

	int inode, ibegblk, iendblk, ngroups;
	int iproc, jproc, jbegblk, jendblk, icheck, nz, nlistcpu, nsets, igroup;
	int jgroup, jjbegblk, jjendblk, iset, ibeggrp, iendgrp, myidloc;
	int nchildsloc, ichild0, ichild1, iproc0, iproc1, nlistfct;
	int *pchildsloc;

	for (ilev=nlevtot-1;ilev>=0;ilev--) {

		inode = i_lev_arr[ilev];

		ibegblk = pnodes[inode].GetIndbeg ();
		iendblk = pnodes[inode].GetIndend ();

// Create the whole set of continuous groups of blocks (set of consecutive blocks on the same cpu)

		iblk = ibegblk;
		iproc = _blk2cpu[iblk];
		ngroups = 0;
		iagroup[0] = 0;
		grp2cpu[ngroups] = iproc;
		iagroup[ngroups+1] = iblk+1-ibegblk;
		while (iblk < iendblk) {
			iblk++;
			jproc = _blk2cpu[iblk];
			if (jproc == iproc) {
				iagroup[ngroups+1] = iblk+1-ibegblk;
			} else {
				ngroups++;
				iproc = jproc;
				grp2cpu[ngroups] = iproc;
				iagroup[ngroups+1] = iblk+1-ibegblk;
			};
		};
		ngroups++;

// Mark all groups (not used, Schur or compute)

		for (i=0;i<ngroups;i++) {
			iproc = grp2cpu[i];
			if (iproc == myid) {
				grptype[i] = 1;
			} else {
				jbegblk = iagroup[i]+ibegblk;
				jendblk = iagroup[i+1]-1+ibegblk;
				icheck = -1;
				for (iblk=jbegblk;iblk<=jendblk;iblk++) {
					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk <= iblk) {
							if (_blk2cpu[jjblk] == myid) icheck = 1;
						};
					};
				};
				if (icheck > 0) {
					grptype[i] = 0;
				} else {
					grptype[i] = -1;
				};
			};
		};

// For each working group determine the list of cpu's

		iagrp2cpu[0] = 0;
		nz = 0;

		for (i=0;i<ngroups;i++) {
			icyclecpu++;
			jbegblk = iagroup[i]+ibegblk;
			jendblk = iagroup[i+1]-1+ibegblk;
			iproc = grp2cpu[i];
			listcpu[0] = iproc;
			imaskcpu[iproc] = icyclecpu;
			nlistcpu = 1;
			for (iblk=jbegblk;iblk<=jendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk <= iblk) {
						jproc = _blk2cpu[jjblk];
						if (imaskcpu[jproc] != icyclecpu) {
							listcpu[nlistcpu] = jproc;
							nlistcpu++;
							imaskcpu[jproc] = icyclecpu;
						};
					};
				};
			};
			if (nlistcpu > 1) qsort (listcpu+1,nlistcpu-1,sizeof(int),compint);
			for (j=0;j<nlistcpu;j++) jagrp2cpu[nz+j] = listcpu[j];
			nz += nlistcpu;
			iagrp2cpu[i+1] = nz;
		};

// Find the sets of sequential groups that should be processed in parallel

		nsets = 0;

		igroup = 0;
		iaset2grp[0] = 0;
		iaset2grp[nsets+1] = igroup+1;
		jbegblk = iagroup[igroup]+ibegblk;
		jendblk = iagroup[igroup+1]-1+ibegblk;
		while (igroup < ngroups-1) {
			jgroup = igroup+1;
			jjbegblk = iagroup[jgroup]+ibegblk;
			jjendblk = iagroup[jgroup+1]-1+ibegblk;
			icheck = 1;
			for (iblk=jjbegblk;iblk<=jjendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk >= jbegblk && jjblk <= jendblk) icheck = -1;
				};
			};
			if (icheck == 1) {
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			} else {
				nsets++;
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jbegblk = iagroup[igroup]+ibegblk;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			};
		};
		nsets++;

// Main cycle over the sets

		for (iset=0;iset<nsets;iset++) {

			ibeggrp = iaset2grp[iset];
			iendgrp = iaset2grp[iset+1]-1;

// Make binary tree type Schur exchanges and summation for all groups of blocks in the set

			for (igroup=ibeggrp;igroup<=iendgrp;igroup++) {
				if (grptype[igroup] >= 0) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistschur = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listschur[nlistschur] = i;
						nlistschur++;
					};

// Plan sends and receives

					myidloc = -1;

					nlistcpu = 0;
					for (j=iagrp2cpu[igroup];j<iagrp2cpu[igroup+1];j++) {
						jproc = jagrp2cpu[j];
						if (jproc == myid) myidloc = nlistcpu;
						listcpu[nlistcpu] = jproc;
						nlistcpu++;
					};

// Do exchanges and summation according to the plan

					inodeloc = pcpuidend[myidloc];

					while (true) {
						nchildsloc = pnodes[inodeloc].GetNchilds ();
						pchildsloc = pnodes[inodeloc].GetChilds ();
						iproc = pnodes[inodeloc].GetNodecpu ();
						if (iproc != myidloc) break;
						if (nchildsloc > 1) {
							ichild0 = pchildsloc[0];
							ichild1 = pchildsloc[1];
							iproc0 = pnodes[ichild0].GetNodecpu ();
							iproc1 = pnodes[ichild1].GetNodecpu ();
							if (iproc0 == myidloc) {
								if (iproc1 < nlistcpu) {
									jproc = listcpu[iproc1];
									ReceiveBlocks (true, jproc, nlistschur, listschur);
								};
							} else if (iproc1 == myidloc) {
								if (iproc0 < nlistcpu) {
									jproc = listcpu[iproc0];
									ReceiveBlocks (true, jproc, nlistschur, listschur);
								};
							};
						};
						ifather = pnodes[inodeloc].FatherId ();
						if (ifather == inodeloc) break;
						iproc0 = pnodes[ifather].GetNodecpu ();
						if (iproc0 != myidloc) {
							jproc = listcpu[iproc0];
							SendBlocks (1, &jproc, nlistschur, listschur);
						};
						inodeloc = ifather;
					};

				};
			};

// Wait for completion of sends

			WaitSends ();

// Perform solution for own blocks

			for (igroup=ibeggrp;igroup<=iendgrp;igroup++) {

				if (grptype[igroup] == 1) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistfct = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listfct[nlistfct] = i;
						nlistfct++;
					};

					SolveLBlocks (_lutype, nlistfct, listfct);

				};

			};

		};

	};

	StorePx (_px);

// Free work arrays

	delete [] i_lev_arr;
	delete [] iagroup;
	delete [] grp2cpu;
	delete [] grptype;
	delete [] iagrp2cpu;
	delete [] jagrp2cpu;
	delete [] imaskcpu;
	delete [] listcpu;
	delete [] iaset2grp;
	delete [] listfct;
	delete [] listschur;
	delete [] imaskschur;

};

// Author: Kharchenko S.A.
// CMvmSchur: Implement abstract block U solve with Schur complement parallelization
//========================================================================================
void CMvmSchur::AbstractBlockSolveU (char _lutype, // Implement abstract block U solve with Schur complement parallelization
													CTree &_tree, CSMatrix &_ablksp, int _nblks, int *_blk2cpu,
													CSVector &_x, CSVector &_px) {

	const char *funcname = "AbstractBlockSolveU";

	int myid = _tree.GetMyid ();
	int nproc = _tree.GetNproc ();

// Open the tree structure

	int nlevtot = _tree.GetNlev ();
//	int nnodestot = _tree.GetNnodes ();
	int *pcpuidend = _tree.GetCpuidend ();
	CNode *pnodes = _tree.GetNode ();

// Open the maximal factorization block sparsity structure

	int *piablksp = _ablksp.GetIa ();
	int *pjablksp = _ablksp.GetJa ();

// Get description of the tree for the current cpu

	int *i_lev_arr;

	i_lev_arr = new int [nlevtot];
	if (!i_lev_arr) MemoryFail (funcname);

	int ilev, ifather, inodeloc;

	inodeloc = pcpuidend[myid];

	while (true) {
		ilev = pnodes[inodeloc].GetNodelv ();
		i_lev_arr[ilev] = inodeloc;
		ifather = pnodes[inodeloc].FatherId ();
		if (ifather == inodeloc) break;
		inodeloc = ifather;
	};

// Allocate the work arrays

	int *iagroup;
	int *grp2cpu;
	int *grptype;
	int *iagrp2cpu;
	int *jagrp2cpu;
	int *imaskcpu;
	int *listcpu;
	int *iaset2grp;
	int *listfct;
	int *listschur;
	int *imaskschur;

	iagroup = new int [_nblks+1];
	if (!iagroup) MemoryFail (funcname);
	grp2cpu = new int [_nblks];
	if (!grp2cpu) MemoryFail (funcname);
	grptype = new int [_nblks];
	if (!grptype) MemoryFail (funcname);
	iagrp2cpu = new int [_nblks+1];
	if (!iagrp2cpu) MemoryFail (funcname);
	jagrp2cpu = new int [_nblks*nproc];
	if (!jagrp2cpu) MemoryFail (funcname);
	imaskcpu = new int [nproc];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);
	iaset2grp = new int [_nblks+1];
	if (!iaset2grp) MemoryFail (funcname);
	listfct = new int [_nblks];
	if (!listfct) MemoryFail (funcname);
	listschur = new int [_nblks];
	if (!listschur) MemoryFail (funcname);
	imaskschur = new int [_nblks];
	if (!imaskschur) MemoryFail (funcname);

	int icyclecpu = -1;

	int i;

	for (i=0;i<nproc;i++) imaskcpu[i] = icyclecpu;

	int icycleschur = -1;

	for (i=0;i<_nblks;i++) imaskschur[i] = icycleschur;

// Init Schur blocks

	int j, iblk, jjblk;

	icycleschur++;

	int nlistschur = 0;

	for (iblk=0;iblk<_nblks;iblk++) {
		if (_blk2cpu[iblk] == myid) {
			if (imaskschur[iblk] != icycleschur) {
				listschur[nlistschur] = iblk;
				nlistschur++;
				imaskschur[iblk] = icycleschur;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskschur[jjblk] != icycleschur) {
						listschur[nlistschur] = jjblk;
						nlistschur++;
						imaskschur[jjblk] = icycleschur;
					};
				};
			};
		};
	};

// Init Schur blocks

	InitPxBlocks (_x);

// Main cycle over the nodes

	int inode, ibegblk, iendblk, ngroups;
	int iproc, jproc, jbegblk, jendblk, icheck, nz, nlistcpu, nsets, igroup;
	int jgroup, jjbegblk, jjendblk, iset, ibeggrp, iendgrp;
	int nlistfct;

	for (ilev=0;ilev<nlevtot;ilev++) {

		inode = i_lev_arr[ilev];

		ibegblk = pnodes[inode].GetIndbeg ();
		iendblk = pnodes[inode].GetIndend ();

// Create the whole set of continuous groups of blocks (set of consecutive blocks on the same cpu)

		iblk = ibegblk;
		iproc = _blk2cpu[iblk];
		ngroups = 0;
		iagroup[0] = 0;
		grp2cpu[ngroups] = iproc;
		iagroup[ngroups+1] = iblk+1-ibegblk;
		while (iblk < iendblk) {
			iblk++;
			jproc = _blk2cpu[iblk];
			if (jproc == iproc) {
				iagroup[ngroups+1] = iblk+1-ibegblk;
			} else {
				ngroups++;
				iproc = jproc;
				grp2cpu[ngroups] = iproc;
				iagroup[ngroups+1] = iblk+1-ibegblk;
			};
		};
		ngroups++;

// Mark all groups (not used, Schur or compute)

		for (i=0;i<ngroups;i++) {
			iproc = grp2cpu[i];
			if (iproc == myid) {
				grptype[i] = 1;
			} else {
				jbegblk = iagroup[i]+ibegblk;
				jendblk = iagroup[i+1]-1+ibegblk;
				icheck = -1;
				for (iblk=jbegblk;iblk<=jendblk;iblk++) {
					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk <= iblk) {
							if (_blk2cpu[jjblk] == myid) icheck = 1;
						};
					};
				};
				if (icheck > 0) {
					grptype[i] = 0;
				} else {
					grptype[i] = -1;
				};
			};
		};

// For each working group determine the list of cpu's

		iagrp2cpu[0] = 0;
		nz = 0;

		for (i=0;i<ngroups;i++) {
			icyclecpu++;
			jbegblk = iagroup[i]+ibegblk;
			jendblk = iagroup[i+1]-1+ibegblk;
			iproc = grp2cpu[i];
			listcpu[0] = iproc;
			imaskcpu[iproc] = icyclecpu;
			nlistcpu = 1;
			for (iblk=jbegblk;iblk<=jendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk <= iblk) {
						jproc = _blk2cpu[jjblk];
						if (imaskcpu[jproc] != icyclecpu) {
							listcpu[nlistcpu] = jproc;
							nlistcpu++;
							imaskcpu[jproc] = icyclecpu;
						};
					};
				};
			};
			if (nlistcpu > 1) qsort (listcpu+1,nlistcpu-1,sizeof(int),compint);
			for (j=0;j<nlistcpu;j++) jagrp2cpu[nz+j] = listcpu[j];
			nz += nlistcpu;
			iagrp2cpu[i+1] = nz;
		};

// Find the sets of sequential groups that should be processed in parallel

		nsets = 0;

		igroup = 0;
		iaset2grp[0] = 0;
		iaset2grp[nsets+1] = igroup+1;
		jbegblk = iagroup[igroup]+ibegblk;
		jendblk = iagroup[igroup+1]-1+ibegblk;
		while (igroup < ngroups-1) {
			jgroup = igroup+1;
			jjbegblk = iagroup[jgroup]+ibegblk;
			jjendblk = iagroup[jgroup+1]-1+ibegblk;
			icheck = 1;
			for (iblk=jjbegblk;iblk<=jjendblk;iblk++) {
				for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
					jjblk = pjablksp[j];
					if (jjblk >= jbegblk && jjblk <= jendblk) icheck = -1;
				};
			};
			if (icheck == 1) {
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			} else {
				nsets++;
				igroup = jgroup;
				iaset2grp[nsets+1] = igroup+1;
				jbegblk = iagroup[igroup]+ibegblk;
				jendblk = iagroup[igroup+1]-1+ibegblk;
			};
		};
		nsets++;

// Main cycle over the sets

		for (iset=nsets-1;iset>=0;iset--) {

			ibeggrp = iaset2grp[iset];
			iendgrp = iaset2grp[iset+1]-1;

// Perform solution for own blocks

			for (igroup=iendgrp;igroup>=ibeggrp;igroup--) {

				if (grptype[igroup] == 1) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistfct = 0;
					for (i=jendblk;i>=jbegblk;i--) {
						listfct[nlistfct] = i;
						nlistfct++;
					};

					SolveUBlocks (_lutype, nlistfct, listfct);

				};

			};

// Make solution exchanges for all groups of blocks in the set

			for (igroup=iendgrp;igroup>=ibeggrp;igroup--) {
				if (grptype[igroup] >= 0) {

					jbegblk = iagroup[igroup]+ibegblk;
					jendblk = iagroup[igroup+1]-1+ibegblk;

					nlistschur = 0;
					for (i=jbegblk;i<=jendblk;i++) {
						listschur[nlistschur] = i;
						nlistschur++;
					};

// Plan sends and receives

					nlistcpu = 0;
					for (j=iagrp2cpu[igroup];j<iagrp2cpu[igroup+1];j++) {
						jproc = jagrp2cpu[j];
						listcpu[nlistcpu] = jproc;
						nlistcpu++;
					};

					if (nlistcpu > 1) {

						if (grptype[igroup] == 1) {
							SendBlocks (nlistcpu-1, listcpu+1, nlistschur, listschur);
						} else {
							ReceiveBlocks (false, listcpu[0], nlistschur, listschur);
						};

					};
				};
			};

// Wait for completion of sends

			WaitSends ();

		};

	};

	StorePx (_px);

// Free work arrays

	delete [] i_lev_arr;
	delete [] iagroup;
	delete [] grp2cpu;
	delete [] grptype;
	delete [] iagrp2cpu;
	delete [] jagrp2cpu;
	delete [] imaskcpu;
	delete [] listcpu;
	delete [] iaset2grp;
	delete [] listfct;
	delete [] listschur;
	delete [] imaskschur;

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CMvmSchurR::CMvmSchurR()
//========================================================================================
CMvmSchurR::CMvmSchurR (): CMvmSchur () { // Memory allocation zero data constructor

	const char *funcname = "CMvmSchurR_00";

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [0];
	if (!ibsvect) MemoryFail (funcname);
	listblk = new int [0];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [0];
	if (!listblkext) MemoryFail (funcname);
	iasndblk = new int [0];
	if (!iasndblk) MemoryFail (funcname);
	jasndblk = new int [0];
	if (!jasndblk) MemoryFail (funcname);
	iarcvblk = new int [0];
	if (!iarcvblk) MemoryFail (funcname);
	jarcvblk = new int [0];
	if (!jarcvblk) MemoryFail (funcname);

	nlistblk = 0;
	nlistblkext = 0;

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	pvect = new CSVector;
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVector;
	if (!pvect1) MemoryFail (funcname);

	PrepareWorkData ();

};

// Author: Kharchenko S.A.
// Description: Memory allocation zero data constructor
// CMvmSchurR::CMvmSchurR()
//========================================================================================
CMvmSchurR::CMvmSchurR (std::ofstream *_pfout, // Memory allocation zero data constructor
								int _nblks, int *_blks, int *_blk2cpu,
								CTree *_ptree, CSMatrix *_pablkstr, 
								void *_pgmtral, void *_pgmtrau, 
								void *_pgmtrl, void *_pgmtru,
								FUNC_MULT *_pmvmal, FUNC_MULT *_pmvmau,
								FUNC_SOLVE *_psolvel, FUNC_MULT *_psolveu,
								CSlvParam *_pparam
								): CMvmSchur () {

	const char *funcname = "CMvmSchurR_01";

	pfout = _pfout;
	nblks = _nblks;
	pblks = _blks;
	pblk2cpu = _blk2cpu;
	ptree = _ptree;
	pablkstr = _pablkstr;
	pgmtral = _pgmtral;
	pgmtrau = _pgmtrau;
	pgmtrl = _pgmtrl;
	pgmtru = _pgmtru;
	pmvmal = _pmvmal;
	pmvmau = _pmvmau;
	psolvel = _psolvel;
	psolveu = _psolveu;

	pparam = _pparam;

	icycleflt = -1;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [0];
	if (!ibsvect) MemoryFail (funcname);
	listblk = new int [0];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [0];
	if (!listblkext) MemoryFail (funcname);
	iasndblk = new int [0];
	if (!iasndblk) MemoryFail (funcname);
	jasndblk = new int [0];
	if (!jasndblk) MemoryFail (funcname);
	iarcvblk = new int [0];
	if (!iarcvblk) MemoryFail (funcname);
	jarcvblk = new int [0];
	if (!jarcvblk) MemoryFail (funcname);

	nlistblk = 0;
	nlistblkext = 0;

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	pvect = new CSVector;
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVector;
	if (!pvect1) MemoryFail (funcname);

	PrepareWorkData ();

};

// Author: Kharchenko S.A.
// Description: Destructor
// CMvmSchurR::~CMvmSchurR()
//========================================================================================
CMvmSchurR::~CMvmSchurR () { // Destructor
//	cout << " On entry to CMvmSchurR destructor " << endl;

	DestroyWorkData ();

	delete [] imaskflt;
	delete [] ibsvect;
	delete [] listblk;
	delete [] listblkext;
	delete [] iasndblk;
	delete [] jasndblk;
	delete [] iarcvblk;
	delete [] jarcvblk;

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

	delete pvect;
	delete pvect1;

//	cout << " On return from CMvmSchurR destructor " << endl;
};

// Author: Kharchenko S.A.
// Description: Prepare solve working data
// CMvmSchurR::PrepareWorkData()
//========================================================================================
void CMvmSchurR::PrepareWorkData () { // Prepare solve working data

	const char *funcname = "PrepareWorkData";

// Open the tree

	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Open the sparsity

	int *piablksp = pablkstr->GetIa ();
	int *pjablksp = pablkstr->GetJa ();

// Init the mask array

	delete [] imaskflt;
	delete [] ibsvect;

	imaskflt = new int [nblks];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [nblks];
	if (!ibsvect) MemoryFail (funcname);

	int i;
	icycleflt = -1;

	for (i=0;i<nblks;i++) imaskflt[i] = icycleflt;
	for (i=0;i<nblks;i++) ibsvect[i] = -1;

// Compute the maximal list of blocks to be stored on current cpu

	int *listloc;

	listloc = new int [nblks];
	if (!listloc) MemoryFail (funcname);

	int j, iblk, jjblk, k, kkblk;

	icycleflt++;

	int nlistschur = 0;
	nlistblk = 0;

	for (iblk=0;iblk<nblks;iblk++) {
		if (pblk2cpu[iblk] == myid) {
			nlistblk++;
			if (imaskflt[iblk] != icycleflt) {
				listloc[nlistschur] = iblk;
				nlistschur++;
				imaskflt[iblk] = icycleflt;
			};
			for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
				jjblk = pjablksp[j];
				if (jjblk >= iblk) {
					if (imaskflt[jjblk] != icycleflt) {
						listloc[nlistschur] = jjblk;
						nlistschur++;
						imaskflt[jjblk] = icycleflt;
					};
				} else if (jjblk < iblk) {
					for (k=piablksp[jjblk];k<piablksp[jjblk+1];k++) {
						kkblk = pjablksp[k];
						if (kkblk >= iblk) {
							if (imaskflt[kkblk] != icycleflt) {
								listloc[nlistschur] = kkblk;
								nlistschur++;
								imaskflt[kkblk] = icycleflt;
							};
						};
					};
				};
			};
		};
	};

// Add tree data

	int *pcpuidend = ptree->GetCpuidend ();
	CNode *pnodes = ptree->GetNode ();

	int inode = pcpuidend[myid];

	int indbeg, indend, ifather;

	while (true) {
		indbeg = pnodes[inode].GetIndbeg ();
		indend = pnodes[inode].GetIndend ();
		for (i=indbeg;i<=indend;i++) {
			if (imaskflt[i] != icycleflt) {
				listloc[nlistschur] = i;
				nlistschur++;
				imaskflt[i] = icycleflt;
			};
		};
		ifather = pnodes[inode].FatherId ();
		if (ifather == inode) break;
		inode = ifather;
	};

// Sort the list

	qsort (listloc, nlistschur, sizeof(int), compint);

	delete [] listblk;
	delete [] listblkext;

	nlistblkext = nlistschur;

	listblk = new int [nlistblk];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [nlistblkext];
	if (!listblkext) MemoryFail (funcname);

	for (i=0;i<nlistblkext;i++) listblkext[i] = listloc[i];

	nlistblk = 0;

	for (iblk=0;iblk<nblks;iblk++) {
		if (pblk2cpu[iblk] == myid) {
			listblk[nlistblk] = iblk;
			nlistblk++;
		};
	};

// Create sends/receives list

	delete [] iasndblk;
	delete [] jasndblk;

	iasndblk = new int [nproc+1];
	if (!iasndblk) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) {
		iasndblk[i] = 0;
	};

	int iproc, jproc, ni;

	for (iproc=0;iproc<nproc;iproc++) {

		if (iproc != myid) {

			icycleflt++;

			for (iblk=0;iblk<nblks;iblk++) {
				if (pblk2cpu[iblk] == iproc) {

					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
							jproc = pblk2cpu[jjblk];
							if (jproc == myid) {
								iasndblk[iproc+1]++;
							};
							imaskflt[jjblk] = icycleflt;
						};
					};

				};
			};

		};

	};

	for (i=0;i<nproc;i++) iasndblk[i+1] = iasndblk[i]+iasndblk[i+1];

	int nzjasnd = iasndblk[nproc];

	jasndblk = new int [nzjasnd];
	if (!jasndblk) MemoryFail (funcname);

	int *iptr;

	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) iptr[i] = iasndblk[i];

	for (iproc=0;iproc<nproc;iproc++) {

		if (iproc != myid) {

			icycleflt++;

			for (iblk=0;iblk<nblks;iblk++) {
				if (pblk2cpu[iblk] == iproc) {

					for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
						jjblk = pjablksp[j];
						if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
							jproc = pblk2cpu[jjblk];
							if (jproc == myid) {
								k = iptr[iproc];
								jasndblk[k] = jjblk;
								iptr[iproc]++;
							};
							imaskflt[jjblk] = icycleflt;
						};
					};

				};
			};

		};
	};

	for (i=0;i<nproc;i++) {
		ni = iasndblk[i+1]-iasndblk[i];
		if (ni != 0) qsort (jasndblk+iasndblk[i], ni, sizeof(int), compint);
	};

	delete [] iarcvblk;
	delete [] jarcvblk;

	iarcvblk = new int [nproc+1];
	if (!iarcvblk) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) {
		iarcvblk[i] = 0;
	};

	icycleflt++;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
			jjblk = pjablksp[j];
			if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
				jproc = pblk2cpu[jjblk];
				if (jproc != myid) {
					iarcvblk[jproc+1]++;
				};
				imaskflt[jjblk] = icycleflt;
			};
		};
	};

	for (i=0;i<nproc;i++) iarcvblk[i+1] = iarcvblk[i]+iarcvblk[i+1];

	int nzjarcv= iarcvblk[nproc];

	jarcvblk = new int [nzjarcv];
	if (!jarcvblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) iptr[i] = iarcvblk[i];

	icycleflt++;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		for (j=piablksp[iblk];j<piablksp[iblk+1];j++) {
			jjblk = pjablksp[j];
			if (jjblk >= iblk && imaskflt[jjblk] != icycleflt) {
				jproc = pblk2cpu[jjblk];
				if (jproc != myid) {
					k = iptr[jproc];
					jarcvblk[k] = jjblk;
					iptr[jproc]++;
				};
				imaskflt[jjblk] = icycleflt;
			};
		};
	};

	for (i=0;i<nproc;i++) {
		ni = iarcvblk[i+1]-iarcvblk[i];
		if (ni != 0) qsort (jarcvblk+iarcvblk[i], ni, sizeof(int), compint);
	};

// Create the maximal vector work array

	int nloc = 0;

	for (i=0;i<nlistschur;i++) {
		iblk = listloc[i];
		ni = pblks[iblk+1]-pblks[iblk];
		ibsvect[iblk] = nloc;
		nloc += ni;
	};

	delete pvect;
	delete pvect1;

	pvect = new CSVector (nloc);
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVector (nloc);
	if (!pvect1) MemoryFail (funcname);

// Allocate send receive work arrays

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

// Allocate send/receive structures

	int nblkmax = 0;

	int nprocloc = ptree->GetNproc ();
	int nnodesloc = ptree->GetNnodes ();

	int iblkbeg, iblkend;

	for (i=0;i<nnodesloc;i++) {
		iblkbeg = pnodes[i].GetIndbeg ();
		iblkend = pnodes[i].GetIndend ();
		ni = iblkend-iblkbeg+1;
		if (ni>nblkmax) nblkmax = ni;
	};

	nsendsmax = nprocloc*nblkmax+10;
	nsends = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	for (i=0;i<nsendsmax;i++) imasksend[i] = -1;

// Free work arrays

	delete [] listloc;
	delete [] iptr;

};

// Author: Kharchenko S.A.
// Description: Exchange X data
// CMvmSchurR::ExchangeX()
//========================================================================================
void CMvmSchurR::ExchangeX (CSVector &_x) { // Exchange X data

	const char *funcname = "ExchangeX";

//	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Init X part

	double *px = _x.GetVect ();
	double *pdvect = pvect->GetVect ();

	pvect->SetSVect (0.0e0);

	int i, iblk, ni, ibs, j;

	int nloc = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pblks[iblk+1]-pblks[iblk];
		for (j=0;j<ni;j++) pdvect[ibs+j] = px[nloc+j];
		nloc += ni;
	};

// Prepare sends and receives

	int *psize;
	double **psendarr;

	psize = new int [nproc];
	if (!psize) MemoryFail (funcname);
	psendarr = new double * [nproc];
	if (!psendarr) MemoryFail (funcname);

	int jjblk, niloc, k;

	for (i=0;i<nproc;i++) {
		ni = 0;
		for (j=iasndblk[i];j<iasndblk[i+1];j++) {
			jjblk = jasndblk[j];
			ni += pblks[jjblk+1]-pblks[jjblk];
		};
		psendarr[i] = new double [ni];
		if (!psendarr[i]) MemoryFail (funcname);
		ni = 0;
		for (j=iasndblk[i];j<iasndblk[i+1];j++) {
			jjblk = jasndblk[j];
			ibs = ibsvect[jjblk];
			niloc = pblks[jjblk+1]-pblks[jjblk];
			for (k=0;k<niloc;k++) psendarr[i][ni+k] = pdvect[ibs+k];
			ni += niloc;
		};
		psize[i] = ni * sizeof(double);
	};

// Send/receive all necessary data

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
		ObjSizeSend[i] = psize[i];
		ObjSend[i] = (char *)psendarr[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (ptree->GetComm (),
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CMvmSchurR::ExchangeX: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] psendarr[i];

// Store the result

	int jproc;
	double *pdarr;

	for (i=0;i<NObjRecv;i++) {
		jproc = CpuIDRecv[i];
		pdarr = (double *)ObjRecv[i];
		ni = 0;
		for (j=iarcvblk[jproc];j<iarcvblk[jproc+1];j++) {
			jjblk = jarcvblk[j];
			ibs = ibsvect[jjblk];
			niloc = pblks[jjblk+1]-pblks[jjblk];
			for (k=0;k<niloc;k++) pdvect[ibs+k] = pdarr[ni+k];
			ni += niloc;
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

// Free work arrays

	delete [] psize;
	delete [] psendarr;

};

// Author: Kharchenko S.A.
// Description: Init necessary Ax by zeroes
// CMvmSchurR::InitAxByZeroes()
//========================================================================================
void CMvmSchurR::InitAxByZeroes () { // Init necessary Ax by zeroes

	pvect1->SetSVect (0.0e0);

};

// Author: Kharchenko S.A.
// Description: Multiply by local block rows/columns of A
// CMvmSchurR::PerformLocalMultiplications()
//========================================================================================
void CMvmSchurR::PerformLocalMultiplications (char _transp) { // Multiply by local block rows/columns of A

//	int myid = ptree->GetMyid ();

	int i, iblk;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		if (_transp == 'N' || _transp == 'n') {
			(*pmvmau) (pgmtrau, iblk, ibsvect, *pvect, ibsvect, *pvect1);
			(*pmvmal) (pgmtral, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else {
			(*pmvmau) (pgmtral, iblk, ibsvect, *pvect, ibsvect, *pvect1);
			(*pmvmal) (pgmtrau, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		};
	};

};

// Author: Kharchenko S.A.
// Description: Multiply by the current block row of Au
// CSMatrixR::MvmBlockRowAu2Index()
//========================================================================================
void CSMatrixR::AbsMvmBlockRowAu2Index (void *_obj, int _iblk, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax) {

	CSMatrixR *pmtrarralu = ((CGSMatrixR *)_obj)->GetMtrarr ();

	pmtrarralu[_iblk].MvmBlockRowAu2Index (_bsx, _x, _bsax, _ax);

};

// Author: Kharchenko S.A.
// Description: Multiply by the current block column of Al without diagonal supernode
// CSMatrixR::MvmBlockColumnAl2Index()
//========================================================================================
void CSMatrixR::AbsMvmBlockColumnAl2Index (void *_obj, int _iblk, // Multiply by the current block column of Al without diagonal supernodeL
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax) {

	CSMatrixR *pmtrarralu = ((CGSMatrixR *)_obj)->GetMtrarr ();

	pmtrarralu[_iblk].MvmBlockColumnAl2Index (false, _bsx, _x, _bsax, _ax);

};

// Author: Kharchenko S.A.
// Description: Multiply by the current block row of Au
// CParSchurMvmSlvRS::MvmBlockRowAu2Index()
//========================================================================================
void CParSchurMvmSlvRS::AbsMvmBlockRowAu (void *_obj, int _iblk, // Multiply by the current block row of Au
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax) {

	CParSchurMvmSlvRS *pobj = (CParSchurMvmSlvRS *)_obj;

	char optypeloc = pobj->GetOptype ();
	int *ibssuploc = pobj->GetIbssup ();

	CGSMatrixRS *pgmtrauloc;

	if (optypeloc == 'N') {
		pgmtrauloc = pobj->pgmtrau;
	} else {
		pgmtrauloc = pobj->pgmtral;
	};

	int *psprndsloc = pgmtrauloc->GetSprndr ();

	CSMatrixRS *pmtrarralu = ((CGSMatrixRS *)pgmtrauloc)->GetMtrarr ();

	pmtrarralu[_iblk].MvmBlockRowAu (psprndsloc, ibssuploc, _x, ibssuploc, _ax);

};

// Author: Kharchenko S.A.
// Description: Multiply by the current block column of Al without diagonal supernode
// CParSchurMvmSlvRS::MvmBlockColumnAl()
//========================================================================================
void CParSchurMvmSlvRS::AbsMvmBlockColumnAl (void *_obj, int _iblk, // Multiply by the current block column of Al without diagonal supernodeL
								int *_bsx, const CSVector &_x, 
								int *_bsax, CSVector &_ax) {

	CParSchurMvmSlvRS *pobj = (CParSchurMvmSlvRS *)_obj;

	char optypeloc = pobj->GetOptype ();
	int *ibssuploc = pobj->GetIbssup ();

	CGSMatrixRS *pgmtralloc;

	if (optypeloc == 'N') {
		pgmtralloc = pobj->pgmtral;
	} else {
		pgmtralloc = pobj->pgmtrau;
	};

	int *psprndsloc = pgmtralloc->GetSprndr ();

	CSMatrixRS *pmtrarralu = ((CGSMatrixRS *)pgmtralloc)->GetMtrarr ();

	pmtrarralu[_iblk].MvmBlockColumnAl (false, psprndsloc, ibssuploc, _x, ibssuploc, _ax);

};

// Author: Kharchenko S.A.
// Description: Exchange X data
// CMvmSchurR::ExchangeAX()
//========================================================================================
void CMvmSchurR::ExchangeAX (CSVector &_ax) { // Exchange AX data

	const char *funcname = "ExchangeAX";

//	int myid = ptree->GetMyid ();
	int nproc = ptree->GetNproc ();

// Prepare sends and receives

	double *pdvect = pvect1->GetVect ();

	int *psize;
	double **psendarr;

	psize = new int [nproc];
	if (!psize) MemoryFail (funcname);
	psendarr = new double * [nproc];
	if (!psendarr) MemoryFail (funcname);

	int i, iblk, ni, ibs, j;
	int jjblk, niloc, k;

	for (i=0;i<nproc;i++) {
		ni = 0;
		for (j=iarcvblk[i];j<iarcvblk[i+1];j++) {
			jjblk = jarcvblk[j];
			ni += pblks[jjblk+1]-pblks[jjblk];
		};
		psendarr[i] = new double [ni];
		if (!psendarr[i]) MemoryFail (funcname);
		ni = 0;
		for (j=iarcvblk[i];j<iarcvblk[i+1];j++) {
			jjblk = jarcvblk[j];
			ibs = ibsvect[jjblk];
			niloc = pblks[jjblk+1]-pblks[jjblk];
			for (k=0;k<niloc;k++) psendarr[i][ni+k] = pdvect[ibs+k];
			ni += niloc;
		};
		psize[i] = ni * sizeof(double);
	};

// Send/receive all necessary data

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
		ObjSizeSend[i] = psize[i];
		ObjSend[i] = (char *)psendarr[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (ptree->GetComm (),
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CMvmSchurR::ExchangeX: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] psendarr[i];

// Store the result

	int jproc;
	double *pdarr;

	for (i=0;i<NObjRecv;i++) {
		jproc = CpuIDRecv[i];
		pdarr = (double *)ObjRecv[i];
		ni = 0;
		for (j=iasndblk[jproc];j<iasndblk[jproc+1];j++) {
			jjblk = jasndblk[j];
			ibs = ibsvect[jjblk];
			niloc = pblks[jjblk+1]-pblks[jjblk];
			for (k=0;k<niloc;k++) pdvect[ibs+k] += pdarr[ni+k];
			ni += niloc;
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

// Free work arrays

	delete [] psize;
	delete [] psendarr;

// Store AX part

	pvect->SetSVect (0.0e0);

	int nloc = 0;

	double *pax = _ax.GetVect ();

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pblks[iblk+1]-pblks[iblk];
		for (j=0;j<ni;j++) pax[nloc+j] = pdvect[ibs+j];
		nloc += ni;
	};

};


// Author: Kharchenko S.A.
// Description: Init the set of Px blocks
// CMvmSchurR::InitPxBlocks()
//========================================================================================
void CMvmSchurR::InitPxBlocks (CSVector &_x) { // Init the set of Px blocks

//	const char *funcname = "InitPxBlocks";

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Init X part

	double *px = _x.GetVect ();
	double *pdvect = pvect->GetVect ();

	pvect->SetSVect (0.0e0);
	pvect1->SetSVect (0.0e0);

	int i, iblk, ni, ibs, j;

	int nloc = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pblks[iblk+1]-pblks[iblk];
		for (j=0;j<ni;j++) pdvect[ibs+j] = px[nloc+j];
		nloc += ni;
	};

};

// Author: Kharchenko S.A.
// Description: Solve L with the set of blocks
// CMvmSchurR::SolveLBlocks()
//========================================================================================
void CMvmSchurR::SolveLBlocks (char _lutype, int _nlistfct, int *_listfct) { // Solve L with the set of blocks

//	int myid = ptree->GetMyid ();

	int i, iblk;

//	*pfout << " Before slv L " << endl;
//	OutArr (*pfout, " List of blocks slv ",_nlistfct,_listfct);
//	*pfout << " pvect on entry " << *pvect << endl;
//	*pfout << " pvect1 on entry " << *pvect1 << endl;
	for (i=0;i<_nlistfct;i++) {
		iblk = _listfct[i];
		if (_lutype == 'L' || _lutype == 'l') {
			(*psolvel) (pgmtrl, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else if (_lutype == 'A' || _lutype == 'a') {
			(*psolvel) (pgmtral, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else if (_lutype == 'U' || _lutype == 'u') {
			(*psolvel) (pgmtru, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		};
	};
//	*pfout << " pvect1 on return " << *pvect1 << endl;

};

// Author: Kharchenko S.A.
// Description: Solve triangular system with block column of L
// CSMatrixR::SolveBlockColumnL2Index()
//========================================================================================
void CSMatrixR::AbsSolveBlockColumnL2Index (void *_obj, int _iblk, // Solve triangular system with block column of L
								int *_bsx, const CSVector &_x, 
								int *_bslx, CSVector &_lx) {

	CSMatrixR *pmtrarrlu = ((CGSMatrixR *)_obj)->GetMtrarr ();

	pmtrarrlu[_iblk].SolveBlockColumnL2Index (_bsx, _x, _bslx, _lx);

};

// Author: Kharchenko S.A.
// Description: Solve triangular system with block column of L
// CParSchurMvmSlvRS::AbsSolveBlockColumnL()
//========================================================================================
void CParSchurMvmSlvRS::AbsSolveBlockColumnL (void *_obj, int _iblk, // Solve triangular system with block column of L
								int *_bsx, const CSVector &_x, 
								int *_bslx, CSVector &_lx) {

	CParSchurMvmSlvRS *pobj = (CParSchurMvmSlvRS *)_obj;

	char optypeloc = pobj->GetOptype ();
	int *ibssuploc = pobj->GetIbssup ();

	CGSMatrixRS *pgmtrlloc;

	if (optypeloc == 'L') {
		pgmtrlloc = pobj->pgmtrl;
	} else {
		pgmtrlloc = pobj->pgmtru;
	};

	int *psprndsloc = pgmtrlloc->GetSprndr ();

	CSMatrixRS *pmtrarrlu = ((CGSMatrixRS *)pgmtrlloc)->GetMtrarr ();

	pmtrarrlu[_iblk].SolveBlockColumnL (psprndsloc, ibssuploc, _x, ibssuploc, _lx);

};

// Author: Kharchenko S.A.
// Description: Solve U with the set of blocks
// CMvmSchurR::SolveUBlocks()
//========================================================================================
void CMvmSchurR::SolveUBlocks (char _lutype, int _nlistfct, int *_listfct) { // Solve U with the set of blocks

//	int myid = ptree->GetMyid ();

	int i, iblk;

	for (i=0;i<_nlistfct;i++) {
		iblk = _listfct[i];
		if (_lutype == 'L' || _lutype == 'l') {
			(*psolveu) (pgmtrl, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		} else {
			(*psolveu) (pgmtru, iblk, ibsvect, *pvect, ibsvect, *pvect1);
		};
	};

};

// Author: Kharchenko S.A.
// Description: Solve triangular system with the block row of U
// CSMatrixR::SolveBlockRowU2Index()
//========================================================================================
void CSMatrixR::AbsSolveBlockRowU2Index (void *_obj, int _iblk, // Solve triangular system with the block row of U
								int *_bsx, const CSVector &_x, 
								int *_bsux, CSVector &_ux) {

	CSMatrixR *pmtrarrlu = ((CGSMatrixR *)_obj)->GetMtrarr ();

	pmtrarrlu[_iblk].SolveBlockRowU2Index (_bsx, _x, _bsux, _ux);

};

// Author: Kharchenko S.A.
// Description: Solve triangular system with the block row of U
// CParSchurMvmSlvRS::SolveBlockRowU()
//========================================================================================
void CParSchurMvmSlvRS::AbsSolveBlockRowU (void *_obj, int _iblk, // Solve triangular system with the block row of U
								int *_bsx, const CSVector &_x, 
								int *_bsux, CSVector &_ux) {

	CParSchurMvmSlvRS *pobj = (CParSchurMvmSlvRS *)_obj;

	char optypeloc = pobj->GetOptype ();
	int *ibssuploc = pobj->GetIbssup ();

	CGSMatrixRS *pgmtruloc;

	if (optypeloc == 'L') {
		pgmtruloc = pobj->pgmtrl;
	} else {
		pgmtruloc = pobj->pgmtru;
	};

	int *psprndsloc = pgmtruloc->GetSprndr ();

	CSMatrixRS *pmtrarrlu = ((CGSMatrixRS *)pgmtruloc)->GetMtrarr ();

	pmtrarrlu[_iblk].SolveBlockRowU (psprndsloc, ibssuploc, _x, ibssuploc, _ux);

};

// Author: Kharchenko S.A.
// Description: Send the set of Px blocks
// CMvmSchurR::SendBlocks()
//========================================================================================
void CMvmSchurR::SendBlocks (int _nlistcpu, int *_listcpu, int _nlistschur, int *_listschur) { // Send the set of Px blocks

	const char *funcname = "SendBlocks";

	if (_nlistschur == 0) return;

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Prepare the send buffer

	int i, iblk, niloc;

	int ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		niloc = pblks[iblk+1]-pblks[iblk];
		ni += niloc;
	};

	double *darr;

	darr = new double [ni];
	if (!darr) MemoryFail (funcname);

	double *pdvect = pvect1->GetVect ();

	int ibs, j;

	ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		ibs = ibsvect[iblk];
		niloc = pblks[iblk+1]-pblks[iblk];
		for (j=0;j<niloc;j++) darr[ni+j] = pdvect[ibs+j];
		ni += niloc;
	};
//	OutArr (*pfout," Send to CPUs ",_nlistcpu,_listcpu);
//	OutArr (*pfout," Send list ",_nlistschur,_listschur);
//	OutArr (*pfout," IbsVect ",nblks,ibsvect);
//	OutArr (*pfout," Blks ",nblks+1,pblks);
//	OutArr (*pfout," Send arr ",ni,darr);

// Send the set of blocks to the list of cpu's

	int jproc;

	for (i=0;i<_nlistcpu;i++) {
		jproc = _listcpu[i];

		if (nsends >= nsendsmax) throw " CMvmSchurR::SendBlocks: insufficient number of sends ";

		CMPIExchange::ISend (ptree->GetComm (), jproc, 1,
									ni*sizeof(double), (char *)darr, sendreqvarr[nsends]);

		imasksend[nsends] = 1;
		if (i==0) {
			psenddata[nsends] = (char *)darr;
		} else {
			psenddata[nsends] = 0;
		};
		nsends++;
	};

};

// Author: Kharchenko S.A.
// Description: Receive the set of Px blocks
// CMvmSchurR::ReceiveBlocks()
//========================================================================================
void CMvmSchurR::ReceiveBlocks (bool _add, int _jproc, int _nlistschur, int *_listschur) { // Receive the set of Px blocks

	const char *funcname = "ReceiveBlocks";

	if (_nlistschur == 0) return;

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Receive the set of blocks

	int i, iblk, niloc;

	int ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		niloc = pblks[iblk+1]-pblks[iblk];
		ni += niloc;
	};

	double *darr;

	darr = new double [ni];
	if (!darr) MemoryFail (funcname);

	CMPIStatus status;

	CMPIExchange::Recv (ptree->GetComm (), _jproc, 1,
								ni*sizeof(double), (char *)darr, status);

//	OutArr (*pfout," Recv from ",1,&_jproc);
//	OutArr (*pfout," Recv list ",_nlistschur,_listschur);
//	OutArr (*pfout," Recv arr ",ni,darr);

// Store the data

	double *pdvect = pvect1->GetVect ();

	int ibs, j;

	ni = 0;

	for (i=0;i<_nlistschur;i++) {
		iblk = _listschur[i];
		ibs = ibsvect[iblk];
		niloc = pblks[iblk+1]-pblks[iblk];
		if (_add) {
			for (j=0;j<niloc;j++) pdvect[ibs+j] += darr[ni+j];
		} else {
			for (j=0;j<niloc;j++) pdvect[ibs+j] = darr[ni+j];
		};
		ni += niloc;
	};

// Free work data

	delete [] darr;

};

// Author: Kharchenko S.A.
// CMvmSchurR: Wait for completion of sends
//========================================================================================
void CMvmSchurR::WaitSends () { // Wait for completion of sends

	const char *funcname = "WaitSends";

	int nloc = nsends;

	CMPIRequest *reqvarrloc;
	CMPIStatus *statarrloc;

	reqvarrloc = new CMPIRequest [nloc];
	if (!reqvarrloc) MemoryFail (funcname);
	statarrloc = new CMPIStatus [nloc];
	if (!statarrloc) MemoryFail (funcname);

	nloc = 0;

	int i;

	for (i=0;i<nsends;i++) {
		if (imasksend[i] == 1) {
			reqvarrloc[nloc] = sendreqvarr[i];
			nloc++;
		};
	};

	CMPIExchange::WaitAll (nloc, reqvarrloc, statarrloc);

	for (i=0;i<nsends;i++) {
		if (psenddata[i] != 0) {
			delete [] psenddata[i];
		};
	};

	delete [] reqvarrloc;
	delete [] statarrloc;

	nsends = 0;

};

// Author: Kharchenko S.A.
// Description: Store the result
// CMvmSchurR::StorePx()
//========================================================================================
void CMvmSchurR::StorePx (CSVector &_px) { // Store the result

//	const char *funcname = "StorePx";

//	int myid = ptree->GetMyid ();
//	int nproc = ptree->GetNproc ();

// Init X part

	double *px = _px.GetVect ();
	double *pdvect = pvect1->GetVect ();

	int i, iblk, ni, ibs, j;

	int nloc = 0;

	for (i=0;i<nlistblk;i++) {
		iblk = listblk[i];
		ibs = ibsvect[iblk];
		ni = pblks[iblk+1]-pblks[iblk];
		for (j=0;j<ni;j++) px[nloc+j] = pdvect[ibs+j];
		nloc += ni;
	};

};

// Author: Kharchenko S.A.
// Description: Destroy solve working data
// CMvmSchurR::DestroyWorkData()
//========================================================================================
void CMvmSchurR::DestroyWorkData () { // Destroy solve working data

	const char *funcname = "DestroyWorkData";

	delete [] imaskflt;
	delete [] ibsvect;
	delete [] listblk;
	delete [] listblkext;
	delete [] iasndblk;
	delete [] jasndblk;
	delete [] iarcvblk;
	delete [] jarcvblk;

	imaskflt = new int [0];
	if (!imaskflt) MemoryFail (funcname);
	ibsvect = new int [0];
	if (!ibsvect) MemoryFail (funcname);
	listblk = new int [0];
	if (!listblk) MemoryFail (funcname);
	listblkext = new int [0];
	if (!listblkext) MemoryFail (funcname);
	iasndblk = new int [0];
	if (!iasndblk) MemoryFail (funcname);
	jasndblk = new int [0];
	if (!jasndblk) MemoryFail (funcname);
	iarcvblk = new int [0];
	if (!iarcvblk) MemoryFail (funcname);
	jarcvblk = new int [0];
	if (!jarcvblk) MemoryFail (funcname);

	nlistblk = 0;
	nlistblkext = 0;

	delete [] imasksend;
	delete [] psenddata;
	delete [] sendreqvarr;

	nsends = 0;
	nsendsmax = 0;

	imasksend = new int [nsendsmax];
	if (!imasksend) MemoryFail (funcname);
	psenddata = new char * [nsendsmax];
	if (!psenddata) MemoryFail (funcname);
	sendreqvarr = new CMPIRequest [nsendsmax];
	if (!sendreqvarr) MemoryFail (funcname);

	delete pvect;
	delete pvect1;

	pvect = new CSVector;
	if (!pvect) MemoryFail (funcname);
	pvect1 = new CSVector;
	if (!pvect1) MemoryFail (funcname);

};

// Author: Kharchenko S.A.
// Description: Destroy solve working data
// CParSchurMvmSlvRS::CreateIbssup()
//========================================================================================
void CParSchurMvmSlvRS::CreateIbssup () { // Create ibssup data for local operations

	const char *funcname = "CreateIbssup";

	if (ibssup != 0) {
		delete [] ibssup;
	};

	int *pibsvect = pmvm->GetIbsvect ();
	int nblksloc = ptree->GetNblkstree ();
	int *pblksloc = ptree->GetBlkstree ();
	int nsuploc = pgmtral->GetNsupr ();
	int *sprndsloc = pgmtral->GetSprndr ();

	ibssup = new int [nsuploc];
	if (!ibssup) MemoryFail (funcname);

	int i, iblk, ibs, j;

	for (i=0;i<nsuploc;i++) ibssup[i] = -1;

	for (iblk=0;iblk<nblksloc;iblk++) {
		if (pibsvect[iblk]>= 0) {
			ibs = pibsvect[iblk];
			for (j=pblksloc[iblk];j<pblksloc[iblk+1];j++) {
				ibssup[j] = ibs;
				ibs += sprndsloc[j+1]-sprndsloc[j];
			};
		};
	};

};
