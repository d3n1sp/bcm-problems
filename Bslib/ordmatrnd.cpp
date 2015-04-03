#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "globals.h"
#include "smatrix.h"
#include "gsmatrix.h"
#include "tree.h"
#include "slvparam.h"
#include "ExchangeMPI.h"

using namespace std;

#ifdef __FV_2004__
namespace flowvision {
namespace EqnSolver {
#endif
/*
// Author: Kharchenko S.A.
// CSMatrix: Compute Nested Dissection ordering of the matrix in parallel via ParMeTiS with Schur partitioning
//========================================================================================
void CSMatrix::OrderNDParMetisSchur (CMPIComm &_comm, int _nparts, // Compute Nested Dissection ordering of the matrix in parallel via ParMeTiS with Schur partitioning
													CTree &_tree, 
													int *&_order, int &_nblks, int *&_blks, int *&_blk2cpu) const {

	const char *funcname = "OrderNDParMetisSchur";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Perform initial ParMetis partitioning

	CTree treeloc;
	int *orderloc;
	int nblksloc;
	int *blksloc, *blk2cpuglb;

	OrderNDParMetis (_comm, _nparts,
							treeloc, 
							orderloc, nblksloc, blksloc, blk2cpuglb);

// Create a smaller binary tree

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nproc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nproc+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nproc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nproc;i++) procweight[i] = 1.0e0;

	CTree treesm (nproc, 2, procweight, memory2cpu);

	treesm.SetMyid (myid);

	delete [] memory2cpu;
	delete [] procweight;

	treesm.SortNodes ();

// Include extended binary tree into a smaller one

	CSMatrix *pmtr = treeloc.GetAblktree ();

	treeloc.IncludeBinaryTreeSchur (treesm, nblksloc, blksloc, *pmtr);

// Determine new local processors partitioning

	int *blk2cpuloc;

	blk2cpuloc = new int [nblksloc+1];
	if (!blk2cpuloc) MemoryFail (funcname);

	treesm.Block2Cpu (blk2cpuloc);

// Create the set of intervals

	int *ibegarrloc;
	int *iendarrloc;

	ibegarrloc = new int [nblksloc];
	if (!ibegarrloc) MemoryFail (funcname);
	iendarrloc = new int [nblksloc];
	if (!iendarrloc) MemoryFail (funcname);

	int npartsloc = 0;
	int nloc = 0;

	int i, j;

	for (i=0;i<nblksloc;i++) {
		if (blk2cpuloc[i] == myid && blksloc[i+1] > blksloc[i]) {
			ibegarrloc[npartsloc] = blksloc[i];
			iendarrloc[npartsloc] = blksloc[i+1]-1;
			npartsloc++;
			nloc += blksloc[i+1]-blksloc[i];
		};
	};

	int *listtree;

	listtree = new int [nloc];
	if (!listtree) MemoryFail (funcname);

	for (i=0;i<nblksloc;i++) {
		if (blk2cpuloc[i] == myid && blksloc[i+1] > blksloc[i]) {
			for (j=blksloc[i];j<blksloc[i+1];j++) {
				listtree[nloc] = j;
				nloc++;
			};
		};
	};

// Create direct and inverse order for the set of intervals

	int nlistini = nlist;
	int *plistini = list;

	int nlisttree;
	int *orderpart;

//	CSMatrix::DirOrderViaIntervals (_comm,
//												nlistini, plistini, orderloc,
//												npartsloc, ibegarrloc, iendarrloc, 
//												nlisttree, orderpart);
	CSMatrix::InvOrderViaIntervals (_comm,
												nlistini, orderloc, plistini,
												npartsloc, ibegarrloc, iendarrloc, 
												nlisttree, orderpart);

	int *invordertree;

	CSMatrix::InvOrderViaIntervals (_comm,
												nlistini, plistini, orderloc,
												npartsloc, ibegarrloc, iendarrloc, 
												nlisttree, invordertree);

// Perform reordering in parallel

	CSMatrix asord = *this;

	asord.ChangeIndicesViaIntervals (_comm,
												orderloc,
												nlisttree, listtree, invordertree,
												npartsloc, ibegarrloc, iendarrloc, orderpart);

// Perform matrix structure repartitioning

	CSMatrix asnew;

	asnew = asord.GetSubmatrixViaIntervals (_comm, 
															npartsloc, ibegarrloc, iendarrloc);

	asord = mtrdummy;

// Main cycle over the levels of the tree

	CSMatrix aschur;
	CSMatrix aschur_blk;

	CSMatrix *mtrarrloc;
	CSMatrix *mtrblkarrloc;

	mtrarrloc = new CSMatrix [nproc];
	if (!mtrarrloc) MemoryFail (funcname);
	mtrblkarrloc = new CSMatrix [nproc];
	if (!mtrblkarrloc) MemoryFail (funcname);

	int isize;
	char *obj;

	int nnodesext = treeloc.GetNnodes ();
	CNode *pnodes = treeloc.GetNode ();

	int *imaskcpu;
	int *listcpu;

	imaskcpu = new int [nproc];
	if (!imaskcpu) MemoryFail (funcname);
	listcpu = new int [nproc];
	if (!listcpu) MemoryFail (funcname);

	for (i=0;i<nproc;i++) imaskcpu[i] = -1;

	int icycle = -1;

	int *piaglb = pmtr->GetIa ();
	int *pjaglb = pmtr->GetJa ();

	for (inode=0;inode<nnodesext;inode++) {

// Get node information

		iblkbeg = pnodes[inode].GetIndbeg ();
		iblkend = pnodes[inode].GetIndend ();
		nchildsloc = pnodes[inode].GetNchilds ();

		if (nchildsloc == 1) break;

		ibegbrd = blksloc[iblkbeg];
		iendbrd = blksloc[iblkend+1]-1;
		ni = iendbrd-ibegbrd+1;

// Get dependency list

		icheck = -1;

		rootloc = blk2cpuloc[iblkbeg];

		if (rootloc == myid) icheck = 1;

		icycle++;

		nlistcpu = 0;

		imaskcpu[rootloc] = icycle;

		for (iblk=iblkbeg;iblk<=iblkend;iblk++) {
			for (j=piaglb[iblk];j<piaglb[iblk+1];j++) {
				jblk = pjaglb[j];
				iproc = blk2cpuloc[jblk];
				if (jblk < iblkbeg && imaskcpu[iproc] != icycle) {
					listcpu[nlistcpu] = jblk;
					nlistcpu++;
					imaskcpu[iproc] = icycle;
				};
			};
		};

		qsort (listcpu, nlistcpu, sizeof(int), compint);

		for (i=0;i<nlistcpu;i++) {
			if (listcpu[i] == myid) icheck = 1;
		};

		if (icheck == -1) break;

// Compute local component of the Schur complement

		asnew.ComputeSchur (nblksloc, blksloc,
									ibegbrd, iendbrd, ?ncycle,
									aschur, aschur_blk);

// Shift submatrices data

		ishift = -ibegbrd;

		aschur.ShiftSparsity (true, ishift);
		aschur_blk.ShiftSparsity (false, ishift);

// Exchange by the local components if necessary

		if (myid != rootloc) {
			aschur.PackMatrix (isize,obj);
			aschur = mtrdummy;
			CMPIExchange::Send (_comm, rootloc, 1, sizeof(int), (char *)&isize);
			CMPIExchange::Send (_comm, rootloc, 2, isize, obj);
			delete [] obj;
			aschur_blk.PackMatrix (isize,obj);
			aschur_blk = mtrdummy;
			CMPIExchange::Send (_comm, rootloc, 3, sizeof(int), (char *)&isize);
			CMPIExchange::Send (_comm, rootloc, 4, isize, obj);
			delete [] obj;
		} else {
			for (i=0;i<nlistcpu;i++) {
				iproc = listcpu[i];
				CMPIExchange::Recv (_comm, iproc, 1, sizeof(int), (char *)&isize, stat);
				Recv (iproc,isize);
				obj = new char [isize];
				if (!obj) MemoryFail (funcname);
				CMPIExchange::Recv (_comm, rootloc, 2, isize, obj, stat);
				mtrarrloc[i].UnPackMatrix (isize,obj);
				delete [] obj;
				CMPIExchange::Recv (_comm, iproc, 3, sizeof(int), (char *)&isize, stat);
				obj = new char [isize];
				if (!obj) MemoryFail (funcname);
				CMPIExchange::Recv (_comm, rootloc, 4, isize, obj, stat);
				mtrarrloc[i].UnPackMatrix (isize,obj);
				delete [] obj;
			};
		};

// Add local components if necessary

		if (myid == rootloc) {

			adiag = CSMatrix::DiagonalMtr (ni);

			aschurtot = adiag+aschur;

			for (i=0;i<nlistcpu;i++) {
				atmp = aschurtot + mtrarrloc[i];
				aschurtot = atmp;
				atmp = mtrdummy;
				mtrarrloc[i] = mtrdummy;
			};

			aschurblktot = aschur_blk;

			for (i=0;i<nlistcpu;i++) {
				atmp = aschurblktot + mtrblkarrloc[i];
				aschurblktot = atmp;
				atmp = mtrdummy;
				mtrblkarrloc[i] = mtrdummy;
			};

		};

// Compute local ordering and partitioning for the bordering

		if (myid == rootloc) {
			aschurtot.PartNDSchurBnd (aschurblktot, rootloc, blk2cpuglb,
												_param,
												order_schur, nblks_schur, blks_schur, blk2cpu_schur);
		};

// Register partitioning results in a tree

		????;

	};

// Exchange partitioning data

	????;

// Exchange ordering 

	????;

// Reorder the matrix

	????;

// Create new blocks partitioning

	????;

// Free work arrays

	delete [] orderloc;
	delete [] blksloc;
	delete [] blk2cpuloc;
	delete [] ibegarrloc;
	delete [] iendarrloc;
	delete [] listtree;
	delete [] orderpart;
	delete [] invordertree;

};
*/
// Author: Kharchenko S.A.
// CSMatrix: Compute Nested Dissection ordering of the matrix in parallel via ParMeTiS
//========================================================================================
void CSMatrix::OrderNDParMetis (CMPIComm &_comm, int _nparts, // Compute Nested Dissection ordering of the matrix in parallel via ParMeTiS
											CTree &_tree, int *&_order, int &_nblks, int *&_blks, int *&_blk2cpu) const {

	const char *funcname = "OrderNDParMetis";

// Check that the number of cpu's is a power of 2

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

	int np = 1;

	while (np*2 <= nproc) np *= 2;

	if (np != nproc) throw " CSMatrix::OrderNDParMetis: the number of processors actually used is not a power of 2 ";

	if (_nparts < nproc) throw " CSMatrix::OrderNDParMetis: too small number of parts requested ";

	np = 1;

	while (np*2 <= _nparts) np *= 2;

	if (np != _nparts) throw " CSMatrix::OrderNDParMetis: the number of parts requested is not a power of 2 ";

// Collect sizes statistics for the whole set of matrices

	int *blkscpu;

	blkscpu = new int [nproc+1];
	if (!blkscpu) MemoryFail (funcname);

	int i;

	for (i=0;i<nproc+1;i++) blkscpu[i] = 0;

	blkscpu[myid+1] = nlist;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blkscpu, blkscpu);

	for (i=0;i<nproc;i++) blkscpu[i+1] = blkscpu[i]+blkscpu[i+1];

// Check that matrix data are also correct

	int indbeg = list[0];
	int indend = list[nlist-1];

	if (indbeg != blkscpu[myid] || indend != blkscpu[myid+1]-1) throw " CSMatrix::OrderNDParMetis: wrong part of the matrix is considered";

	int ntot = blkscpu[nproc];

// Compute symmetrized matrix in parallel

	CSMatrix as;

//	as = SymmMtr ();
	as = *this;

// Call parmetis

	int nblksloc;
	int *blksloc;

	as.OrderNodeNDParMetis (_comm, _order, nblksloc, blksloc);

	CSMatrix mtrdummy;

// Create binary tree

	double *memory2cpu;
	double *procweight;

	memory2cpu = new double [nproc+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nproc+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nproc;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nproc;i++) procweight[i] = 1.0e0;

	CTree tree (nproc, 2, procweight, memory2cpu);

	tree.SetMyid (myid);

	delete [] memory2cpu;
	delete [] procweight;

	tree.SortNodes ();

	int nnodesloc = tree.GetNnodes ();

	CNode *pnodes = tree.GetNode ();

// Recompute initial blocks array according to the blks arrays

	int nlevloc = tree.GetNlev ();

	int *iplev;
	int *iptrlev;

	iplev = new int [nlevloc+2];
	if (!iplev) MemoryFail (funcname);
	iptrlev = new int [nlevloc+2];
	if (!iptrlev) MemoryFail (funcname);

	int ip = 1;

	iplev[0] = 0;

	for (i=0;i<nlevloc;i++) {
		iplev[i+1] = iplev[i] + ip;
		ip *= 2;
	};

	for (i=0;i<=nlevloc;i++) iptrlev[i]= iplev[i];

	for (i=0;i<nnodesloc;i++) {
		pnodes[i].indbeg = -1;
		pnodes[i].indend = -1;
		pnodes[i].indbegtot = -1;
		pnodes[i].indendtot = -1;
	};

	int inode, ilev, ichild;

	if (nnodesloc == 1) {
		pnodes[0].indbeg = 0;
		pnodes[0].indend = ntot-1;
	} else {
		for (i=0;i<nnodesloc;i++) {
			inode = i;
			if (pnodes[inode].nchilds > 1) {
				ilev = pnodes[inode].nodelv;
				ip = iptrlev[ilev];
				pnodes[inode].indbeg = 0;
				pnodes[inode].indend = blksloc[ip*4+3]-blksloc[ip*4+2]-1;
				if (ilev == nlevloc-2) {
					ichild = pnodes[inode].childs[0];
					pnodes[ichild].indbeg = 0;
					pnodes[ichild].indend = blksloc[ip*4+1]-blksloc[ip*4]-1;
					ichild = pnodes[inode].childs[1];
					pnodes[ichild].indbeg = 0;
					pnodes[ichild].indend = blksloc[ip*4+2]-blksloc[ip*4+1]-1;
				};
				iptrlev[ilev]++;
			};
		};
	};

	int nitot = 0;

	for (inode=0;inode<nnodesloc;inode++) {
		pnodes[inode].indbeg = nitot;
		nitot += (pnodes[inode].indend+1);
		pnodes[inode].indend = nitot-1;
		pnodes[inode].indbegtot = pnodes[inode].indbeg;
		pnodes[inode].indendtot = pnodes[inode].indend;
	};

// Determine indices which describe position of the submatrices for each cpu

	int *ibegcpuarr;
	int *iendcpuarr;

	ibegcpuarr = new int [nproc];
	if (!ibegcpuarr) MemoryFail (funcname);
	iendcpuarr = new int [nproc];
	if (!iendcpuarr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) {
		inode = tree.cpuidend[i];
		ibegcpuarr[i] = pnodes[inode].indbeg;
		iendcpuarr[i] = pnodes[inode].indend;
	};

// Create lists of row indices for each cpu (including own) by binary search

	CIntInt *iarr;

	iarr = new CIntInt [nproc];
	if (!iarr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) {
		iarr[i].intvalue = ibegcpuarr[i];
		iarr[i].int2value = i;
	};

	qsort (iarr,nproc,sizeof(CIntInt),CIntInt::CompareIntInt);

	for (i=0;i<nproc;i++) {
		if (iarr[i].int2value != i) throw " CSMatrix::OrderNDParMetis: wrong order array ";
	};

	int *list2cpu;

	list2cpu = new int [nlist];
	if (!list2cpu) MemoryFail (funcname);

	for (i=0;i<nlist;i++) list2cpu[i] = -1;

	int icpubeg, icpuend, icpu, irow;

	for (i=0;i<nlist;i++) {
		icpubeg = 0;
		icpuend = nproc-1;
		irow = _order[i];
		while (icpubeg != icpuend) {
			icpu = (icpuend+icpubeg) / 2;
			if (icpu < icpubeg) icpu = icpubeg;
			if (icpu > icpuend) icpu = icpuend;
			if (irow >= ibegcpuarr[icpu] && irow <= iendcpuarr[icpu]) {
				icpubeg = icpu;
				icpuend = icpu;
			} else if (irow < ibegcpuarr[icpu]) {
				icpuend = icpu-1;
				if (icpuend < icpubeg) icpuend = icpubeg;
			} else if (irow > iendcpuarr[icpu]) {
				icpubeg = icpu+1;
				if (icpubeg > icpuend) icpubeg = icpuend;
			};
		};
		icpu = icpubeg;
		if (irow >= ibegcpuarr[icpu] && irow <= iendcpuarr[icpu]) {
			list2cpu[i] = icpu;
		};
	};

// Collect submatrices corresponding to the cpu's

	int *nlistarr;
	int *nzjaarr;

	nlistarr = new int [nproc];
	if (!nlistarr) MemoryFail (funcname);
	nzjaarr = new int [nproc];
	if (!nzjaarr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) nlistarr[i] = 0;
	for (i=0;i<nproc;i++) nzjaarr[i] = 0;

	int iproc;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlistarr[iproc]++;
			nzjaarr[iproc] += (as.ia[i+1]-as.ia[i]);
		};
	};

	CSMatrix *mtrarrloc;

	mtrarrloc = new CSMatrix [nproc];
	if (!mtrarrloc) MemoryFail (funcname);

	int nlistloc, nzjaloc;

	for (i=0;i<nproc;i++) {
		nlistloc = nlistarr[i];
		nzjaloc = nzjaarr[i];
		CSMatrix mtrloc (nlistloc,nzjaloc);
		mtrarrloc[i] = mtrloc;
	};

	int *iptrcpu;

	iptrcpu = new int [nproc];
	if (!iptrcpu) MemoryFail (funcname);

	for (i=0;i<nproc;i++) mtrarrloc[i].ia[0] = 0;
	for (i=0;i<nproc;i++) iptrcpu[i] = 0;

	int k, nz, j, jj;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			k = iptrcpu[iproc];
			mtrarrloc[iproc].list[k] = as.list[i];
			nz = mtrarrloc[iproc].ia[k];
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				mtrarrloc[iproc].ja[nz] = as.ja[j];
				nz++;
			};
			mtrarrloc[iproc].ia[k+1] = nz;
			iptrcpu[iproc]++;
		};
	};

	as = mtrdummy;

	delete [] list2cpu;

// Exchange the data

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
		mtrarrloc[i].PackMatrix (ObjSizeSend[i],ObjSend[i]);
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
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

	if (NObjRecv != nproc) throw " CSMatrix::OrderNDParMetis: wrong number of the received objects ";

	for(i = 0; i < NObjRecv; i++) mtrarrloc[i].UnPackMatrix (ObjSizeRecv[i],ObjRecv[i]);

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Collect submatrices of the current cpu

	int *listblk;

	listblk = new int [nproc];
	if (!listblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) listblk[i] = i;

	CSMatrix mtrnew;

	CGSMatrix::CombineMatricesSort (ntot, nproc, listblk,
												mtrarrloc,
												mtrnew);

	for (i=0;i<nproc;i++) mtrarrloc[i] = mtrdummy;

// For each column number determine the cpu index

	int nlisttot = mtrnew.nlist;
	int nzjatot = mtrnew.nzja;

	int *ja2cpu;

	ja2cpu = new int [nzjatot];
	if (!ja2cpu) MemoryFail (funcname);

	for (j=0;j<nzjatot;j++) {
		jj = mtrnew.ja[j];
		icpubeg = 0;
		icpuend = nproc-1;
		while (icpubeg != icpuend) {
			icpu = (icpuend+icpubeg) / 2;
			if (icpu < icpubeg) icpu = icpubeg;
			if (icpu > icpuend) icpu = icpuend;
			if (jj >= blkscpu[icpu] && jj < blkscpu[icpu+1]) {
				icpubeg = icpu;
				icpuend = icpu;
			} else if (jj < blkscpu[icpu]) {
				icpuend = icpu-1;
				if (icpuend < icpubeg) icpuend = icpubeg;
			} else if (jj >= blkscpu[icpu+1]) {
				icpubeg = icpu+1;
				if (icpubeg > icpuend) icpubeg = icpuend;
			};
		};
		ja2cpu[j] = icpubeg;
	};

// Create cpu substructures

	int *iacpu;
	int *iptr;
	int *jaindcpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);
	jaindcpu = new int [nzjatot];
	if (!jaindcpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	for (j=0;j<nzjatot;j++) {
		icpu = ja2cpu[j];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	for (j=0;j<nzjatot;j++) {
		icpu = ja2cpu[j];
		k = iptr[icpu];
		jaindcpu[k] = j;
		iptr[icpu]++;
	};

	delete [] ja2cpu;

// Create and init global imask

	int nlistmx = 0;

	int ni;

	for (i=0;i<nproc;i++) {
		ni = blkscpu[i+1]-blkscpu[i];
		if (ni > nlistmx) nlistmx = ni;
	};

	int *imask;
	int *listloc;

	imask = new int [nlistmx];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nlistmx];
	if (!listloc) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nlistmx;i++) {
		imask[i] = icycle;
	};

// Create the lists of indices for each cpu

	int *nlistcpu;
	int **plistcpu;

	nlistcpu = new int [nproc];
	if (!nlistcpu) MemoryFail (funcname);
	plistcpu = new int * [nproc*2];
	if (!plistcpu) MemoryFail (funcname);

	int ind, ishft, jjloc;

	for (i=0;i<nproc;i++) {
		icycle++;
		ishft = blkscpu[i];
		nz = 0;
		for (j=iacpu[i];j<iacpu[i+1];j++) {
			ind = jaindcpu[j];
			jj = mtrnew.ja[ind];
			jjloc = jj-ishft;
			if (imask[jjloc] != icycle) {
				listloc[nz] = jj;
				nz++;
				imask[jjloc] = icycle;
			};
		};
		nlistcpu[i] = nz;
		plistcpu[i] = new int [nz];
		if (!plistcpu[i]) MemoryFail (funcname);
		plistcpu[nproc+i] = new int [nz];
		if (!plistcpu[nproc+i]) MemoryFail (funcname);
		for (j=0;j<nz;j++) plistcpu[i][j] = listloc[j];
	};

// Exchange lists

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = nlistcpu[i]*sizeof(int);
		ObjSend[i] = (char *) plistcpu[i];
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	ishft = blkscpu[myid];

	int jjnew;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / sizeof(int);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[j];
			jjloc = jj-ishft;
			jjnew = _order[jjloc];
			piarr[j] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store new indices numbers

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / sizeof(int);
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			plistcpu[nproc+icpu][j] = piarr[j];
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Mark diagonal elements

	int *iptrdia;

	iptrdia = new int [nlisttot];
	if (!iptrdia) MemoryFail (funcname);

	for (i=0;i<nlisttot;i++) {
		int icheck = -1;
		irow = mtrnew.list[i];
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			jj = mtrnew.ja[j];
			if (jj == irow) {
				iptrdia[i] = j;
				icheck = 1;
			};
		};
		if (icheck == -1) throw " CSMatrix::OrderNDParMetis: diagonal element not found";
	};

// Perform renumbering of the columns

	for (icpu=0;icpu<nproc;icpu++) {

		ishft = blkscpu[icpu];

// Init list data

		ni = nlistcpu[icpu];

		for (i=0;i<ni;i++) {
			jj = plistcpu[icpu][i];
			jjloc = jj-ishft;
			jjnew = plistcpu[icpu+nproc][i];
			listloc[jjloc] = jjnew;
		};

// Renumber the columns

		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ind = jaindcpu[j];
			jj = mtrnew.ja[ind];
			jjloc = jj-ishft;
			jjnew = listloc[jjloc];
			mtrnew.ja[ind] = jjnew;
		};

	};

	delete [] jaindcpu;

	for (i=0;i<2*nproc;i++) delete [] plistcpu[i];

// Renumber the rows

	for (i=0;i<nlisttot;i++) {
		j = iptrdia[i];
		mtrnew.list[i] = mtrnew.ja[j];
	};

	delete [] iptrdia;

// Sort matrix data

	mtrnew.SortListAndColumns ();

// Filter off block diagonal data and shift the column indices

	int *ianew;

	ianew = new int [mtrnew.nlist+1];
	if (!ianew) MemoryFail (funcname);

	int ibegdia = ibegcpuarr[myid];
	int ienddia = iendcpuarr[myid];

	nz = 0;
	ianew[0] = 0;

	for (i=0;i<mtrnew.nlist;i++) {
		mtrnew.list[i] -= ibegdia;
		for (j=mtrnew.ia[i];j<mtrnew.ia[i+1];j++) {
			jj = mtrnew.ja[j];
			if (jj>=ibegdia && jj<=ienddia) {
				jj -= ibegdia;
				mtrnew.ja[nz] = jj;
				nz++;
			};
		};
		ianew[i+1] = nz;
	};

	for (i=0;i<=mtrnew.nlist;i++) mtrnew.ia[i] = ianew[i];

	mtrnew.nzja = nz;

	delete [] ianew;

// Create local and extended tree of processors

	int nchilds = 2;

	int nprocext = _nparts / nproc;

	np = 1;

	while (np*2 <= nprocext) np *= 2;

	if (np != nprocext) throw " CSMatrix::OrderNDParMetis: the number of additional parts is not a power of 2 ";

	memory2cpu = new double [nprocext+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [nprocext+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<nprocext;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<nprocext;i++) procweight[i] = 1.0e0;

	CTree treelocext (nprocext, nchilds, procweight, memory2cpu);

	treelocext.SetMyid (myid);

	delete [] memory2cpu;
	delete [] procweight;

// Sort the tree

	treelocext.SortNodes ();

// Compute partitioning of the submatrix

//	char strbuff[256];
//	sprintf (strbuff,"%s%i%s","OrdMtr_",myid,".ps");
//	mtrnew.A2Ps (100,strbuff,0,&i);

	mtrnew.PartBinaryTreeND (treelocext);

// Create packed object

	int isize;
	char *ptree;

	treelocext.PackTree (isize,ptree);

// Collect all the tree's produced by different cpu's

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = isize;
		ObjSend[i] = ptree;
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	delete [] ptree;

	CTree *treearr;

	treearr = new CTree [nproc];
	if (!treearr) MemoryFail (funcname);

	for (i=0;i<nproc;i++) {
		iproc = CpuIDRecv[i];
		treearr[iproc].UnPackTree (ObjSizeRecv[i],ObjRecv[i]);
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Create new tree

	memory2cpu = new double [_nparts+1];
	if (!memory2cpu) MemoryFail (funcname);
	procweight = new double [_nparts+1];
	if (!procweight) MemoryFail (funcname);

	for (i=0;i<_nparts;i++) memory2cpu[i] = 1.0e0;
	for (i=0;i<_nparts;i++) procweight[i] = 1.0e0;

	nchilds = 2;

	CTree treeloc (_nparts, nchilds, procweight, memory2cpu);

	delete [] memory2cpu;
	delete [] procweight;

	_tree = treeloc;

	_tree.SortNodes ();

	for (i=0;i<_tree.nnodes;i++) {
		_tree.nodes[i].indbeg = -1;
		_tree.nodes[i].indend = -1;
		_tree.nodes[i].indbegtot = -1;
		_tree.nodes[i].indendtot = -1;
	};

	int *ndold2new;

	ndold2new = new int [tree.nnodes];
	if (!ndold2new) MemoryFail (funcname);

	for (i=0;i<tree.nnodes;i++) ndold2new[i] = -1;

	ndold2new[tree.nnodes-1] = _tree.nnodes-1;

	int ichild0new, ichild1new, inodeloc, inodeold, inodenew, inodeshft;
	int ichild0, ichild1;

	for (inodeold=tree.nnodes-1;inodeold>=0;inodeold--) {
		inodenew = ndold2new[inodeold];
		if (tree.nodes[inodeold].nchilds > 1) {
			ichild0 = tree.nodes[inodeold].childs[0];
			ichild1 = tree.nodes[inodeold].childs[1];
			ichild0new = _tree.nodes[inodenew].childs[0];
			ichild1new = _tree.nodes[inodenew].childs[1];
			ndold2new[ichild0] = ichild0new;
			ndold2new[ichild1] = ichild1new;
		};
	};
	for (inodeold=0;inodeold<tree.nnodes;inodeold++) {
		inodenew = ndold2new[inodeold];
		_tree.nodes[inodenew].indbeg = tree.nodes[inodeold].indbeg;
		_tree.nodes[inodenew].indend = tree.nodes[inodeold].indend;
		_tree.nodes[inodenew].indbegtot = tree.nodes[inodeold].indbegtot;
		_tree.nodes[inodenew].indendtot = tree.nodes[inodeold].indendtot;
	};
	for (icpu=0;icpu<nproc;icpu++) {
		inodeold = tree.cpuidend[icpu];
		inodenew = ndold2new[inodeold];
		inodeshft = inodenew-(treearr[icpu].nnodes-1);
		ishft = ibegcpuarr[icpu];
		for (inode=0;inode<treearr[icpu].nnodes;inode++) {
			inodeloc = inode+inodeshft;
			_tree.nodes[inodeloc].indbeg = treearr[icpu].nodes[inode].indbeg+ishft;
			_tree.nodes[inodeloc].indend = treearr[icpu].nodes[inode].indend+ishft;
			_tree.nodes[inodeloc].indbegtot = treearr[icpu].nodes[inode].indbegtot+ishft;
			_tree.nodes[inodeloc].indendtot = treearr[icpu].nodes[inode].indendtot+ishft;
		};
	};

// Init blks partitionings and distribution

	_tree.PartitionTreeSchur (_nblks, _blks, _blk2cpu);
	_tree.StorePartitioning (_nblks, _blks, _blk2cpu);
	_tree.CpuListsForNodes ();

// Store the block sparsity

	CSMatrix *paspblk = _tree.GetAblktree ();

	*paspblk = _tree.MaximalBlockSparsityForSortedTree (_nblks);

// Free work arrays

	delete [] blkscpu;
	delete [] blksloc;
	delete [] iplev;
	delete [] iptrlev;
	delete [] ibegcpuarr;
	delete [] iendcpuarr;
	delete [] iarr;
	delete [] nlistarr;
	delete [] nzjaarr;
	delete [] mtrarrloc;
	delete [] iptrcpu;
	delete [] listblk;
	delete [] iacpu;
	delete [] iptr;
	delete [] imask;
	delete [] listloc;
	delete [] nlistcpu;
	delete [] plistcpu;
	delete [] treearr;
	delete [] ndold2new;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute new ND-type ordering and partioning via ParMeTiS
//========================================================================================
void CSMatrix::OrderNodeNDParMetis (CMPIComm &_comm, // Compute new ND-type ordering and partioning via ParMeTiS
												int *&_order, int &_nblks, int *&_blks) const {

	const char *funcname = "OrderNodeNDParMetis";

// Check that nproc and nparts are both the powers of 2

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

	int np = 1;

	int nlev = 0;

	while (np*2 <= nproc) {
		np *= 2;
		nlev++;
	};

	if (np != nproc) throw " CSMatrix::OrderNodeNDParMetis: the number of cpu's is not a power of 2 ";

// Allocate the memory that supports the levels computations

	int *blkstot;

	blkstot = new int [4*nlev+8];
	if (!blkstot) MemoryFail (funcname);

	int i;

	for (i=0;i<4*nlev+8;i++) blkstot[i] = 0;

	int **pporderlev, *nblkslev, **ppblkslev;

	pporderlev = new int * [nlev+2];
	if (!pporderlev) MemoryFail (funcname);
	nblkslev = new int [nlev+2];
	if (!nblkslev) MemoryFail (funcname);
	ppblkslev = new int * [nlev+2];
	if (!ppblkslev) MemoryFail (funcname);

	CMPIComm *commlev;

	commlev = new CMPIComm [nlev+1];
	if (!commlev) MemoryFail (funcname);

// Create the whole set of communicators

	int *listcpu;

	listcpu = new int [nproc+1];
	if (!listcpu) MemoryFail (funcname);

	commlev[0] = _comm;

	int ilev, nprocloc, myidloc, nproc1;
	CMPIComm comm1, comm2;

	for (ilev=1;ilev<=nlev;ilev++) {
		nprocloc = commlev[ilev-1].GetNproc ();
		myidloc = commlev[ilev-1].GetMyid ();
		nproc1 = nprocloc / 2;
		for (i=0;i<nproc1;i++) listcpu[i] = i;
		comm1 = commlev[ilev-1].CreateComm (nproc1, listcpu);
		for (i=nproc1;i<nproc;i++) listcpu[i-nproc1] = i;
		comm2 = commlev[ilev-1].CreateComm (nproc1, listcpu);
		if (myidloc < nproc1) {
			commlev[ilev] = comm1;
		} else {
			commlev[ilev] = comm2;
		};
	};

// Init zero level data

	int *blksloc;

	blksloc = new int [nproc+1];
	if (!blksloc) MemoryFail (funcname);

	for (i=0;i<nproc+1;i++) blksloc[i] = 0;

	blksloc[myid+1] = nlist;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksloc, blksloc);

	for (i=0;i<nproc;i++) blksloc[i+1] = blksloc[i]+blksloc[i+1];

	int indbeg = list[0];
	int indend = list[nlist-1];

	if (indbeg != blksloc[myid] || indend != blksloc[myid+1]-1) throw " CSMatrix::OrderNodeNDParMetis: wrong part of the matrix is considered";

	nblkslev[0] = nproc;
	ppblkslev[0] = blksloc;

	int ni = indend-indbeg+1;

// Main forward cycle over the levels

	int nblksloc, igroup, ibegloc, iendloc, ibegsub, iendsub, nloc, myid1, nlistloc, nz;
	int j, jj;
	int *plist, *pia, *pja, *ialoc;

	CSMatrix mtrlev;
	CSMatrix mtrdummy;

	mtrlev = *this;

	for (ilev=0;ilev<nlev;ilev++) {

// Compute local ordering and partitioning, and reorder the matrix indices

		mtrlev.SplitNDParMetis (commlev[ilev],
										pporderlev[ilev], nblksloc, blkstot+ilev*4);

// Detect the group number for the current processor on the current level

		nprocloc = commlev[ilev].GetNproc ();
		myidloc = commlev[ilev].GetMyid ();

		nproc1 = nprocloc / 2;

		if (myidloc < nproc1) {
			igroup = 0;
		} else {
			igroup = 1;
		};

// Determine local initial partitioning

		ibegloc = blkstot[ilev*4+igroup];
		iendloc = blkstot[ilev*4+igroup+1]-1;

		nloc = iendloc-ibegloc+1;

		ni = nloc / nproc1;

		blksloc = new int [nproc1+1];
		if (!blksloc) MemoryFail (funcname);

		for (i=0;i<nproc1;i++) blksloc[i] = i*ni;
		blksloc[nproc1] = nloc;

		nblkslev[ilev+1] = nproc1;
		ppblkslev[ilev+1] = blksloc;

// Exchange matrix data according to the partitioning

		myid1 = commlev[ilev+1].GetMyid ();
		ibegsub = blksloc[myid1];
		iendsub = blksloc[myid1+1]-1;

		ibegsub += ibegloc;
		iendsub += ibegloc;

		CSMatrix mtrnew;

		mtrnew = mtrlev.GetSubmatrixViaIntervals (commlev[ilev],
																1, &ibegsub, &iendsub);

// Filter and shift the rows/columns indices

		nlistloc = mtrnew.GetNlist ();
		plist = mtrnew.GetList ();
		pia = mtrnew.GetIa ();
		pja = mtrnew.GetJa ();

		ialoc = new int [nlistloc+1];
		if (!ialoc) MemoryFail (funcname);

		ialoc[0] = 0;
		nz = 0;

		for (i=0;i<nlistloc;i++) {
			plist[i] -= ibegloc;
			for (j=pia[i];j<pia[i+1];j++) {
				jj = pja[j];
				if (jj >= ibegloc && jj <= iendloc) {
					pja[nz] = jj-ibegloc;
					nz++;
				};
			};
			ialoc[i+1] = nz;
		};

		for (i=0;i<=nlistloc;i++) pia[i] = ialoc[i];

		mtrnew.SetM (nloc);
		mtrnew.SetN (nloc);
		mtrnew.SetNsupr (nloc);
		mtrnew.SetNsupcBase (nloc);
		mtrnew.SetNzja (nz);

		delete [] ialoc;

		mtrlev = mtrnew;

		mtrnew = mtrdummy;

	};

// Reorder by ND on current cpu the last submatrix

	nloc = mtrlev.GetN ();

	pporderlev[nlev] = new int [nloc];
	if (!pporderlev[nlev]) MemoryFail (funcname);
	ppblkslev[nlev+1] = new int [2];
	if (!ppblkslev[nlev+1]) MemoryFail (funcname);

	mtrlev.OrderPrfMtr (-1, pporderlev[nlev]);

	nblkslev[nlev+1] = 1;
	ppblkslev[nlev+1][0] = 0;
	ppblkslev[nlev+1][1] = nloc;

// Backward cycle with the ordering restoration

	if (nproc > 1) {

		for (ilev=nlev-1;ilev>=0;ilev--) {

			CSMatrix::RestoreOrderNd (ilev, nlev,
												commlev, 
												blkstot, pporderlev, nblkslev, ppblkslev);

		};

	};

	ni = indend-indbeg+1;

	_order = new int [ni];
	if (!_order) MemoryFail (funcname);

	for (i=0;i<ni;i++) _order[i] = pporderlev[0][i];

// Compute the final blks array

	_nblks = 2*nproc-1;

	_blks = new int [_nblks*4];
	if (!_blks) MemoryFail (funcname);

	for (i=0;i<_nblks*4;i++) _blks[i] = 0;

	for (i=0;i<nlev;i++) {
		myidloc = commlev[i].GetMyid ();
		if (myidloc != 0) {
			for (j=0;j<4;j++) {
				blkstot[i*4+j] = 0;
			};
		};
	};

	int inodesum = 0;

	int inode, np1;

	np = 1;

	for (i=0;i<nlev;i++) {

		np1 = nproc/np;

		IntCompDummy (np1);

		inode = inodesum + myid / np1;

		IntCompDummy (inode);

		for (j=0;j<4;j++) _blks[inode*4+j] = blkstot[i*4+j];

		inodesum = inodesum + np;

		IntCompDummy (inodesum);

		np = np * 2;

		IntCompDummy (np);

	};

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nblks*4, _blks, _blks);

// Free work arrays

	delete [] blkstot;

	for (i=0;i<=nlev;i++) {
		delete [] pporderlev[i];
		delete [] ppblkslev[i];
	};

	delete [] ppblkslev[nlev+1];
	delete [] pporderlev;
	delete [] nblkslev;
	delete [] ppblkslev;
	delete [] commlev;
	delete [] listcpu;

};

// Author: Kharchenko S.A.
// CSMatrix: Split the matrix into two parts and create the corresponding ordering via ParMeTiS
//========================================================================================
void CSMatrix::SplitNDParMetis (CMPIComm &_comm, // Split the matrix into two parts and create the corresponding ordering via ParMeTiS
											int *&_order, int &_nblks, int *_blks) {

	const char *funcname = "SplitNDParMetis";

// Collect sizes statistics for the whole set of matrices

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

	int *blkscpu;
	int *blk2cpu;

	blkscpu = new int [nproc+1];
	if (!blkscpu) MemoryFail (funcname);
	blk2cpu = new int [nproc];
	if (!blk2cpu) MemoryFail (funcname);

	int i;

	for (i=0;i<nproc+1;i++) blkscpu[i] = 0;
	for (i=0;i<nproc;i++) blk2cpu[i] = i;

	blkscpu[myid+1] = nlist;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blkscpu, blkscpu);

	for (i=0;i<nproc;i++) blkscpu[i+1] = blkscpu[i]+blkscpu[i+1];

// Check that matrix data are also correct

	int indbeg = list[0];
	int indend = list[nlist-1];

	if (indbeg != blkscpu[myid] || indend != blkscpu[myid+1]-1) throw " CSMatrix::SplitNDParMetis: wrong part of the matrix is considered";

//	int ntot = blkscpu[nproc];

// Compute symmetrized matrix in parallel

	CSMatrix as;

	as = SymmMtr (_comm, nproc, blkscpu, blk2cpu);

	if (false) {

// Check that sparsity is symmetric

		int ntotloc = blkscpu[nproc];

		int nploc = 0;
		int ibegloc = 0;
		int iendloc = -1;

		if (myid == 0) {
			nploc = 1;
			ibegloc = 0;
			iendloc = ntotloc-1;
		};

		CSMatrix mtrloc;

		mtrloc = as.GetSubmatrixViaIntervals (_comm,
																	nploc, &ibegloc, &iendloc);

		if (myid == 0) {

			int *plistloc, *pialoc, *pjaloc;

			plistloc = mtrloc.GetList ();
			pialoc = mtrloc.GetIa ();
			pjaloc = mtrloc.GetJa ();

			int j, jj, k, kk, icheck, jcheck;

			for (i=0;i<ntotloc;i++) {
				if (plistloc[i] != i) throw " CSMatrix::SplitNDParMetis: wrong array list ";
				for (j=pialoc[i];j<pialoc[i+1];j++) {
					jj = pjaloc[j];
					if (jj < 0 || jj >= ntotloc) {
						throw " CSMatrix::SplitNDParMetis: wrong sparsity ";
					};
				};
				icheck = -1;
				for (j=pialoc[i];j<pialoc[i+1];j++) {
					jj = pjaloc[j];
					if (jj == i) icheck = 0;
					jcheck = -1;
					for (k=pialoc[jj];k<pialoc[jj+1];k++) {
						kk = pjaloc[k];
						if (kk == i) jcheck = 0;
					};
					if (jcheck != 0) {
						throw " CSMatrix::SplitNDParMetis: matrix is not symmetric ";
					};
				};
				if (icheck != 0) {
					throw " CSMatrix::SplitNDParMetis: diagonal not found ";
				};
			};
		};

	};

// Eliminate diagonal

	CSMatrix asfltd = as;

	int j, jj, irow;

	int nzjaloc = 0;

	asfltd.ia[0] = 0;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			jj = as.ja[j];
			if (irow != jj) {
				asfltd.ja[nzjaloc] = jj;
				nzjaloc++;
			};
		};
		asfltd.ia[i+1] = nzjaloc;
	};

	asfltd.nzja = nzjaloc;
	asfltd.nsupc = asfltd.n;
	asfltd.nsupr = asfltd.n;

// Allocate the memory

	int *partition;

	partition = new int [nlist];
	if (!partition) MemoryFail (funcname);

// Call ParMETIS

	int options[10];
	options[0] = 0;

	int numflag = 0;

	int ncon = 1;
	int nparts = 2;

	float *tpwgts, *ubvec;

	tpwgts = new float [ncon*nparts];
	if (!tpwgts) MemoryFail (funcname);
	ubvec = new float [ncon];
	if (!ubvec) MemoryFail (funcname);

	for (i=0;i<nparts;i++) tpwgts[i] = (float)1.0 / ((float)nparts);

	ubvec[0] = (float) 1.05;

	int edgecut;

   CParMETIS::PartKway (_comm, blkscpu, asfltd.ia, asfltd.ja, 
                        NULL, NULL, 0,
                        &numflag, &ncon, &nparts,
                        tpwgts, ubvec,
                        options, &edgecut, partition);

	CSMatrix mtrdummy;

	asfltd = mtrdummy;

// Determine all separator nodes

// Mark the number of processor for each connection

	int nzjatot = as.GetNzja ();

	int *ja2cpu;

	ja2cpu = new int [nzjatot];
	if (!ja2cpu) MemoryFail (funcname);

	int icpubeg, icpuend, icpu;

	for (j=0;j<nzjatot;j++) {
		jj = as.ja[j];
		if (jj >= blkscpu[myid] && jj < blkscpu[myid+1]) {
			icpubeg = myid;
			icpuend = myid;
		} else {
			icpubeg = 0;
			icpuend = nproc-1;
			while (icpubeg != icpuend) {
				icpu = (icpuend+icpubeg) / 2;
				if (icpu < icpubeg) icpu = icpubeg;
				if (icpu > icpuend) icpu = icpuend;
				if (jj >= blkscpu[icpu] && jj < blkscpu[icpu+1]) {
					icpubeg = icpu;
					icpuend = icpu;
				} else if (jj < blkscpu[icpu]) {
					icpuend = icpu-1;
					if (icpuend < icpubeg) icpuend = icpubeg;
				} else if (jj >= blkscpu[icpu+1]) {
					icpubeg = icpu+1;
					if (icpubeg > icpuend) icpubeg = icpuend;
				};
			};
		};
		ja2cpu[j] = icpubeg;
	};

// Create cpu substructures

	int *iacpu;
	int *iptr;
	int *jaindcpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);
	jaindcpu = new int [nzjatot];
	if (!jaindcpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	for (j=0;j<nzjatot;j++) {
		icpu = ja2cpu[j];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int k;

	for (j=0;j<nzjatot;j++) {
		icpu = ja2cpu[j];
		k = iptr[icpu];
		jaindcpu[k] = j;
		iptr[icpu]++;
	};

	delete [] ja2cpu;

// Create and init global imask

	int nlistmx = 0;

	int ni;

	for (i=0;i<nproc;i++) {
		ni = blkscpu[i+1]-blkscpu[i];
		if (ni > nlistmx) nlistmx = ni;
	};

	int *imask;
	int *listloc;

	imask = new int [nlistmx];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nlistmx];
	if (!listloc) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nlistmx;i++) {
		imask[i] = icycle;
	};

// Create the lists of indices for each cpu

	int *nlistcpu;
	int **plistcpu;

	nlistcpu = new int [nproc];
	if (!nlistcpu) MemoryFail (funcname);
	plistcpu = new int * [nproc*2];
	if (!plistcpu) MemoryFail (funcname);

	int ind, ishft, jjloc, nz;

	for (i=0;i<nproc;i++) {
		icycle++;
		ishft = blkscpu[i];
		nz = 0;
		for (j=iacpu[i];j<iacpu[i+1];j++) {
			ind = jaindcpu[j];
			jj = as.ja[ind];
			jjloc = jj-ishft;
			if (imask[jjloc] != icycle) {
				listloc[nz] = jj;
				nz++;
				imask[jjloc] = icycle;
			};
		};
		nlistcpu[i] = nz;
		plistcpu[i] = new int [nz];
		if (!plistcpu[i]) MemoryFail (funcname);
		plistcpu[nproc+i] = new int [nz];
		if (!plistcpu[nproc+i]) MemoryFail (funcname);
		for (j=0;j<nz;j++) plistcpu[i][j] = listloc[j];
	};

// Exchange lists

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = nlistcpu[i]*sizeof(int);
		ObjSend[i] = (char *) plistcpu[i];
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
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	ishft = blkscpu[myid];

	int jjnew;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / sizeof(int);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[j];
			jjloc = jj-ishft;
			jjnew = partition[jjloc];
			piarr[j] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store new indices numbers

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / sizeof(int);
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			plistcpu[nproc+icpu][j] = piarr[j];
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// For each sparsity element compute its partitioning value

	int *ja2part;

	ja2part = new int [nzjatot];
	if(!ja2part) MemoryFail(funcname);

	for (icpu=0;icpu<nproc;icpu++) {

		ishft = blkscpu[icpu];

// Init list data

		ni = nlistcpu[icpu];

		for (i=0;i<ni;i++) {
			jj = plistcpu[icpu][i];
			jjloc = jj-ishft;
			jjnew = plistcpu[icpu+nproc][i];
			listloc[jjloc] = jjnew;
		};

// Repartition

		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ind = jaindcpu[j];
			jj = as.ja[ind];
			jjloc = jj-ishft;
			jjnew = listloc[jjloc];
			ja2part[ind] = jjnew;
		};

	};

// Mark each node of the submatrix ("0","1" - partitions, "-1" - separator)

	int *nd2part;

	nd2part = new int [nlist];
	if(!nd2part) MemoryFail(funcname);

	int ivalue, ipart;

	for (i=0;i<nlist;i++) {
		ivalue = partition[i];
		ipart = ivalue;
		for (j=as.ia[i];j<as.ia[i+1];j++) {
			if (ja2part[j] != ivalue) ipart = -1;
		};
		nd2part[i] = ipart;
	};

// Initial the number for each portion

	int npart0 = 0;
	int npart1 = 0;
	int nsep = 0;

	for (i=0;i<nlist;i++) {
		if (nd2part[i] == 0) npart0++;
		if (nd2part[i] == 1) npart1++;
		if (nd2part[i] ==-1) nsep++;
	};

//	cout << " Initial npart0 = " << npart0 << " npart1 = " << npart1 << " nsep = " << nsep << endl;

// Modify if necessary the nd2part array to minimize the bordering

	if (true) {

// Cycle over the partitions

		int ip, ip1;

		for (ip=0;ip<2;ip++) {

// Collect all necessary current ipart values from other processors

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

			for(i = 0; i < NObjSend; i++) {
				ObjTypeSend[i] = 1;
				ObjIDSend[i] = i;
				CpuIDSend[i] = i;
				ObjSizeSend[i] = nlistcpu[i]*sizeof(int);
				ObjSend[i] = (char *) plistcpu[i];
			};

			info = CMPIExchange::DataExchangeMPI (_comm,
																NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
																ObjSizeSend, ObjSend,
																NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
																ObjSizeRecv, ObjRecv);
			if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

			delete [] ObjTypeSend;
			delete [] ObjIDSend;
			delete [] CpuIDSend;
			delete [] ObjSizeSend;
			delete [] ObjSend;

// Prepare reply

			ishft = blkscpu[myid];

			for (i=0;i<NObjRecv;i++) {
				ni = ObjSizeRecv[i] / sizeof(int);
				piarr = (int *) ObjRecv[i];
				for (j=0;j<ni;j++) {
					jj = piarr[j];
					jjloc = jj-ishft;
					jjnew = nd2part[jjloc];
					piarr[j] = jjnew;
				};
			};

// Send indices back

			info = CMPIExchange::DataExchangeMPI (_comm,
																NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
																ObjSizeRecv, ObjRecv,
																NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
																ObjSizeSend, ObjSend);
			if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

			delete [] ObjTypeRecv;
			delete [] ObjIDRecv;
			delete [] CpuIDRecv;
			delete [] ObjSizeRecv;
			for(i = 0; i < NObjRecv; i++) {
				delete [] ObjRecv[i];
			};
			delete [] ObjRecv;

// Store new indices numbers

			for (i=0;i<NObjSend;i++) {
				icpu = CpuIDSend[i];
				ni = ObjSizeSend[i] / sizeof(int);
				piarr = (int *) ObjSend[i];
				for (j=0;j<ni;j++) {
					plistcpu[nproc+icpu][j] = piarr[j];
				};
			};

// Free receive structures

			delete [] ObjTypeSend;
			delete [] ObjIDSend;
			delete [] CpuIDSend;
			delete [] ObjSizeSend;
			for(i = 0; i < NObjSend; i++) {
				delete [] ObjSend[i];
			};
			delete [] ObjSend;

// Scan the matrix once again and modify current partition value for each connection

			for (icpu=0;icpu<nproc;icpu++) {

				ishft = blkscpu[icpu];

// Init list data

				ni = nlistcpu[icpu];

				for (i=0;i<ni;i++) {
					jj = plistcpu[icpu][i];
					jjloc = jj-ishft;
					jjnew = plistcpu[icpu+nproc][i];
					listloc[jjloc] = jjnew;
				};

// Repartition

				for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
					ind = jaindcpu[j];
					jj = as.ja[ind];
					jjloc = jj-ishft;
					jjnew = listloc[jjloc];
					ja2part[ind] = jjnew;
				};

			};

// Modify nd2part if possible

			ip1 = (ip+1)%2;

			if (ip == 0) {

				int ncount = 0;

				for (i=0;i<nlist;i++) {
					ivalue = partition[i];
					if (ivalue == ip && nd2part[i] == -1) {
						ipart = ivalue;
						for (j=as.ia[i];j<as.ia[i+1];j++) {
							if (ja2part[j] == ip1) ipart = -1;
						};
						if (ipart >= 0) ncount++;
					};
				};

				int ncountmax = ncount / 2;

				ncount = 0;

				for (i=0;i<nlist;i++) {
					ivalue = partition[i];
					if (ivalue == ip && nd2part[i] == -1 && ncount < ncountmax) {
						ipart = ivalue;
						for (j=as.ia[i];j<as.ia[i+1];j++) {
							if (ja2part[j] == ip1) ipart = -1;
						};
						if (ipart >= 0) {
							nd2part[i] = ipart;
							ncount++;
						};
					};
				};

			} else {

				for (i=0;i<nlist;i++) {
					ivalue = partition[i];
					if (ivalue == ip && nd2part[i] == -1) {
						ipart = ivalue;
						for (j=as.ia[i];j<as.ia[i+1];j++) {
							if (ja2part[j] == ip1) ipart = -1;
						};
						nd2part[i] = ipart;
					};
				};

			};

		};

	};

// Count the number for each portion

	npart0 = 0;
	npart1 = 0;
	nsep = 0;

	for (i=0;i<nlist;i++) {
		if (nd2part[i] == 0) npart0++;
		if (nd2part[i] == 1) npart1++;
		if (nd2part[i] ==-1) nsep++;
	};

//	cout << " Final npart0 = " << npart0 << " npart1 = " << npart1 << " nsep = " << nsep << endl;

// Create arrays that support ordering computations

	int *blkssep;

	blkssep = new int [nproc*3+1];
	if(!blkssep) MemoryFail(funcname);

	for (i=0;i<=nproc*3;i++) blkssep[i] = 0;

	blkssep[myid+1] = npart0;
	blkssep[nproc+myid+1] = npart1;
	blkssep[2*nproc+myid+1] = nsep;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, 3*nproc+1, blkssep, blkssep);

	for (i=0;i<nproc*3;i++) blkssep[i+1] = blkssep[i]+blkssep[i+1];

// Create ordering

	int ibs0 = blkssep[myid];
	int ibs1 = blkssep[nproc+myid];
	int ibssep = blkssep[2*nproc+myid];

	_order = new int [nlist];
	if(!_order) MemoryFail(funcname);

	for (i=0;i<nlist;i++) {
		if (nd2part[i] == 0) {
			_order[i] = ibs0;
			ibs0++;
		};
		if (nd2part[i] == 1) {
			_order[i] = ibs1;
			ibs1++;
		};
		if (nd2part[i] == -1) {
			_order[i] = ibssep;
			ibssep++;
		};
	};

// Assign blks array

	int npart0tot = blkssep[nproc];
	int npart1tot = blkssep[2*nproc]-blkssep[nproc];
	int nseptot = blkssep[3*nproc]-blkssep[2*nproc];

	_nblks = 3;
	_blks[0] = 0;
	_blks[1] = npart0tot;
	_blks[2] = npart0tot+npart1tot;
	_blks[3] = npart0tot+npart1tot+nseptot;

// Renumber all rows/column indices

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = nlistcpu[i]*sizeof(int);
		ObjSend[i] = (char *) plistcpu[i];
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	ishft = blkscpu[myid];

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / sizeof(int);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[j];
			jjloc = jj-ishft;
			jjnew = _order[jjloc];
			piarr[j] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store new indices numbers

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / sizeof(int);
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			plistcpu[nproc+icpu][j] = piarr[j];
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Renumber

	for (icpu=0;icpu<nproc;icpu++) {

		ishft = blkscpu[icpu];

// Init list data

		ni = nlistcpu[icpu];

		for (i=0;i<ni;i++) {
			jj = plistcpu[icpu][i];
			jjloc = jj-ishft;
			jjnew = plistcpu[icpu+nproc][i];
			listloc[jjloc] = jjnew;
		};

// Renumber the columns

		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ind = jaindcpu[j];
			jj = as.ja[ind];
			jjloc = jj-ishft;
			jjnew = listloc[jjloc];
			as.ja[ind] = jjnew;
		};

	};

	for (i=0;i<nlist;i++) as.list[i] = _order[i];

	*this = as;

// Sort matrix data

	SortListAndColumns ();

// Free work arrays

	delete [] partition;
	delete [] blkscpu;
	delete [] blk2cpu;
	delete [] tpwgts;
	delete [] ubvec;
	delete [] iacpu;
	delete [] iptr;
	delete [] jaindcpu;
	delete [] imask;
	delete [] listloc;
	delete [] nlistcpu;

	for (i=0;i<2*nproc;i++) {
		delete [] plistcpu[i];
	};
	delete [] plistcpu;

	delete [] ja2part;
	delete [] nd2part;
	delete [] blkssep;

};

// Author: Kharchenko S.A.
// CSMatrix: Get the submatrix from other processors described by the list of intervals
//========================================================================================
CSMatrix CSMatrix::GetSubmatrixViaIntervals (CMPIComm &_comm, // Get the matrix part from other processors described by the list of intervals
															int _nparts, int *_ibegarr, int *_iendarr) {

	const char *funcname = "GetSubmatrixViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	int i;

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::GetSubmatrixViaIntervals: intervals from different cpu's intersect ";
		};
	};

// For local list data find to which processor the data should be sent

	int *list2cpu;

	list2cpu = new int [nlist];
	if(!list2cpu) MemoryFail(funcname);

	for (i=0;i<nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ip, irow, iproc;

	if (nptot > 0) {

		ipbeg = 0;

		for (i=0;i<nlist;i++) {
			irow = list[i];
			ipend = nptot-1;
			ip = ipbeg;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[i] = icputot[ip];
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						list2cpu[i] = icputot[ip];
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;
	int *nzja2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);
	nzja2cpu = new int [nproc];
	if(!nzja2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;
	for (i=0;i<nproc;i++) nzja2cpu[i] = 0;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
			nzja2cpu[iproc] += (ia[i+1]-ia[i]);
		};
	};

	CSMatrix *mtrarr;

	mtrarr = new CSMatrix [nproc];
	if(!mtrarr) MemoryFail(funcname);

	for (i=0;i<nproc;i++) {
		CSMatrix temp (nlist2cpu[i],nzja2cpu[i]);
		mtrarr[i] = temp;
	};

	int *nl2cpu;
	int *nz2cpu;

	nl2cpu = new int [nproc];
	if(!nl2cpu) MemoryFail(funcname);
	nz2cpu = new int [nproc];
	if(!nz2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nl2cpu[i] = 0;
	for (i=0;i<nproc;i++) nz2cpu[i] = 0;

	int *plist, *pia, *pja;

	for (i=0;i<nproc;i++) {
		pia = mtrarr[i].GetIa ();
		pia[0] = 0;
	};

	int ilist, nzloc;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			plist = mtrarr[iproc].GetList ();
			pia = mtrarr[iproc].GetIa ();
			pja = mtrarr[iproc].GetJa ();
			ilist = nl2cpu[iproc];
			nzloc = nz2cpu[iproc];
			plist[ilist] = list[i];
			for (j=ia[i];j<ia[i+1];j++) {
				pja[nzloc] = ja[j];
				nzloc++;
			};
			pia[ilist+1] = nzloc;
			nl2cpu[iproc]++;
			nz2cpu[iproc] = nzloc;
		};
	};

	for (i=0;i<nproc;i++) {
		mtrarr[i].SetM (m);
		mtrarr[i].SetN (n);
		mtrarr[i].SetNsupr (nsupr);
		mtrarr[i].SetNsupcBase (nsupc);
		mtrarr[i].SetNlist (nlist2cpu[i]);
		mtrarr[i].SetNzja (nzja2cpu[i]);
	};

// Exchange the data

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
		mtrarr[i].PackMatrix (ObjSizeSend[i],ObjSend[i]);
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
	if(info) throw " CSMatrix::GetSubmatrixViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

	if (NObjRecv != nproc) throw " CSMatrix::GetSubmatrixViaIntervals: wrong number of the received objects ";

	for(i = 0; i < NObjRecv; i++) {
		mtrarr[i].UnPackMatrix (ObjSizeRecv[i],ObjRecv[i]);
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Collect submatrices of the current cpu

	int *listblk;

	listblk = new int [nproc];
	if (!listblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) listblk[i] = i;

	CSMatrix mtrnew;

	CGSMatrix::CombineMatricesSort (n, nproc, listblk,
												mtrarr,
												mtrnew);

	int nloc = mtrnew.GetN ();

	mtrnew.SetNsupr (nloc);
	mtrnew.SetNsupcBase (nloc);

	CSMatrix mtrdummy;

	for (i=0;i<nproc;i++) mtrarr[i] = mtrdummy;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] list2cpu;
	delete [] nlist2cpu;
	delete [] nzja2cpu;
	delete [] mtrarr;
	delete [] nl2cpu;
	delete [] nz2cpu;
	delete [] listblk;

// Return the result

	return mtrnew;

};

// Author: Kharchenko S.A.
// CSMatrix: Get the submatrix from other processors described by the list of intervals (2index version)
//========================================================================================
CSMatrix CSMatrix::GetSubmatrixViaIntervals2Index (CMPIComm &_comm, // Get the matrix part from other processors described by the list of intervals (2index version)
																	int _nparts, int *_ibegarr, int *_iendarr) {

	const char *funcname = "GetSubmatrixViaIntervals2Index";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	int i;

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::GetSubmatrixViaIntervals: intervals from different cpu's intersect ";
		};
	};

// For local list data find to which processor the data should be sent

	int *list2cpu;

	list2cpu = new int [nlist];
	if(!list2cpu) MemoryFail(funcname);

	for (i=0;i<nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ip, irow, iproc;

	if (nptot > 0) {

		ipbeg = 0;

		for (i=0;i<nlist;i++) {
			irow = list2[i];
			ipend = nptot-1;
			ip = ipbeg;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[i] = icputot[ip];
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						list2cpu[i] = icputot[ip];
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;
	int *nzja2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);
	nzja2cpu = new int [nproc];
	if(!nzja2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;
	for (i=0;i<nproc;i++) nzja2cpu[i] = 0;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
			nzja2cpu[iproc] += (ia[i+1]-ia[i]);
		};
	};

	CSMatrix *mtrarr;

	mtrarr = new CSMatrix [nproc];
	if(!mtrarr) MemoryFail(funcname);

	for (i=0;i<nproc;i++) {
		CSMatrix temp (nlist2cpu[i],nlist2cpu[i],nzja2cpu[i],nzja2cpu[i]);
		mtrarr[i] = temp;
	};

	int *nl2cpu;
	int *nz2cpu;

	nl2cpu = new int [nproc];
	if(!nl2cpu) MemoryFail(funcname);
	nz2cpu = new int [nproc];
	if(!nz2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nl2cpu[i] = 0;
	for (i=0;i<nproc;i++) nz2cpu[i] = 0;

	int *plist, *pia, *pja;

	for (i=0;i<nproc;i++) {
		pia = mtrarr[i].GetIa ();
		pia[0] = 0;
	};

	int ilist, nzloc;
	int *plist2, *pja2;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			plist = mtrarr[iproc].GetList ();
			plist2 = mtrarr[iproc].GetList2 ();
			pia = mtrarr[iproc].GetIa ();
			pja = mtrarr[iproc].GetJa ();
			pja2 = mtrarr[iproc].GetJa2 ();
			ilist = nl2cpu[iproc];
			nzloc = nz2cpu[iproc];
			plist[ilist] = list[i];
			plist2[ilist] = list2[i];
			for (j=ia[i];j<ia[i+1];j++) {
				pja[nzloc] = ja[j];
				pja2[nzloc] = ja2[j];
				nzloc++;
			};
			pia[ilist+1] = nzloc;
			nl2cpu[iproc]++;
			nz2cpu[iproc] = nzloc;
		};
	};

	for (i=0;i<nproc;i++) {
		mtrarr[i].SetM (m);
		mtrarr[i].SetN (n);
		mtrarr[i].SetNsupr (nsupr);
		mtrarr[i].SetNsupcBase (nsupc);
		mtrarr[i].SetNlist (nlist2cpu[i]);
		mtrarr[i].SetNlist2 (nlist2cpu[i]);
		mtrarr[i].SetNzja (nzja2cpu[i]);
		mtrarr[i].SetNzja2 (nzja2cpu[i]);
	};

// Exchange the data

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
		mtrarr[i].PackMatrix (ObjSizeSend[i],ObjSend[i]);
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
	if(info) throw " CSMatrix::GetSubmatrixViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

	if (NObjRecv != nproc) throw " CSMatrix::GetSubmatrixViaIntervals: wrong number of the received objects ";

	for(i = 0; i < NObjRecv; i++) {
		mtrarr[i].UnPackMatrix (ObjSizeRecv[i],ObjRecv[i]);
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Collect submatrices of the current cpu

	int *listblk;

	listblk = new int [nproc];
	if (!listblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) listblk[i] = i;

	CSMatrix mtrnew;

	CGSMatrix::CombineMatricesSort2Index (n, nproc, listblk,
														mtrarr,
														mtrnew);

	int nloc = mtrnew.GetN ();

	mtrnew.SetNsupr (nloc);
	mtrnew.SetNsupcBase (nloc);

	CSMatrix mtrdummy;

	for (i=0;i<nproc;i++) mtrarr[i] = mtrdummy;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] list2cpu;
	delete [] nlist2cpu;
	delete [] nzja2cpu;
	delete [] mtrarr;
	delete [] nl2cpu;
	delete [] nz2cpu;
	delete [] listblk;

// Return the result

	return mtrnew;

};

// Author: Kharchenko S.A.
// CSMatrix: Restore local ordering at the previous level
//========================================================================================
void CSMatrix::RestoreOrderNd (int _ilev, int _nlev, // Restore local ordering at the previous level
											CMPIComm *_commlev, 
											int *_blks, int **_pporderlev, int *_nblkslev, int **_ppblkslev) {

	const char *funcname = "RestoreOrderNd";

// Shift local indices

	int nprocloc = _commlev[_ilev].GetNproc ();
	int myidloc  = _commlev[_ilev].GetMyid ();
	int nproc1   = _commlev[_ilev+1].GetNproc ();
	int myidloc1 = _commlev[_ilev+1].GetMyid ();

	int igroup;

	if (myidloc < nproc1) {
		igroup = 0;
	} else {
		igroup = 1;
	};

	int ishft = _blks[_ilev*4+igroup];
	int niloc1 = _ppblkslev[_ilev+1][myidloc1+1]-_ppblkslev[_ilev+1][myidloc1];

	int i;

// Determine the whole set of intervals

	int *blksind;

	blksind = new int [nprocloc+2];
	if(!blksind) MemoryFail(funcname);

	for (i=0;i<nprocloc+2;i++) blksind[i] = 0;

	blksind[myidloc+1] = niloc1;

	if (myidloc == 0) {
		if (_ilev < _nlev) {
			blksind[nprocloc+1] = _blks[_ilev*4+3]-_blks[_ilev*4+2];
		};
	};

	CMPIExchange::ExchangeArrayMPI (_commlev[_ilev], INTEGERVALUE, ADD, nprocloc+2, blksind, blksind);

	for (i=0;i<nprocloc+1;i++) blksind[i+1] = blksind[i]+blksind[i+1];

// Find interval numbers of each index

	int indbeg = _ppblkslev[_ilev][myidloc];
	int indend = _ppblkslev[_ilev][myidloc+1]-1;
	int niloc = indend-indbeg+1;

	int *order2cpu;

	order2cpu = new int [niloc];
	if(!order2cpu) MemoryFail(funcname);

	int irow, iblk, ipbeg, ipend, ip;

	for (i=0;i<niloc;i++) {
		irow = _pporderlev[_ilev][i];
		iblk = -1;
		ipbeg = 0;
		ipend = nprocloc+1;
		while (true) {
			if (ipbeg > ipend) break;
			ip = (ipbeg+ipend) / 2;
			if (ip < ipbeg) ip = ipbeg;
			if (ip > ipend) ip = ipend;
			if (irow >= blksind[ip] && irow <= blksind[ip+1]-1) {
				iblk = ip;
				break;
			};
			if (irow < blksind[ip]) ipend = ip-1;
			if (irow >= blksind[ip+1]) ipbeg = ip+1;
		};
		if (iblk >= nprocloc) iblk = -1;
		order2cpu[i] = iblk;
	};

// Prepare the lists for sends

	int *nlist2cpu;
	int **plist2cpu;
	int *iptr2cpu;

	nlist2cpu = new int [nprocloc];
	if(!nlist2cpu) MemoryFail(funcname);
	plist2cpu = new int * [nprocloc];
	if(!plist2cpu) MemoryFail(funcname);
	iptr2cpu = new int [nprocloc];
	if(!iptr2cpu) MemoryFail(funcname);

	for (i=0;i<nprocloc;i++) nlist2cpu[i] = 0;

	int iproc;

	for (i=0;i<niloc;i++) {
		iproc = order2cpu[i];
		if (iproc >= 0) nlist2cpu[iproc]++;
	};

	for (i=0;i<nprocloc;i++) {
		plist2cpu[i] = new int [2*nlist2cpu[i]];
		if(!plist2cpu[i]) MemoryFail(funcname);
	};

	for (i=0;i<nprocloc;i++) iptr2cpu[i] = 0;

	int k, ishftloc;

	for (i=0;i<niloc;i++) {
		iproc = order2cpu[i];
		if (iproc >= 0) {
			k = iptr2cpu[iproc];
			plist2cpu[iproc][k] = i;
			ishftloc = 0;
			if (iproc >= nproc1) ishftloc = _blks[_ilev*4+1];
			plist2cpu[iproc][nlist2cpu[iproc]+k] = _pporderlev[_ilev][i]-ishftloc;
			iptr2cpu[iproc]++;
		};
	};

// Exchange the data

	int NObjSend;
	int* ObjTypeSend;
	int* ObjIDSend;
	int* CpuIDSend;
	int* ObjSizeSend;
	char** ObjSend;

	NObjSend = nprocloc;

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = 2*nlist2cpu[i]*sizeof(int);
		ObjSend[i] = (char *) plist2cpu[i];
	};

	int NObjRecv;
	int* ObjTypeRecv;
	int* ObjIDRecv;
	int* CpuIDRecv;
	int* ObjSizeRecv;
	char** ObjRecv;

	int info;

	info = CMPIExchange::DataExchangeMPI (_commlev[_ilev],
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::RestoreOrderNd: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	int jshft = _ppblkslev[_ilev+1][myidloc1];

	int ni, j, irowloc, irownew;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (2*sizeof(int));
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			irow = piarr[ni+j];
			irowloc = irow - jshft;
			irownew = _pporderlev[_ilev+1][irowloc];
			piarr[ni+j] = irownew+ishft;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_commlev[_ilev],
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::RestoreOrderNd: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store new ordering

	int ilist;

	for (i=0;i<NObjSend;i++) {
		ni = ObjSizeSend[i] / (2*sizeof(int));
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			ilist = piarr[j];
			irow = piarr[ni+j];
			_pporderlev[_ilev][ilist] = irow;
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

// Free work arrays

	delete [] blksind;
	delete [] order2cpu;
	delete [] nlist2cpu;
	for (i=0;i<nprocloc;i++) {
		delete [] plist2cpu[i];
	};
	delete [] plist2cpu;
	delete [] iptr2cpu;

};

// Author: Kharchenko S.A.
// CSMatrixR: Change rows/column indices
//========================================================================================
void CSMatrixR::ChangeIndices (CMPIComm &_comm, int *_order) { // Change rows/column indices

	const char *funcname = "ChangeIndices";

// Collect sizes statistics for the whole set of matrices

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

	int *blkscpu;

	blkscpu = new int [nproc+1];
	if (!blkscpu) MemoryFail (funcname);

	int i;

	for (i=0;i<nproc+1;i++) blkscpu[i] = 0;

	blkscpu[myid+1] = nlist;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blkscpu, blkscpu);

	for (i=0;i<nproc;i++) blkscpu[i+1] = blkscpu[i]+blkscpu[i+1];

// Check that matrix data are also correct

	int indbeg = list[0];
	int indend = list[nlist-1];

	if (indbeg != blkscpu[myid] || indend != blkscpu[myid+1]-1) throw " CSMatrix::ChangeIndices: wrong part of the matrix is considered";

//	int ntot = blkscpu[nproc];

// Mark the number of processor for each connection

	int nzjatot = GetNzja ();

	int *ja2cpu;

	ja2cpu = new int [nzjatot];
	if (!ja2cpu) MemoryFail (funcname);

	int icpubeg, icpuend, icpu;
	int j, jj;

	for (j=0;j<nzjatot;j++) {
		jj = ja[j];
		if (jj >= blkscpu[myid] && jj < blkscpu[myid+1]) {
			icpubeg = myid;
			icpuend = myid;
		} else {
			icpubeg = 0;
			icpuend = nproc-1;
			while (icpubeg != icpuend) {
				icpu = (icpuend+icpubeg) / 2;
				if (icpu < icpubeg) icpu = icpubeg;
				if (icpu > icpuend) icpu = icpuend;
				if (jj >= blkscpu[icpu] && jj < blkscpu[icpu+1]) {
					icpubeg = icpu;
					icpuend = icpu;
				} else if (jj < blkscpu[icpu]) {
					icpuend = icpu-1;
					if (icpuend < icpubeg) icpuend = icpubeg;
				} else if (jj >= blkscpu[icpu+1]) {
					icpubeg = icpu+1;
					if (icpubeg > icpuend) icpubeg = icpuend;
				};
			};
		};
		ja2cpu[j] = icpubeg;
	};

// Create cpu substructures

	int *iacpu;
	int *iptr;
	int *jaindcpu;

	iacpu = new int [nproc+1];
	if (!iacpu) MemoryFail (funcname);
	iptr = new int [nproc];
	if (!iptr) MemoryFail (funcname);
	jaindcpu = new int [nzjatot];
	if (!jaindcpu) MemoryFail (funcname);

	for (i=0;i<=nproc;i++) iacpu[i] = 0;

	for (j=0;j<nzjatot;j++) {
		icpu = ja2cpu[j];
		iacpu[icpu+1]++;
	};

	for (i=0;i<nproc;i++) iacpu[i+1] = iacpu[i]+iacpu[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iacpu[i];

	int k;

	for (j=0;j<nzjatot;j++) {
		icpu = ja2cpu[j];
		k = iptr[icpu];
		jaindcpu[k] = j;
		iptr[icpu]++;
	};

	delete [] ja2cpu;

// Create and init global imask

	int nlistmx = 0;

	int ni;

	for (i=0;i<nproc;i++) {
		ni = blkscpu[i+1]-blkscpu[i];
		if (ni > nlistmx) nlistmx = ni;
	};

	int *imask;
	int *listloc;

	imask = new int [nlistmx];
	if (!imask) MemoryFail (funcname);
	listloc = new int [nlistmx];
	if (!listloc) MemoryFail (funcname);

	int icycle = -1;

	for (i=0;i<nlistmx;i++) {
		imask[i] = icycle;
	};

// Create the lists of indices for each cpu

	int *nlistcpu;
	int **plistcpu;

	nlistcpu = new int [nproc];
	if (!nlistcpu) MemoryFail (funcname);
	plistcpu = new int * [nproc*2];
	if (!plistcpu) MemoryFail (funcname);

	int ind, ishft, jjloc, nz;

	for (i=0;i<nproc;i++) {
		icycle++;
		ishft = blkscpu[i];
		nz = 0;
		for (j=iacpu[i];j<iacpu[i+1];j++) {
			ind = jaindcpu[j];
			jj = ja[ind];
			jjloc = jj-ishft;
			if (imask[jjloc] != icycle) {
				listloc[nz] = jj;
				nz++;
				imask[jjloc] = icycle;
			};
		};
		nlistcpu[i] = nz;
		plistcpu[i] = new int [nz];
		if (!plistcpu[i]) MemoryFail (funcname);
		plistcpu[nproc+i] = new int [nz];
		if (!plistcpu[nproc+i]) MemoryFail (funcname);
		for (j=0;j<nz;j++) plistcpu[i][j] = listloc[j];
	};

// Exchange lists

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

	for(i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = nlistcpu[i]*sizeof(int);
		ObjSend[i] = (char *) plistcpu[i];
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
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

// Prepare reply

	ishft = blkscpu[myid];

	int jjnew;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / sizeof(int);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			jj = piarr[j];
			jjloc = jj-ishft;
			jjnew = _order[jjloc];
			piarr[j] = jjnew;
		};
	};

// Send indices back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::OrderNDParMetis: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store new indices numbers

	for (i=0;i<NObjSend;i++) {
		icpu = CpuIDSend[i];
		ni = ObjSizeSend[i] / sizeof(int);
		piarr = (int *) ObjSend[i];
		for (j=0;j<ni;j++) {
			plistcpu[nproc+icpu][j] = piarr[j];
		};
	};

// Free receive structures

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Renumber

	for (icpu=0;icpu<nproc;icpu++) {

		ishft = blkscpu[icpu];

// Init list data

		ni = nlistcpu[icpu];

		for (i=0;i<ni;i++) {
			jj = plistcpu[icpu][i];
			jjloc = jj-ishft;
			jjnew = plistcpu[icpu+nproc][i];
			listloc[jjloc] = jjnew;
		};

// Renumber the columns

		for (j=iacpu[icpu];j<iacpu[icpu+1];j++) {
			ind = jaindcpu[j];
			jj = ja[ind];
			jjloc = jj-ishft;
			jjnew = listloc[jjloc];
			ja[ind] = jjnew;
		};

	};

	for (i=0;i<nlist;i++) list[i] = _order[i];

// Sort matrix data

	SortListAndColumns ();

// Free work arrays

	delete [] blkscpu;
	delete [] iacpu;
	delete [] iptr;
	delete [] jaindcpu;
	delete [] imask;
	delete [] listloc;
	delete [] nlistcpu;

	for (i=0;i<2*nproc;i++) {
		delete [] plistcpu[i];
	};
	delete [] plistcpu;

};

// Author: Kharchenko S.A.
// CSMatrixR: Get the submatrix from other processors described by the list of intervals
//========================================================================================
CSMatrixR CSMatrixR::GetSubmatrixViaIntervals (CMPIComm &_comm, // Get the matrix part from other processors described by the list of intervals
																int _nparts, int *_ibegarr, int *_iendarr) {

	const char *funcname = "GetSubmatrixViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	int i;

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::GetSubmatrixViaIntervals: intervals from different cpu's intersect ";
		};
	};

// For local list data find to which processor the data should be sent

	int *list2cpu;

	list2cpu = new int [nlist];
	if(!list2cpu) MemoryFail(funcname);

	for (i=0;i<nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ip, irow, iproc;

	if (nptot > 0) {

		ipbeg = 0;

		for (i=0;i<nlist;i++) {
			irow = list[i];
			ipend = nptot-1;
			ip = ipbeg;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[i] = icputot[ip];
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						list2cpu[i] = icputot[ip];
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;
	int *nzja2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);
	nzja2cpu = new int [nproc];
	if(!nzja2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;
	for (i=0;i<nproc;i++) nzja2cpu[i] = 0;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
			nzja2cpu[iproc] += (ia[i+1]-ia[i]);
		};
	};

	CSMatrixR *mtrarr;

	mtrarr = new CSMatrixR [nproc];
	if(!mtrarr) MemoryFail(funcname);

	for (i=0;i<nproc;i++) {
		CSMatrixR temp (nlist2cpu[i],nzja2cpu[i]);
		mtrarr[i] = temp;
	};

	int *nl2cpu;
	int *nz2cpu;

	nl2cpu = new int [nproc];
	if(!nl2cpu) MemoryFail(funcname);
	nz2cpu = new int [nproc];
	if(!nz2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nl2cpu[i] = 0;
	for (i=0;i<nproc;i++) nz2cpu[i] = 0;

	int *plist, *pia, *pja;
	double *pa;

	for (i=0;i<nproc;i++) {
		pia = mtrarr[i].GetIa ();
		pia[0] = 0;
	};

	int ilist, nzloc;

	for (i=0;i<nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			plist = mtrarr[iproc].GetList ();
			pia = mtrarr[iproc].GetIa ();
			pja = mtrarr[iproc].GetJa ();
			pa = mtrarr[iproc].GetA ();
			ilist = nl2cpu[iproc];
			nzloc = nz2cpu[iproc];
			plist[ilist] = list[i];
			for (j=ia[i];j<ia[i+1];j++) {
				pja[nzloc] = ja[j];
				pa[nzloc] = a[j];
				nzloc++;
			};
			pia[ilist+1] = nzloc;
			nl2cpu[iproc]++;
			nz2cpu[iproc] = nzloc;
		};
	};

	for (i=0;i<nproc;i++) {
		mtrarr[i].SetM (m);
		mtrarr[i].SetN (n);
		mtrarr[i].SetNsupr (nsupr);
		mtrarr[i].SetNsupcBase (nsupc);
		mtrarr[i].SetNlist (nlist2cpu[i]);
		mtrarr[i].SetNzja (nzja2cpu[i]);
		mtrarr[i].SetNza (nzja2cpu[i]);
		mtrarr[i].SetNzatot (nzja2cpu[i]);
	};

// Exchange the data

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
		mtrarr[i].PackMatrix (ObjSizeSend[i],ObjSend[i]);
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
	if(info) throw " CSMatrix::GetSubmatrixViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

	if (NObjRecv != nproc) throw " CSMatrix::GetSubmatrixViaIntervals: wrong number of the received objects ";

	for(i = 0; i < NObjRecv; i++) {
		mtrarr[i].UnPackMatrix (ObjSizeRecv[i],ObjRecv[i]);
	};

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Collect submatrices of the current cpu

	int *listblk;

	listblk = new int [nproc];
	if (!listblk) MemoryFail (funcname);

	for (i=0;i<nproc;i++) listblk[i] = i;

	CSMatrixR mtrnew;

	CGSMatrixRS::CombineMatricesSort (n, nproc, listblk,
												mtrarr,
												mtrnew);

	int nloc = mtrnew.GetN ();

	mtrnew.SetNsupr (nloc);
	mtrnew.SetNsupcBase (nloc);

	CSMatrixR mtrdummy;

	for (i=0;i<nproc;i++) mtrarr[i] = mtrdummy;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] list2cpu;
	delete [] nlist2cpu;
	delete [] nzja2cpu;
	delete [] mtrarr;
	delete [] nl2cpu;
	delete [] nz2cpu;
	delete [] listblk;

// Return the result

	return mtrnew;

};

// Author: Kharchenko S.A.
// CSMatrixR: Get the subvectors from other processors described by the list of intervals
//========================================================================================
void CSMatrixR::GetSubvectorsViaIntervals (CMPIComm &_comm, // Get the subvectors from other processors described by the list of intervals
															int _nparts, int *_ibegarr, int *_iendarr, 
															int _nrhs, int _nlist, int *_list, double *_vectini,
															int &_nlistnew, int *&_listnew, double *&_vectnew) {

	const char *funcname = "GetSubvectorsViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	int i;

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::GetSubvectorsViaIntervals: intervals from different cpu's intersect ";
		};
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	_nlistnew = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = _nlistnew;
			_nlistnew += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

	_listnew = new int [_nlistnew];
	if(!_listnew) MemoryFail(funcname);
	_vectnew = new double [_nlistnew*_nrhs];
	if(!_vectnew) MemoryFail(funcname);

	_nlistnew = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			for (j=ibegtot[i];j<=iendtot[i];j++) {
				_listnew[_nlistnew] = j;
				_nlistnew++;
			};
		};
	};

// Sort the list

	CIntInt *iivarr;

	iivarr = new CIntInt [_nlist];
	if(!iivarr) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) {
		iivarr[i].intvalue = _list[i];
		iivarr[i].int2value = i;
	};

	qsort (iivarr, _nlist, sizeof(CIntInt), CIntInt::CompareIntInt);

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [_nlist];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [_nlist];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ip, irow, iproc, ind;

	if (nptot > 0) {

		ipbeg = 0;

		for (i=0;i<_nlist;i++) {
			ind = iivarr[i].int2value;
			irow = _list[ind];
			ipend = nptot-1;
			ip = ipbeg;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[ind] = icputot[ip];
				interval2cpu[ind] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[ind] = icputot[ipbeg];
				interval2cpu[ind] = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsvect;
	int *isizevect;
	double **pdvect;

	ibsvect = new int [nproc];
	if(!ibsvect) MemoryFail(funcname);
	isizevect = new int [nproc];
	if(!isizevect) MemoryFail(funcname);
	pdvect = new double * [nproc];
	if(!pdvect) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsvect[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		isizevect[i] = niloc*(_nrhs+2);
		pdvect[i] = new double [niloc*(_nrhs+2)];
		if(!pdvect[i]) MemoryFail(funcname);
	};

	double *pv;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		ibs = ibsvect[iproc];
		pv = pdvect[iproc];
		pv[ibs] = (double) _list[i];
		ibs++;
		pv[ibs] = (double) interval2cpu[i];
		ibs++;
		for (j=0;j<_nrhs;j++) {
			pv[ibs] = _vectini[j*_nlist+i];
			ibs++;
		};
		ibsvect[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = sizeof(double) * isizevect[i];
		ObjSend[i] = (char *)pdvect[i];
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
	if(info) throw " CSMatrix::GetSubmatrixViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] pdvect[i];

	if (NObjRecv != nproc) throw " CSMatrix::GetSubvectorViaIntervals: wrong number of the received objects ";

// Store received data in the output array

	int ni, iv, k;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(double)*(_nrhs+2));
		ibs = 0;
		pv = (double *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = (int) pv[ibs];
			ibs++;
			iv = (int) pv[ibs];
			ibs++;
			ip = ibsinterval[iv]+ind-ibegtot[iv];
			for (k=0;k<_nrhs;k++) {
				_vectnew[k*_nlistnew+ip] = pv[ibs];
				ibs++;
			};
		};
	};

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] ibsinterval;
	delete [] iivarr;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsvect;
	delete [] isizevect;
	delete [] pdvect;

};

// Author: Kharchenko S.A.
// CSMatrixR: Get the subvectors from other processors described by the list of entry intervals and list of final indices
//========================================================================================
void CSMatrixR::GetSubvectorsViaEntryIntervals (CMPIComm &_comm, // Get the subvectors from other processors described by the list of entry intervals and list of final indices
																int _nparts, int *_ibegarr, int *_iendarr, 
																int _nrhs, double *_vectini,
																int _nlistnew, int *_listnew, double *_vectnew) {

	const char *funcname = "GetSubvectorsViaEntryIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that on entry all intervals are sorted

	int i;

	for (i=0;i<_nparts-1;i++) {
		if (_ibegarr[i]<=_iendarr[i] && _ibegarr[i+1] <= _iendarr[i+1]) {
			if (_ibegarr[i] >= _iendarr[i+1]) throw " CSMatrix::GetSubvectorsViaEntryIntervals: intervals on entry are not sorted locally ";
		};
	};

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::GetSubvectorsViaEntryIntervals: intervals from different cpu's intersect ";
		};
	};

// Check that intervals cover the whole line

	if (nptot > 0) {
		if (ibegtot[0] != 0) throw " CSMatrix::GetSubvectorsViaEntryIntervals: intervals start is incorrect ";
		for (i=0;i<nptot-1;i++) {
			if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
				if (iendtot[i]+1 != ibegtot[i+1]) throw " CSMatrix::GetSubvectorsViaEntryIntervals: there are wholes in the intervals ";
			};
		};
	};

// Create the partitioning of the initial vector data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	int nlistini = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = nlistini;
			nlistini += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

// Sort the final list

	CIntInt *iivarr;

	iivarr = new CIntInt [_nlistnew];
	if(!iivarr) MemoryFail(funcname);

	for (i=0;i<_nlistnew;i++) {
		iivarr[i].intvalue = _listnew[i];
		iivarr[i].int2value = i;
	};

	qsort (iivarr, _nlistnew, sizeof(CIntInt), CIntInt::CompareIntInt);

// For local list data find from which processor the data should be received

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [_nlistnew];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [_nlistnew];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<_nlistnew;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ip, irow, iproc, ind;

	if (nptot > 0) {

		ipbeg = 0;

		for (i=0;i<_nlistnew;i++) {
			ind = iivarr[i].int2value;
			irow = _listnew[ind];
			ipend = nptot-1;
			ip = ipbeg;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[ind] = icputot[ip];
				interval2cpu[ind] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[ind] = icputot[ipbeg];
				interval2cpu[ind] = ipbeg;
			};
		};

	};

// Prepare the data for initial send

	int *nlist2cpu;
	int **plists;
	int *iptr;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);
	iptr = new int [nproc];
	if(!iptr) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<_nlistnew;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	for (i=0;i<nproc;i++) {
		plists[i] = new int [nlist2cpu[i]*2];
		if(!plists[i]) MemoryFail(funcname);
	};

	for (i=0;i<nproc;i++) iptr[i] = 0;

	int k;

	for (i=0;i<_nlistnew;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			k = iptr[iproc];
			plists[iproc][k*2] = _listnew[i];
			plists[iproc][k*2+1] = interval2cpu[i];
			iptr[iproc]++;
		};
	};

// Exchange the data

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
		ObjSizeSend[i] = 2 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrix::GetSubmatrixViaEntryIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrix::GetSubvectorViaEntryIntervals: wrong number of the received objects ";

// Prepare the data for new send

	int *isizevect;
	double **pdvect;

	isizevect = new int [nproc];
	if(!isizevect) MemoryFail(funcname);
	pdvect = new double * [nproc];
	if(!pdvect) MemoryFail(funcname);

	int niloc;

	int *piarr;
	double *pv;

	for (i=0;i<NObjRecv;i++) {
		iproc = CpuIDRecv[i];
		niloc = ObjSizeRecv[i] / (sizeof(int)*2);
		isizevect[iproc] = niloc*_nrhs;
		pdvect[iproc] = new double [niloc*_nrhs];
		if(!pdvect[iproc]) MemoryFail(funcname);
		piarr = (int *)ObjRecv[i];
		pv = pdvect[iproc];
		for (k=0;k<niloc;k++) {
			ind = piarr[k*2];
			ip = piarr[k*2+1];
			ibs = ibsinterval[ip]+ind-ibegtot[ip];
			for (j=0;j<_nrhs;j++) {
				pv[k*_nrhs+j] = _vectini[j*nlistini+ibs];
			};
		};
	};

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Exchange the data

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
		ObjSizeSend[i] = sizeof(double) * isizevect[i];
		ObjSend[i] = (char *)pdvect[i];
	};

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv);
	if(info) throw " CSMatrix::GetSubmatrixViaEntryIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] pdvect[i];

	if (NObjRecv != nproc) throw " CSMatrix::GetSubvectorViaEntryIntervals: wrong number of the received objects ";

// Store received data in the output array

	int *cpu2recv;

	cpu2recv = new int [nproc];
	if(!cpu2recv) MemoryFail(funcname);

	for (i=0;i<NObjRecv;i++) {
		iproc = CpuIDRecv[i];
		cpu2recv[iproc] = i;
	};

	for (i=0;i<nproc;i++) iptr[i] = 0;

	int irecv;

	for (i=0;i<_nlistnew;i++) {
		iproc = list2cpu[i];
		irecv = cpu2recv[iproc];
		pv = (double *)ObjRecv[irecv];
		k = iptr[iproc];
		for (j=0;j<_nrhs;j++) {
			_vectnew[j*_nlistnew+i] = pv[k*_nrhs+j];
		};
		iptr[iproc]++;
	};

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] ibsinterval;
	delete [] iivarr;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] plists;
	delete [] iptr;
	delete [] isizevect;
	delete [] pdvect;
	delete [] cpu2recv;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute new reordered list for the ordering presented for intervals
//========================================================================================
void CSMatrix::OrderViaIntervals (CMPIComm &_comm, // Compute new reordered list for the ordering presented for intervals
												int _nlist, int *_listini,
												int _nparts, int *_ibegarr, int *_iendarr, 
												int *_order,
												int *_listnew) {

	const char *funcname = "OrderViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that on entry all intervals are sorted

	int i;

	for (i=0;i<_nparts-1;i++) {
		if (_ibegarr[i]<=_iendarr[i] && _ibegarr[i+1] <= _iendarr[i+1]) {
			if (_ibegarr[i] >= _iendarr[i+1]) throw " CSMatrix::OrderViaIntervals: intervals on entry are not sorted locally ";
		};
	};

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::OrderViaIntervals: intervals from different cpu's intersect ";
		};
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	int nlistnew = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = nlistnew;
			nlistnew += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

// Sort the list

	CIntInt *iivarr;

	iivarr = new CIntInt [_nlist];
	if(!iivarr) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) {
		iivarr[i].intvalue = _listini[i];
		iivarr[i].int2value = i;
	};

	qsort (iivarr, _nlist, sizeof(CIntInt), CIntInt::CompareIntInt);

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [_nlist];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [_nlist];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ip, irow, iproc, ind;

	if (nptot > 0) {

		ipbeg = 0;

		for (i=0;i<_nlist;i++) {
			ind = iivarr[i].int2value;
			irow = _listini[ind];
			ipend = nptot-1;
			ip = ipbeg;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[ind] = icputot[ip];
				interval2cpu[ind] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[ind] = icputot[ipbeg];
				interval2cpu[ind] = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsarr;
	int **plists;

	ibsarr = new int [nproc];
	if(!ibsarr) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsarr[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		plists[i] = new int [niloc*2];
		if(!plists[i]) MemoryFail(funcname);
	};

	int *piarr;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		ibs = ibsarr[iproc];
		piarr = plists[iproc];
		piarr[ibs] = _listini[i];
		ibs++;
		piarr[ibs] = interval2cpu[i];
		ibs++;
		ibsarr[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = 2 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrix::OrderViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrix::OrderViaIntervals: wrong number of the received objects ";

// Recompute the result

	int ni, iv, indnew;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*2);
		piarr = (int *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*2];
			iv = piarr[j*2+1];
			ip = ibsinterval[iv]+ind-ibegtot[iv];
			indnew = _order[ip];
			piarr[j*2+1] = indnew;
		};
	};

// Exchange the data back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::OrderViaIntervals: Error in DataExchangeMPI";

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Store received data in the output array

	int *cpu2recv;
	int *iptr;

	cpu2recv = new int [nproc];
	if(!cpu2recv) MemoryFail(funcname);
	iptr = new int [nproc];
	if(!iptr) MemoryFail(funcname);

	for (i=0;i<NObjSend;i++) {
		iproc = CpuIDSend[i];
		cpu2recv[iproc] = i;
	};

	for (i=0;i<nproc;i++) iptr[i] = 0;

	int irecv, k;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		irecv = cpu2recv[iproc];
		piarr = (int *)ObjSend[irecv];
		k = iptr[iproc];
		ind = piarr[k*2];
		indnew = piarr[k*2+1];
		if (ind != _listini[i]) {
			throw " CSMatrix::OrderViaIntervals: wrong obtained index ";
		};
		_listnew[i] = indnew;
		iptr[iproc]++;
	};

// Free send data

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for (i=0;i<NObjSend;i++) delete [] ObjSend[i];
	delete [] ObjSend;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] ibsinterval;
	delete [] iivarr;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsarr;
	delete [] plists;
	delete [] cpu2recv;
	delete [] iptr;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute direct order for the list described by the intervals
//========================================================================================
void CSMatrix::DirOrderViaIntervals (CMPIComm &_comm, // Compute inverse order for the list described by the intervals
													int _nlist, int *_listini, int *_listord,
													int _nparts, int *_ibegarr, int *_iendarr, 
													int &_nlistpart, int *&_orderpart) {

	const char *funcname = "DirOrderViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that on entry all intervals are sorted

	int i;

	for (i=0;i<_nparts-1;i++) {
		if (_ibegarr[i]<=_iendarr[i] && _ibegarr[i+1] <= _iendarr[i+1]) {
			if (_ibegarr[i] >= _iendarr[i+1]) throw " CSMatrix::InvOrderViaIntervals: intervals on entry are not sorted locally ";
		};
	};

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::InvOrderViaIntervals: intervals from different cpu's intersect ";
		};
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	_nlistpart = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = _nlistpart;
			_nlistpart += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

	_orderpart = new int [_nlistpart];
	if(!_orderpart) MemoryFail(funcname);

// Sort the list

	CIntInt *iivarr;

	iivarr = new CIntInt [_nlist];
	if(!iivarr) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) {
		iivarr[i].intvalue = _listini[i];
		iivarr[i].int2value = i;
	};

	qsort (iivarr, _nlist, sizeof(CIntInt), CIntInt::CompareIntInt);

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [_nlist];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [_nlist];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ipprev, ip, irow, iproc, ind;

	if (nptot > 0) {

		ipprev = 0;

		for (i=0;i<_nlist;i++) {
			ind = iivarr[i].int2value;
			irow = _listini[ind];
			ipbeg = 0;
			ipend = nptot-1;
			ip = ipprev;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[ind] = icputot[ip];
				interval2cpu[ind] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[ind] = icputot[ipbeg];
				interval2cpu[ind] = ipbeg;
				ipprev = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsarr;
	int **plists;

	ibsarr = new int [nproc];
	if(!ibsarr) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsarr[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		plists[i] = new int [niloc*3];
		if(!plists[i]) MemoryFail(funcname);
	};

	int *piarr;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		ibs = ibsarr[iproc];
		piarr = plists[iproc];
		piarr[ibs] = _listini[i];
		ibs++;
		piarr[ibs] = _listord[i];
		ibs++;
		piarr[ibs] = interval2cpu[i];
		ibs++;
		ibsarr[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = 3 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrix::InvOrderViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrix::InvOrderViaIntervals: wrong number of the received objects ";

// Store the result

	int ni, iv, indnew;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*3);
		piarr = (int *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*3];
			indnew = piarr[j*3+1];
			iv = piarr[j*3+2];
			ip = ibsinterval[iv]+ind-ibegtot[iv];
			_orderpart[ip] = indnew;
		};
	};

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] ibsinterval;
	delete [] iivarr;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsarr;
	delete [] plists;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute inverse order for the list described by the intervals
//========================================================================================
void CSMatrix::InvOrderViaIntervals (CMPIComm &_comm, // Compute inverse order for the list described by the intervals
													int _nlist, int *_listini, int *_listord,
													int _nparts, int *_ibegarr, int *_iendarr, 
													int &_nlistpart, int *&_invorderpart) {

	const char *funcname = "InvOrderViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that on entry all intervals are sorted

	int i;

	for (i=0;i<_nparts-1;i++) {
		if (_ibegarr[i]<=_iendarr[i] && _ibegarr[i+1] <= _iendarr[i+1]) {
			if (_ibegarr[i] >= _iendarr[i+1]) throw " CSMatrix::InvOrderViaIntervals: intervals on entry are not sorted locally ";
		};
	};

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	int nptot = blksnp[nproc];

	int *ibegtot, *iendtot, *icputot;

	ibegtot = new int [nptot];
	if(!ibegtot) MemoryFail(funcname);
	iendtot = new int [nptot];
	if(!iendtot) MemoryFail(funcname);
	icputot = new int [nptot];
	if(!icputot) MemoryFail(funcname);

	for (i=0;i<nptot;i++) ibegtot[i] = 0;
	for (i=0;i<nptot;i++) iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, ibegtot, ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nptot, iendtot, iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		iiarr[i].intvalue = ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = ibegtot[j];
	};
	for (i=0;i<nptot;i++) ibegtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = iendtot[j];
	};
	for (i=0;i<nptot;i++) iendtot[i] = iwork[i];

	for (i=0;i<nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = icputot[j];
	};
	for (i=0;i<nptot;i++) icputot[i] = iwork[i];

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (ibegtot[i]<=iendtot[i] && ibegtot[i+1] <= iendtot[i+1]) {
			if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::InvOrderViaIntervals: intervals from different cpu's intersect ";
		};
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	_nlistpart = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = _nlistpart;
			_nlistpart += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

	_invorderpart = new int [_nlistpart];
	if(!_invorderpart) MemoryFail(funcname);

// Sort the list

	CIntInt *iivarr;

	iivarr = new CIntInt [_nlist];
	if(!iivarr) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) {
		iivarr[i].intvalue = _listord[i];
		iivarr[i].int2value = i;
	};

	qsort (iivarr, _nlist, sizeof(CIntInt), CIntInt::CompareIntInt);

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [_nlist];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [_nlist];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<_nlist;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ipprev, ip, irow, iproc, ind;

	if (nptot > 0) {

		ipprev = 0;

		for (i=0;i<_nlist;i++) {
			ind = iivarr[i].int2value;
			irow = _listord[ind];
			ipbeg = 0;
			ipend = nptot-1;
			ip = ipprev;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[ind] = icputot[ip];
				interval2cpu[ind] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[ind] = icputot[ipbeg];
				interval2cpu[ind] = ipbeg;
				ipprev = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsarr;
	int **plists;

	ibsarr = new int [nproc];
	if(!ibsarr) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsarr[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		plists[i] = new int [niloc*3];
		if(!plists[i]) MemoryFail(funcname);
	};

	int *piarr;

	for (i=0;i<_nlist;i++) {
		iproc = list2cpu[i];
		ibs = ibsarr[iproc];
		piarr = plists[iproc];
		piarr[ibs] = _listini[i];
		ibs++;
		piarr[ibs] = _listord[i];
		ibs++;
		piarr[ibs] = interval2cpu[i];
		ibs++;
		ibsarr[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = 3 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrix::InvOrderViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrix::InvOrderViaIntervals: wrong number of the received objects ";

// Store the result

	int ni, iv, indnew;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*3);
		piarr = (int *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*3];
			indnew = piarr[j*3+1];
			iv = piarr[j*3+2];
			ip = ibsinterval[iv]+indnew-ibegtot[iv];
			_invorderpart[ip] = ind;
		};
	};

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Free work arrays

	delete [] blksnp;
	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] iiarr;
	delete [] iwork;
	delete [] ibsinterval;
	delete [] iivarr;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsarr;
	delete [] plists;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute new reordered local vector data according to the new global ordering of the local data
//========================================================================================
void CSMatrix::ChangeVectorIndices (int _nrhs, // Compute new reordered local vector data according to the new global ordering of the local data
												int _nlist, int *_list, double *_vect, int *_order) {

	const char *funcname = "ChangeVectorIndices";

// Compute new local ordering

	CIntInt *iiarr;

	iiarr = new CIntInt [_nlist];
	if(!iiarr) MemoryFail(funcname);

	int i;

	for (i=0;i<_nlist;i++) {
		iiarr[i].intvalue = _order[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, _nlist, sizeof(CIntInt), CIntInt::CompareIntInt);

// Reorder

	double *vectloc;

	vectloc = new double [_nlist];
	if(!vectloc) MemoryFail(funcname);

	int irhs, inew;

	for (i=0;i<_nlist;i++) {
		inew = iiarr[i].int2value;
		_list[inew] = _order[i];
	};

	for (irhs=0;irhs<_nrhs;irhs++) {
		for (i=0;i<_nlist;i++) {
			inew = iiarr[i].int2value;
			vectloc[inew] = _vect[irhs*_nlist+i];
		};
		for (i=0;i<_nlist;i++) _vect[irhs*_nlist+i] = vectloc[i];
	};

// Free work arrays

	delete [] iiarr;
	delete [] vectloc;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute symmetrized matrix
//========================================================================================
CSMatrix CSMatrix::SymmMtr (CMPIComm &_comm, // Compute symmetrized matrix
										int _nblks, int *_blks, int *_blk2cpu) const { 

	const char *funcname = "SymmMtr";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// For each column element find its block number

	int nzjaloc = nzja;

	int *ja2blk;

	ja2blk = new int [nzjaloc];
	if (!ja2blk) MemoryFail (funcname);

	int i;

	for (i=0;i<nzjaloc;i++) ja2blk[i] = -1;

	int iblkprev = 0;

	int jj, iblkbeg, iblkend, iblk;

	for (i=0;i<nzjaloc;i++) {
		jj = ja[i];
		if (jj >= _blks[iblkprev] && jj < _blks[iblkprev+1]) {
			ja2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = _nblks-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend) / 2;
				if (jj >= _blks[iblk] && jj < _blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < _blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= _blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			ja2blk[i] = iblkbeg;
		};
		iblkprev = ja2blk[i];
	};

// Compute local list base addresses

	int *list2blk;
	int *ibslist;

	list2blk = new int [nlist];
	if (!list2blk) MemoryFail (funcname);
	ibslist = new int [_nblks];
	if (!ibslist) MemoryFail (funcname);

	int nz = 0;

	int j;

	for (i=0;i<_nblks;i++) {
		if (_blk2cpu[i] == myid) {
			ibslist[i] = nz;
			for (j=_blks[i];j<_blks[i+1];j++) {
				list2blk[nz] = i;
				nz++;
			};
		} else {
			ibslist[i] = -1;
		};
	};

// Create the lists of indices for each processor

	int *nzcpu;

	nzcpu = new int [nproc];
	if (!nzcpu) MemoryFail (funcname);

	for (i=0;i<nproc;i++) nzcpu[i] = 0;

	int jblk, jjloc, jcpu, i1, icheck, k, kk;

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jblk = ja2blk[j];
			jjloc = jj-_blks[jblk];
			jcpu = _blk2cpu[jblk];
			if (jcpu == myid) {
				i1 = ibslist[jblk]+jjloc;
				icheck = -1;
				for (k=ia[i1];k<ia[i1+1];k++) {
					kk = ja[k];
					if (kk == list[i]) icheck = 0;
				};
				if (icheck == -1) nzcpu[myid]++;
			} else {
				nzcpu[jcpu]++;
			};
		};
	};

	int **plists;

	plists = new int * [nproc];
	if (!plists) MemoryFail (funcname);

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nzcpu[i];
		plists[i] = new int [niloc*4];
		if (!plists[i]) MemoryFail (funcname);
		nzcpu[i] = 0;
	};

	int irow, ip;

	for (i=0;i<nlist;i++) {
		irow = list[i];
		iblk = list2blk[i];
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			jblk = ja2blk[j];
			jjloc = jj-_blks[jblk];
			jcpu = _blk2cpu[jblk];
			if (jcpu == myid) {
				i1 = ibslist[jblk]+jjloc;
				icheck = -1;
				for (k=ia[i1];k<ia[i1+1];k++) {
					kk = ja[k];
					if (kk == list[i]) icheck = 0;
				};
				if (icheck == -1) {
					ip = nzcpu[myid];
					plists[myid][ip*4] = irow;
					plists[myid][ip*4+1] = iblk;
					plists[myid][ip*4+2] = jj;
					plists[myid][ip*4+3] = jblk;
					nzcpu[myid]++;
				};
			} else {
				ip = nzcpu[jcpu];
				plists[jcpu][ip*4] = irow;
				plists[jcpu][ip*4+1] = iblk;
				plists[jcpu][ip*4+2] = jj;
				plists[jcpu][ip*4+3] = jblk;
				nzcpu[jcpu]++;
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

	for (i = 0; i < NObjSend; i++) {
		ObjTypeSend[i] = 1;
		ObjIDSend[i] = i;
		CpuIDSend[i] = i;
		ObjSizeSend[i] = sizeof(int) * nzcpu[i] * 4;
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
	if (info) throw " CSMatrix::SymmMtr: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

// Add received data into the sparsity structure

	int *ialoc;
	int *iptr;

	ialoc = new int [nlist+1];
	if(!ialoc) MemoryFail(funcname);
	iptr = new int [nlist];
	if(!iptr) MemoryFail(funcname);

	ialoc[0] = 0;
	for (i=0;i<nlist;i++) ialoc[i+1] = ia[i+1]-ia[i];

	if (NObjRecv != nproc) throw " CSMatrix::SymmMtr: wrong number of the received objects ";

	int ni;
	int *piarr;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*4);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			irow = piarr[j*4];
			iblk = piarr[j*4+1];
			jj   = piarr[j*4+2];
			jblk = piarr[j*4+3];
			jjloc = jj-_blks[jblk];
			i1 = ibslist[jblk]+jjloc;
			icheck = -1;
			for (k=ia[i1];k<ia[i1+1];k++) {
				kk = ja[k];
				if (kk == irow) icheck = 0;
			};
			if (icheck == 0) {
				piarr[j*4] = -1;
				piarr[j*4+1] = -1;
				piarr[j*4+2] = -1;
				piarr[j*4+3] = -1;
			} else {
				ialoc[i1+1]++;
			};
		};
	};

	for (i=0;i<nlist;i++) ialoc[i+1] = ialoc[i]+ialoc[i+1];

	int nzjanew = ialoc[nlist];

	CSMatrix mtrnew (nlist,nzjanew);

	int *plist = mtrnew.list;
	int *pia = mtrnew.ia;
	int *pja = mtrnew.ja;

	for (i=0;i<nlist;i++) plist[i] = list[i];
	for (i=0;i<=nlist;i++) pia[i] = ialoc[i];

	for (i=0;i<nlist;i++) iptr[i] = ialoc[i];

	for (i=0;i<nlist;i++) {
		for (j=ia[i];j<ia[i+1];j++) {
			jj = ja[j];
			k = iptr[i];
			pja[k] = jj;
			iptr[i]++;
		};
	};

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*4);
		piarr = (int *) ObjRecv[i];
		for (j=0;j<ni;j++) {
			irow = piarr[j*4];
			iblk = piarr[j*4+1];
			jj   = piarr[j*4+2];
			jblk = piarr[j*4+3];
			if (irow >= 0) {
				jjloc = jj-_blks[jblk];
				i1 = ibslist[jblk]+jjloc;
				k = iptr[i1];
				pja[k] = irow;
				iptr[i1]++;
			};
		};
	};

// Free received data

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

// Perform sorting of the data

	for (i=0;i<nlist;i++) {
		ni = pia[i+1]-pia[i];
		qsort (pja+pia[i],ni,sizeof(int),compint);
	};

	mtrnew.n = n;
	mtrnew.m = m;
	mtrnew.nsupc = n;
	mtrnew.nsupr = m;
	mtrnew.nlist = nlist;
	mtrnew.nzja = nzjanew;

// Free work arrays

	delete [] ja2blk;
	delete [] list2blk;
	delete [] ibslist;
	delete [] nzcpu;
	delete [] plists;
	delete [] ialoc;
	delete [] iptr;

	return mtrnew;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute global intervals
//========================================================================================
void CSMatrix::ComputeIntervals (CMPIComm &_comm, // Compute global intervals
												int _nparts, int *_ibegarr, int *_iendarr, 
												int &_nptot, int *&_ibegtot, int *&_iendtot, int *&_icputot) {

	const char *funcname = "ComputeIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Check that on entry all intervals are sorted

	int i;

	for (i=0;i<_nparts-1;i++) {
		if (_ibegarr[i] >= _iendarr[i+1]) throw " CSMatrix::ComputeIntervals: intervals on entry are not sorted locally ";
	};

// Create global intervals description

	int *blksnp;

	blksnp = new int [nproc+1];
	if(!blksnp) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) blksnp[i] = 0;

	blksnp[myid+1] = _nparts;

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, nproc+1, blksnp, blksnp);

	for (i=0;i<nproc;i++) blksnp[i+1] = blksnp[i]+blksnp[i+1];

	_nptot = blksnp[nproc];

	_ibegtot = new int [_nptot];
	if(!_ibegtot) MemoryFail(funcname);
	_iendtot = new int [_nptot];
	if(!_iendtot) MemoryFail(funcname);
	_icputot = new int [_nptot];
	if(!_icputot) MemoryFail(funcname);

	for (i=0;i<_nptot;i++) _ibegtot[i] = 0;
	for (i=0;i<_nptot;i++) _iendtot[i] = 0;

	int ibs = blksnp[myid];

	for (i=0;i<_nparts;i++) _ibegtot[ibs+i] = _ibegarr[i];
	for (i=0;i<_nparts;i++) _iendtot[ibs+i] = _iendarr[i];

	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nptot, _ibegtot, _ibegtot);
	CMPIExchange::ExchangeArrayMPI (_comm, INTEGERVALUE, ADD, _nptot, _iendtot, _iendtot);

	int j;

	for (i=0;i<nproc;i++) {
		for (j=blksnp[i];j<blksnp[i+1];j++) {
			_icputot[j] = i;
		};
	};

	CIntInt *iiarr;

	iiarr = new CIntInt [_nptot];
	if(!iiarr) MemoryFail(funcname);

	for (i=0;i<_nptot;i++) {
		iiarr[i].intvalue = _ibegtot[i];
		iiarr[i].int2value = i;
	};

	qsort (iiarr, _nptot, sizeof(CIntInt),CIntInt::CompareIntInt);

	int *iwork;

	iwork = new int [_nptot];
	if(!iwork) MemoryFail(funcname);

	for (i=0;i<_nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = _ibegtot[j];
	};
	for (i=0;i<_nptot;i++) _ibegtot[i] = iwork[i];

	for (i=0;i<_nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = _iendtot[j];
	};
	for (i=0;i<_nptot;i++) _iendtot[i] = iwork[i];

	for (i=0;i<_nptot;i++) {
		j = iiarr[i].int2value;
		iwork[i] = _icputot[j];
	};
	for (i=0;i<_nptot;i++) _icputot[i] = iwork[i];

// Free work arrays

	delete [] blksnp;
	delete [] iiarr;
	delete [] iwork;

};

// Author: Kharchenko S.A.
// CSMatrix: Compute inverse order for the list described by the intervals
//========================================================================================
void CSMatrix::OrderViaIntervalsAndInverseOrder (CMPIComm &_comm, // Compute new reordered list for the inverse ordering presented for intervals
																	int _nlistini, int *_listini,
																	int _nlisttree, int *_listtree, int *_invordertree, 
																	int _nparts, int *_ibegarr, int *_iendarr, int *_orderpart,
																	int *_listord) {

	const char *funcname = "OrderViaIntervalsAndInverseOrder";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Partition the total indices space

	int ntot = 0;

	int i, jj;

	for (i=0;i<_nlistini;i++) {
		jj = _listini[i];
		if (jj > ntot) ntot = jj;
	};

	CMPIExchange::ExchangeArrayMPI (_comm,
												INTEGERVALUE, MAXIMUM, 1, &ntot, &ntot);

	ntot++;

	int ni = ntot / nproc;

	int *blks;

	blks = new int [nproc+1];
	if(!blks) MemoryFail(funcname);

	for (i=0;i<nproc;i++) blks[i] = i*ni;

	blks[nproc] = ntot;

// Create the mask

	int nimax = 0;

	for (i=0;i<nproc;i++) {
		ni = blks[i+1]-blks[i];
		if (ni > nimax) nimax = ni;
	};

	int *imaskloc;
	int *listloc;

	imaskloc = new int [nimax];
	if(!imaskloc) MemoryFail(funcname);
	listloc = new int [nimax];
	if(!listloc) MemoryFail(funcname);

	int icycle = -1;

	for (i=0;i<nimax;i++) imaskloc[i] = icycle;

// Split list data into the blocks

	int *list2blk;

	list2blk = new int [_nlistini];
	if(!list2blk) MemoryFail(funcname);

	int iblkprev = 0;

	int iblk, iblkbeg, iblkend;

	for (i=0;i<_nlistini;i++) {
		jj = _listini[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			list2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			list2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iablk;
	int *iptr;
	int *jablk2ind;

	iablk = new int [nproc+1];
	if(!iablk) MemoryFail(funcname);
	iptr = new int [nproc];
	if(!iptr) MemoryFail(funcname);
	jablk2ind = new int [_nlistini];
	if(!jablk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iablk[i] = 0;

	for (i=0;i<_nlistini;i++) {
		iblk = list2blk[i];
		iablk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iablk[i+1] = iablk[i]+iablk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iablk[i];

	int k;

	for (i=0;i<_nlistini;i++) {
		iblk = list2blk[i];
		k = iptr[iblk];
		jablk2ind[k] = i;
		iptr[iblk]++;
	};

// Split inverse ordering data into blocks

	int *inv2blk;

	inv2blk = new int [_nlisttree];
	if(!inv2blk) MemoryFail(funcname);

	iblkprev = 0;

	for (i=0;i<_nlisttree;i++) {
		jj = _invordertree[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			inv2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			inv2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iainvblk;
	int *jainvblk2ind;

	iainvblk = new int [nproc+1];
	if(!iainvblk) MemoryFail(funcname);
	jainvblk2ind = new int [_nlisttree];
	if(!jablk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iainvblk[i] = 0;

	for (i=0;i<_nlisttree;i++) {
		iblk = inv2blk[i];
		iainvblk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iainvblk[i+1] = iainvblk[i]+iainvblk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iainvblk[i];

	for (i=0;i<_nlisttree;i++) {
		iblk = inv2blk[i];
		k = iptr[iblk];
		jainvblk2ind[k] = i;
		iptr[iblk]++;
	};

// Perform initial filling of the ordering array

	for (i=0;i<_nlistini;i++) _listord[i] = -1;

	int j, ind, irow, irowold, ibs;

	for (iblk=0;iblk<nproc;iblk++) {
		icycle++;
		for (j=iainvblk[iblk];j<iainvblk[iblk+1];j++) {
			ind = jainvblk2ind[j];
			irow = _listtree[ind];
			irowold = _invordertree[ind];
			ibs = irowold-blks[iblk];
			imaskloc[ibs] = icycle;
			listloc[ibs] = irow;
		};
		for (j=iablk[iblk];j<iablk[iblk+1];j++) {
			ind = jablk2ind[j];
			irowold = _listini[ind];
			ibs = irowold-blks[iblk];
			if (imaskloc[ibs] == icycle) {
				_listord[ind] = listloc[ibs];
			};
		};
	};

// Free work arrays

	delete [] blks;
	delete [] listloc;
	delete [] imaskloc;
	delete [] list2blk;
	delete [] iablk;
	delete [] iptr;
	delete [] jablk2ind;
	delete [] inv2blk;
	delete [] iainvblk;
	delete [] jainvblk2ind;

// Create the list of remaining indices

	int nlistfin = 0;

	for (i=0;i<_nlistini;i++) {
		if (_listord[i] == -1) nlistfin++;
	};

	int *listfin;
	int *fin2ind;

	listfin = new int [nlistfin];
	if(!listfin) MemoryFail(funcname);
	fin2ind = new int [nlistfin];
	if(!fin2ind) MemoryFail(funcname);

	nlistfin = 0;

	for (i=0;i<_nlistini;i++) {
		if (_listord[i] == -1) {
			listfin[nlistfin] = _listini[i];
			fin2ind[nlistfin] = i;
			nlistfin++;
		};
	};

// Create global partitionings

	int nptot, *ibegtot, *iendtot, *icputot;

	CSMatrix::ComputeIntervals (_comm,
											_nparts, _ibegarr, _iendarr,
											nptot, ibegtot, iendtot, icputot);

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::OrderViaIntervalsAndInverseOrder: intervals from different cpu's intersect ";
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	int nlistpart = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = nlistpart;
			nlistpart += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [nlistfin];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [nlistfin];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<nlistfin;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ipprev, ip, iproc;

	if (nptot > 0) {

		ipprev = 0;

		for (i=0;i<nlistfin;i++) {
			irow = listfin[i];
			ipbeg = 0;
			ipend = nptot-1;
			ip = ipprev;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[i] = icputot[ip];
				interval2cpu[i] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[i] = icputot[ipbeg];
				interval2cpu[i] = ipbeg;
				ipprev = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<nlistfin;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsarr;
	int **plists;

	ibsarr = new int [nproc];
	if(!ibsarr) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsarr[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		plists[i] = new int [niloc*3];
		if(!plists[i]) MemoryFail(funcname);
	};

	int *piarr;

	for (i=0;i<nlistfin;i++) {
		iproc = list2cpu[i];
		ibs = ibsarr[iproc];
		piarr = plists[iproc];
		piarr[ibs] = listfin[i];
		ibs++;
		piarr[ibs] = fin2ind[i];
		ibs++;
		piarr[ibs] = interval2cpu[i];
		ibs++;
		ibsarr[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = 3 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrix::OrderViaIntervalsAndInverseOrder: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrix::OrderViaIntervalsAndInverseOrder: wrong number of the received objects ";

// Store the result

	int iv, indnew;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*3);
		piarr = (int *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*3];
			iv = piarr[j*3+2];
			ip = ibsinterval[iv]+ind-ibegtot[iv];
			indnew = _orderpart[ip];
			piarr[j*3+2] = indnew;
		};
	};

// Send back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::OrderViaIntervalsAndInverseOrder: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

	if (NObjSend != nproc) throw " CSMatrix::OrderViaIntervalsAndInverseOrder: wrong number of the received objects ";

// Use received data

	int irownew;

	for (i=0;i<NObjSend;i++) {
		ni = ObjSizeSend[i] / (sizeof(int)*3);
		piarr = (int *)ObjSend[i];
		for (j=0;j<ni;j++) {
			irow = piarr[j*3];
			ind = piarr[j*3+1];
			irownew = piarr[j*3+2];
			_listord[ind] = irownew;
		};
	};

// Free received data

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Free work arrays

	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] listfin;
	delete [] fin2ind;
	delete [] ibsinterval;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsarr;
	delete [] plists;

};

// Author: Kharchenko S.A.
// CSMatrix: Change rows/column indices
//========================================================================================
void CSMatrix::ChangeIndicesViaIntervals (CMPIComm &_comm, // Change rows/column indices
														int *_listord,
														int _nlisttree, int *_listtree, int *_invordertree, 
														int _nparts, int *_ibegarr, int *_iendarr, int *_orderpart) {

	const char *funcname = "ChangeIndicesViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Partition the total indices space

	int ntot = 0;

	int i, jj;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj > ntot) ntot = jj;
	};

	CMPIExchange::ExchangeArrayMPI (_comm,
												INTEGERVALUE, MAXIMUM, 1, &ntot, &ntot);

	ntot++;

	int ni = ntot / nproc;

	int *blks;

	blks = new int [nproc+1];
	if(!blks) MemoryFail(funcname);

	for (i=0;i<nproc;i++) blks[i] = i*ni;

	blks[nproc] = ntot;

// Create the mask

	int nimax = 0;

	for (i=0;i<nproc;i++) {
		ni = blks[i+1]-blks[i];
		if (ni > nimax) nimax = ni;
	};

	int *imaskloc;
	int *listloc;

	imaskloc = new int [nimax];
	if(!imaskloc) MemoryFail(funcname);
	listloc = new int [nimax];
	if(!listloc) MemoryFail(funcname);

	int icycle = -1;

	for (i=0;i<nimax;i++) imaskloc[i] = icycle;

// Split list data into the blocks

	int *list2blk;

	list2blk = new int [nlist];
	if(!list2blk) MemoryFail(funcname);

	int iblkprev = 0;

	int iblk, iblkbeg, iblkend;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			list2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			list2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iablk;
	int *iptr;
	int *jablk2ind;

	iablk = new int [nproc+1];
	if(!iablk) MemoryFail(funcname);
	iptr = new int [nproc];
	if(!iptr) MemoryFail(funcname);
	jablk2ind = new int [nlist];
	if(!jablk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iablk[i] = 0;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		iablk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iablk[i+1] = iablk[i]+iablk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iablk[i];

	int k;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		k = iptr[iblk];
		jablk2ind[k] = i;
		iptr[iblk]++;
	};

// Split inverse ordering data into blocks

	int *inv2blk;

	inv2blk = new int [_nlisttree];
	if(!inv2blk) MemoryFail(funcname);

	iblkprev = 0;

	for (i=0;i<_nlisttree;i++) {
		jj = _invordertree[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			inv2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			inv2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iainvblk;
	int *jainvblk2ind;

	iainvblk = new int [nproc+1];
	if(!iainvblk) MemoryFail(funcname);
	jainvblk2ind = new int [_nlisttree];
	if(!jablk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iainvblk[i] = 0;

	for (i=0;i<_nlisttree;i++) {
		iblk = inv2blk[i];
		iainvblk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iainvblk[i+1] = iainvblk[i]+iainvblk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iainvblk[i];

	for (i=0;i<_nlisttree;i++) {
		iblk = inv2blk[i];
		k = iptr[iblk];
		jainvblk2ind[k] = i;
		iptr[iblk]++;
	};

// Split column indices into blocks

	int *col2blk;

	col2blk = new int [nzja];
	if(!col2blk) MemoryFail(funcname);

	iblkprev = 0;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			col2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			col2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iacolblk;
	int *jacolblk2ind;

	iacolblk = new int [nproc+1];
	if(!iacolblk) MemoryFail(funcname);
	jacolblk2ind = new int [nzja];
	if(!jacolblk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacolblk[i] = 0;

	for (i=0;i<nzja;i++) {
		iblk = col2blk[i];
		iacolblk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iacolblk[i+1] = iacolblk[i]+iacolblk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iacolblk[i];

	for (i=0;i<nzja;i++) {
		iblk = col2blk[i];
		k = iptr[iblk];
		jacolblk2ind[k] = i;
		iptr[iblk]++;
	};

// Perform initial filling of the ordering array

	int *janew;

	janew = new int [nzja];
	if(!janew) MemoryFail(funcname);

	for (i=0;i<nzja;i++) janew[i] = -1;

	int j, ind, irow, irowold, ibs;

	for (iblk=0;iblk<nproc;iblk++) {
		icycle++;
		for (j=iablk[iblk];j<iablk[iblk+1];j++) {
			ind = jablk2ind[j];
			irowold = list[ind];
			irow = _listord[ind];
			ibs = irowold-blks[iblk];
			imaskloc[ibs] = icycle;
			listloc[ibs] = irow;
		};
		for (j=iainvblk[iblk];j<iainvblk[iblk+1];j++) {
			ind = jainvblk2ind[j];
			irow = _listtree[ind];
			irowold = _invordertree[ind];
			ibs = irowold-blks[iblk];
			imaskloc[ibs] = icycle;
			listloc[ibs] = irow;
		};
		for (j=iacolblk[iblk];j<iacolblk[iblk+1];j++) {
			ind = jacolblk2ind[j];
			irowold = ja[ind];
			ibs = irowold-blks[iblk];
			if (imaskloc[ibs] == icycle) {
				janew[ind] = listloc[ibs];
			};
		};
	};

// Free work arrays

	delete [] blks;
	delete [] listloc;
	delete [] imaskloc;
	delete [] list2blk;
	delete [] iablk;
	delete [] iptr;
	delete [] jablk2ind;
	delete [] inv2blk;
	delete [] iainvblk;
	delete [] jainvblk2ind;
	delete [] col2blk;
	delete [] iacolblk;
	delete [] jacolblk2ind;

// Create the list of remaining indices

	int nzjafin = 0;

	for (i=0;i<nzja;i++) {
		if (janew[i] == -1) nzjafin++;
	};

	int *listfin;
	int *fin2ind;

	listfin = new int [nzjafin];
	if(!listfin) MemoryFail(funcname);
	fin2ind = new int [nzjafin];
	if(!fin2ind) MemoryFail(funcname);

	nzjafin = 0;

	for (i=0;i<nzja;i++) {
		if (janew[i] == -1) {
			listfin[nzjafin] = ja[i];
			fin2ind[nzjafin] = i;
			nzjafin++;
		};
	};

// Create global partitionings

	int nptot, *ibegtot, *iendtot, *icputot;

	CSMatrix::ComputeIntervals (_comm,
											_nparts, _ibegarr, _iendarr,
											nptot, ibegtot, iendtot, icputot);

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrix::ChangeIndicesViaIntervals: intervals from different cpu's intersect ";
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	int nlistpart = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = nlistpart;
			nlistpart += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [nzjafin];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [nzjafin];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<nzjafin;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ipprev, ip, iproc;

	if (nptot > 0) {

		ipprev = 0;

		for (i=0;i<nzjafin;i++) {
			irow = listfin[i];
			ipbeg = 0;
			ipend = nptot-1;
			ip = ipprev;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[i] = icputot[ip];
				interval2cpu[i] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[i] = icputot[ipbeg];
				interval2cpu[i] = ipbeg;
				ipprev = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<nzjafin;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsarr;
	int **plists;

	ibsarr = new int [nproc];
	if(!ibsarr) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsarr[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		plists[i] = new int [niloc*3];
		if(!plists[i]) MemoryFail(funcname);
	};

	int *piarr;

	for (i=0;i<nzjafin;i++) {
		iproc = list2cpu[i];
		ibs = ibsarr[iproc];
		piarr = plists[iproc];
		piarr[ibs] = listfin[i];
		ibs++;
		piarr[ibs] = fin2ind[i];
		ibs++;
		piarr[ibs] = interval2cpu[i];
		ibs++;
		ibsarr[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = 3 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrix::ChangeIndicesViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrix::ChangeIndicesViaIntervals: wrong number of the received objects ";

// Store the result

	int iv, indnew;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*3);
		piarr = (int *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*3];
			iv = piarr[j*3+2];
			ip = ibsinterval[iv]+ind-ibegtot[iv];
			indnew = _orderpart[ip];
			piarr[j*3+2] = indnew;
		};
	};

// Send back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrix::ChangeIndicesViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

	if (NObjSend != nproc) throw " CSMatrix::ChangeIndicesViaIntervals: wrong number of the received objects ";

// Use received data

	int irownew;

	for (i=0;i<NObjSend;i++) {
		ni = ObjSizeSend[i] / (sizeof(int)*3);
		piarr = (int *)ObjSend[i];
		for (j=0;j<ni;j++) {
			irow = piarr[j*3];
			ind = piarr[j*3+1];
			irownew = piarr[j*3+2];
			janew[ind] = irownew;
		};
	};

// Free received data

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Free work arrays

	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] listfin;
	delete [] fin2ind;
	delete [] ibsinterval;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsarr;
	delete [] plists;

// Finalize renumbering

	for (i=0;i<nzja;i++) ja[i] = janew[i];
	for (i=0;i<nlist;i++) list[i] = _listord[i];

	delete [] janew;

// Sort matrix data

	SortListAndColumns ();

};

// Author: Kharchenko S.A.
// CSMatrix: Change rows/column indices
//========================================================================================
void CSMatrixR::ChangeIndicesViaIntervals (CMPIComm &_comm, // Change rows/column indices
																	int *_listord,
																	int _nlisttree, int *_listtree, int *_invordertree, 
																	int _nparts, int *_ibegarr, int *_iendarr, int *_orderpart) {

	const char *funcname = "ChangeIndicesViaIntervals";

	int nproc = _comm.GetNproc ();
	int myid = _comm.GetMyid ();

// Partition the total indices space

	int ntot = 0;

	int i, jj;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj > ntot) ntot = jj;
	};

	CMPIExchange::ExchangeArrayMPI (_comm,
												INTEGERVALUE, MAXIMUM, 1, &ntot, &ntot);

	ntot++;

	int ni = ntot / nproc;

	int *blks;

	blks = new int [nproc+1];
	if(!blks) MemoryFail(funcname);

	for (i=0;i<nproc;i++) blks[i] = i*ni;

	blks[nproc] = ntot;

// Create the mask

	int nimax = 0;

	for (i=0;i<nproc;i++) {
		ni = blks[i+1]-blks[i];
		if (ni > nimax) nimax = ni;
	};

	int *imaskloc;
	int *listloc;

	imaskloc = new int [nimax];
	if(!imaskloc) MemoryFail(funcname);
	listloc = new int [nimax];
	if(!listloc) MemoryFail(funcname);

	int icycle = -1;

	for (i=0;i<nimax;i++) imaskloc[i] = icycle;

// Split list data into the blocks

	int *list2blk;

	list2blk = new int [nlist];
	if(!list2blk) MemoryFail(funcname);

	int iblkprev = 0;

	int iblk, iblkbeg, iblkend;

	for (i=0;i<nlist;i++) {
		jj = list[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			list2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			list2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iablk;
	int *iptr;
	int *jablk2ind;

	iablk = new int [nproc+1];
	if(!iablk) MemoryFail(funcname);
	iptr = new int [nproc];
	if(!iptr) MemoryFail(funcname);
	jablk2ind = new int [nlist];
	if(!jablk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iablk[i] = 0;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		iablk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iablk[i+1] = iablk[i]+iablk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iablk[i];

	int k;

	for (i=0;i<nlist;i++) {
		iblk = list2blk[i];
		k = iptr[iblk];
		jablk2ind[k] = i;
		iptr[iblk]++;
	};

// Split inverse ordering data into blocks

	int *inv2blk;

	inv2blk = new int [_nlisttree];
	if(!inv2blk) MemoryFail(funcname);

	iblkprev = 0;

	for (i=0;i<_nlisttree;i++) {
		jj = _invordertree[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			inv2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			inv2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iainvblk;
	int *jainvblk2ind;

	iainvblk = new int [nproc+1];
	if(!iainvblk) MemoryFail(funcname);
	jainvblk2ind = new int [_nlisttree];
	if(!jablk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iainvblk[i] = 0;

	for (i=0;i<_nlisttree;i++) {
		iblk = inv2blk[i];
		iainvblk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iainvblk[i+1] = iainvblk[i]+iainvblk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iainvblk[i];

	for (i=0;i<_nlisttree;i++) {
		iblk = inv2blk[i];
		k = iptr[iblk];
		jainvblk2ind[k] = i;
		iptr[iblk]++;
	};

// Split column indices into blocks

	int *col2blk;

	col2blk = new int [nzja];
	if(!col2blk) MemoryFail(funcname);

	iblkprev = 0;

	for (i=0;i<nzja;i++) {
		jj = ja[i];
		if (jj >= blks[iblkprev] && jj < blks[iblkprev+1]) {
			col2blk[i] = iblkprev;
		} else {
			iblkbeg = 0;
			iblkend = nproc-1;
			while (iblkbeg != iblkend) {
				iblk = (iblkbeg+iblkend)/2;
				if (jj >= blks[iblk] && jj < blks[iblk+1]) {
					iblkbeg = iblk;
					iblkend = iblk;
					break;
				} else if (jj < blks[iblk]) {
					iblkend = iblk-1;
				} else if (jj >= blks[iblk+1]) {
					iblkbeg = iblk+1;
				};
			};
			col2blk[i] = iblkbeg;
			iblkprev = iblkbeg;
		};
	};

	int *iacolblk;
	int *jacolblk2ind;

	iacolblk = new int [nproc+1];
	if(!iacolblk) MemoryFail(funcname);
	jacolblk2ind = new int [nzja];
	if(!jacolblk2ind) MemoryFail(funcname);

	for (i=0;i<=nproc;i++) iacolblk[i] = 0;

	for (i=0;i<nzja;i++) {
		iblk = col2blk[i];
		iacolblk[iblk+1]++;
	};

	for (i=0;i<nproc;i++) iacolblk[i+1] = iacolblk[i]+iacolblk[i+1];

	for (i=0;i<nproc;i++) iptr[i] = iacolblk[i];

	for (i=0;i<nzja;i++) {
		iblk = col2blk[i];
		k = iptr[iblk];
		jacolblk2ind[k] = i;
		iptr[iblk]++;
	};

// Perform initial filling of the ordering array

	int *janew;

	janew = new int [nzja];
	if(!janew) MemoryFail(funcname);

	for (i=0;i<nzja;i++) janew[i] = -1;

	int j, ind, irow, irowold, ibs;

	for (iblk=0;iblk<nproc;iblk++) {
		icycle++;
		for (j=iablk[iblk];j<iablk[iblk+1];j++) {
			ind = jablk2ind[j];
			irowold = list[ind];
			irow = _listord[ind];
			ibs = irowold-blks[iblk];
			imaskloc[ibs] = icycle;
			listloc[ibs] = irow;
		};
		for (j=iainvblk[iblk];j<iainvblk[iblk+1];j++) {
			ind = jainvblk2ind[j];
			irow = _listtree[ind];
			irowold = _invordertree[ind];
			ibs = irowold-blks[iblk];
			imaskloc[ibs] = icycle;
			listloc[ibs] = irow;
		};
		for (j=iacolblk[iblk];j<iacolblk[iblk+1];j++) {
			ind = jacolblk2ind[j];
			irowold = ja[ind];
			ibs = irowold-blks[iblk];
			if (imaskloc[ibs] == icycle) {
				janew[ind] = listloc[ibs];
			};
		};
	};

// Free work arrays

	delete [] blks;
	delete [] listloc;
	delete [] imaskloc;
	delete [] list2blk;
	delete [] iablk;
	delete [] iptr;
	delete [] jablk2ind;
	delete [] inv2blk;
	delete [] iainvblk;
	delete [] jainvblk2ind;
	delete [] col2blk;
	delete [] iacolblk;
	delete [] jacolblk2ind;

// Create the list of remaining indices

	int nzjafin = 0;

	for (i=0;i<nzja;i++) {
		if (janew[i] == -1) nzjafin++;
	};

	int *listfin;
	int *fin2ind;

	listfin = new int [nzjafin];
	if(!listfin) MemoryFail(funcname);
	fin2ind = new int [nzjafin];
	if(!fin2ind) MemoryFail(funcname);

	nzjafin = 0;

	for (i=0;i<nzja;i++) {
		if (janew[i] == -1) {
			listfin[nzjafin] = ja[i];
			fin2ind[nzjafin] = i;
			nzjafin++;
		};
	};

// Create global partitionings

	int nptot, *ibegtot, *iendtot, *icputot;

	CSMatrix::ComputeIntervals (_comm,
											_nparts, _ibegarr, _iendarr,
											nptot, ibegtot, iendtot, icputot);

// Check that all intervals are correct (do not intersect)

	for (i=0;i<nptot-1;i++) {
		if (iendtot[i] >= ibegtot[i+1]) throw " CSMatrixR::ChangeIndicesViaIntervals: intervals from different cpu's intersect ";
	};

// Create the output data

	int *ibsinterval;

	ibsinterval = new int [nptot];
	if(!ibsinterval) MemoryFail(funcname);

	int nlistpart = 0;

	for (i=0;i<nptot;i++) {
		if (icputot[i] == myid) {
			ibsinterval[i] = nlistpart;
			nlistpart += (iendtot[i]-ibegtot[i]+1);
		} else {
			ibsinterval[i] = -1;
		};
	};

// For local list data find to which processor the data should be sent

	int *list2cpu;
	int *interval2cpu;

	list2cpu = new int [nzjafin];
	if(!list2cpu) MemoryFail(funcname);
	interval2cpu = new int [nzjafin];
	if(!interval2cpu) MemoryFail(funcname);

	for (i=0;i<nzjafin;i++) list2cpu[i] = -1;

	int ipbeg, ipend, ipprev, ip, iproc;

	if (nptot > 0) {

		ipprev = 0;

		for (i=0;i<nzjafin;i++) {
			irow = listfin[i];
			ipbeg = 0;
			ipend = nptot-1;
			ip = ipprev;
			if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
				list2cpu[i] = icputot[ip];
				interval2cpu[i] = ip;
			} else {
				while (true) {
					if (ipbeg > ipend) break;
					ip = (ipbeg+ipend) / 2;
					if (ip < ipbeg) ip = ipbeg;
					if (ip > ipend) ip = ipend;
					if (irow >= ibegtot[ip] && irow <= iendtot[ip]) {
						ipbeg = ip;
						break;
					};
					if (ipbeg == ipend) break;
					if (irow > iendtot[ip]) ipbeg = ip+1;
					if (irow < ibegtot[ip]) ipend = ip-1;
				};
				list2cpu[i] = icputot[ipbeg];
				interval2cpu[i] = ipbeg;
				ipprev = ipbeg;
			};
		};

	};

// Prepare the data for send

	int *nlist2cpu;

	nlist2cpu = new int [nproc];
	if(!nlist2cpu) MemoryFail(funcname);

	for (i=0;i<nproc;i++) nlist2cpu[i] = 0;

	for (i=0;i<nzjafin;i++) {
		iproc = list2cpu[i];
		if (iproc >= 0) {
			nlist2cpu[iproc]++;
		};
	};

	int *ibsarr;
	int **plists;

	ibsarr = new int [nproc];
	if(!ibsarr) MemoryFail(funcname);
	plists = new int * [nproc];
	if(!plists) MemoryFail(funcname);

	for (i=0;i<nproc;i++) ibsarr[i] = 0;

	int niloc;

	for (i=0;i<nproc;i++) {
		niloc = nlist2cpu[i];
		plists[i] = new int [niloc*3];
		if(!plists[i]) MemoryFail(funcname);
	};

	int *piarr;

	for (i=0;i<nzjafin;i++) {
		iproc = list2cpu[i];
		ibs = ibsarr[iproc];
		piarr = plists[iproc];
		piarr[ibs] = listfin[i];
		ibs++;
		piarr[ibs] = fin2ind[i];
		ibs++;
		piarr[ibs] = interval2cpu[i];
		ibs++;
		ibsarr[iproc] = ibs;
	};

// Exchange the data

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
		ObjSizeSend[i] = 3 * sizeof(int) * nlist2cpu[i];
		ObjSend[i] = (char *)plists[i];
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
	if(info) throw " CSMatrixR::ChangeIndicesViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	delete [] ObjSend;

	for (i=0;i<nproc;i++) delete [] plists[i];

	if (NObjRecv != nproc) throw " CSMatrixR::ChangeIndicesViaIntervals: wrong number of the received objects ";

// Store the result

	int iv, indnew;

	for (i=0;i<NObjRecv;i++) {
		ni = ObjSizeRecv[i] / (sizeof(int)*3);
		piarr = (int *)ObjRecv[i];
		for (j=0;j<ni;j++) {
			ind = piarr[j*3];
			iv = piarr[j*3+2];
			ip = ibsinterval[iv]+ind-ibegtot[iv];
			indnew = _orderpart[ip];
			piarr[j*3+2] = indnew;
		};
	};

// Send back

	info = CMPIExchange::DataExchangeMPI (_comm,
														NObjRecv, ObjTypeRecv, ObjIDRecv, CpuIDRecv,
														ObjSizeRecv, ObjRecv,
														NObjSend, ObjTypeSend, ObjIDSend, CpuIDSend,
														ObjSizeSend, ObjSend);
	if(info) throw " CSMatrixR::ChangeIndicesViaIntervals: Error in DataExchangeMPI";

	delete [] ObjTypeRecv;
	delete [] ObjIDRecv;
	delete [] CpuIDRecv;
	delete [] ObjSizeRecv;
	for(i = 0; i < NObjRecv; i++) {
		delete [] ObjRecv[i];
	};
	delete [] ObjRecv;

	if (NObjSend != nproc) throw " CSMatrixR::ChangeIndicesViaIntervals: wrong number of the received objects ";

// Use received data

	int irownew;

	for (i=0;i<NObjSend;i++) {
		ni = ObjSizeSend[i] / (sizeof(int)*3);
		piarr = (int *)ObjSend[i];
		for (j=0;j<ni;j++) {
			irow = piarr[j*3];
			ind = piarr[j*3+1];
			irownew = piarr[j*3+2];
			janew[ind] = irownew;
		};
	};

// Free received data

	delete [] ObjTypeSend;
	delete [] ObjIDSend;
	delete [] CpuIDSend;
	delete [] ObjSizeSend;
	for(i = 0; i < NObjSend; i++) {
		delete [] ObjSend[i];
	};
	delete [] ObjSend;

// Free work arrays

	delete [] ibegtot;
	delete [] iendtot;
	delete [] icputot;
	delete [] listfin;
	delete [] fin2ind;
	delete [] ibsinterval;
	delete [] list2cpu;
	delete [] interval2cpu;
	delete [] nlist2cpu;
	delete [] ibsarr;
	delete [] plists;

// Finalize renumbering

	for (i=0;i<nzja;i++) ja[i] = janew[i];
	for (i=0;i<nlist;i++) list[i] = _listord[i];

	delete [] janew;

// Sort matrix data

	SortListAndColumns ();

};

// Author: Kharchenko S.A.
// CSMatrix: Compute ordering of the nodes into two parts with the boundary on blocks surfaces
//========================================================================================
void CSMatrix::Compute2LevelOrderingWithBlockSurfaceBoundary (int _nblks, int *_blks, // Compute ordering of the nodes into two parts with the boundary on blocks surfaces
																					int *_orderblk, int *_blkso_nobnd, int *_order, 
																					int &_ind1blk, int &_ind2blk,
																					int &_ind1, int &_ind2) const {

	const char *funcname = "Compute2LevelOrderingWithBlockSurfaceBoundary";

// Check data on entry

	if (_blks[_nblks] != n) {
		throw " CSMatrix::Compute2LevelOrderingWithBlockSurfaceBoundary: incorrect block partitioning ";
	};

// Compute initial Nd ordering

//	A2Ps (10,"MtrIni.ps",0,&_nblks);

	int ordtype = -1;

	OrderPrfMtr (ordtype, _order);

// Reorder and find initial separators

	CSMatrix mtro;

	mtro = this->OrdMtr (_order);

//	mtro.A2Ps (10,"MtrOrd.ps",0,&_nblks);

	int ind1, ind2;

	mtro.FindSeparator (ind1, ind2);

	if (ind1 < 0 || ind2 < 0) {
		ind1 = 0;
		ind2 = 0;
	};

// Compute inverse ordering and ind2blk

	int *iorder;
	int *ind2blk;

	iorder = new int [n];
	if(!iorder) MemoryFail(funcname);
	ind2blk = new int [n];
	if(!ind2blk) MemoryFail(funcname);

	int i, j;

	for (i=0;i<n;i++) iorder[_order[i]] = i;

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			ind2blk[j] = i;
		};
	};

// Mark all blocks according to the separators

	int *markblk;

	markblk = new int [_nblks];
	if(!markblk) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) markblk[i] = -2;

	int inew, iold, iblkold;

	for (inew=0;inew<ind1;inew++) {
		iold = iorder[inew];
		iblkold = ind2blk[iold];
		if (markblk[iblkold] == -2) {
			markblk[iblkold] = 0;
		} else if (markblk[iblkold] == 1) {
			markblk[iblkold] = -1;
		};
	};

	for (inew=ind1;inew<ind2;inew++) {
		iold = iorder[inew];
		iblkold = ind2blk[iold];
		if (markblk[iblkold] == -2) {
			markblk[iblkold] = 1;
		} else if (markblk[iblkold] == 0) {
			markblk[iblkold] = -1;
		};
	};

	for (inew=ind2;inew<n;inew++) {
		iold = iorder[inew];
		iblkold = ind2blk[iold];
		markblk[iblkold] = -1;
	};

// Compute block sparsity of the initial matrix

	int *sp2blkloc;

	sp2blkloc = new int [n];
	if(!sp2blkloc) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) {
		for (j=_blks[i];j<_blks[i+1];j++) {
			sp2blkloc[j] = i;
		};
	};

	CSMatrix aspblk;

	aspblk = this->BlockSparsity (_nblks, sp2blkloc);

	delete [] sp2blkloc;

// Extend the boundary according to the parameter

	int ncycleext = 1;

	int *listblk;

	listblk = new int [_nblks];
	if(!listblk) MemoryFail(funcname);

	int nlistblk = 0;

	for (i=0;i<_nblks;i++) {
		if (markblk[i] < 0) {
			listblk[nlistblk] = i;
			nlistblk++;
		};
	};

	int *piablk = aspblk.GetIa ();
	int *pjablk = aspblk.GetJa ();

	int icycleext, nlistblk1, iblk, jblk;

	for (icycleext=0;icycleext<ncycleext;icycleext++) {
		nlistblk1 = nlistblk;
		for (i=0;i<nlistblk;i++) {
			iblk = listblk[i];
			for (j=piablk[iblk];j<piablk[iblk+1];j++) {
				jblk = pjablk[j];
				if (markblk[jblk] >= 0) {
					listblk[nlistblk1] = jblk;
					nlistblk1++;
					markblk[jblk] = -1;
				};
			};
		};
		nlistblk = nlistblk1;
	};

	delete [] listblk;

// Compute block ordering

	int nblks1=0;
	int nblks2=0;
	int nblksbnd=0;

	for (i=0;i<_nblks;i++) {
		if (markblk[i] == 0) nblks1++;
		if (markblk[i] == 1) nblks2++;
		if (markblk[i] == -1) nblksbnd++;
	};

	int ibs1 = 0;
	int ibs2 = ibs1+nblks1;
	int ibsbnd = ibs2+nblks2;

	int *orderblk;
	int *iorderblk;

	orderblk = new int [_nblks];
	if(!orderblk) MemoryFail(funcname);
	iorderblk = new int [_nblks];
	if(!iorderblk) MemoryFail(funcname);

	for (i=0;i<_nblks;i++) {
		if (markblk[i] == 0) {
			orderblk[i] = ibs1;
			ibs1++;
		};
		if (markblk[i] == 1) {
			orderblk[i] = ibs2;
			ibs2++;
		};
		if (markblk[i] == -1) {
			orderblk[i] = ibsbnd;
			ibsbnd++;
		};
	};

	for (i=0;i<_nblks;i++) iorderblk[orderblk[i]] = i;

// Recompute point ordering and new blocks partitioning

	int *blkso;
	int *orderloc;

	blkso = new int [_nblks+1];
	if(!blkso) MemoryFail(funcname);
	orderloc = new int [n];
	if(!orderloc) MemoryFail(funcname);

	blkso[0] = 0;

	for (i=0;i<_nblks;i++) {
		iold = iorderblk[i];
		blkso[i+1] = (_blks[iold+1]-_blks[iold]);
	};

	for (i=0;i<_nblks;i++) blkso[i+1] = blkso[i]+blkso[i+1];

	int ibs, ibsn, ni;

	for (i=0;i<_nblks;i++) {
		inew = orderblk[i];
		ibs = _blks[i];
		ibsn = blkso[inew];
		ni = _blks[i+1]-_blks[i];
		for (j=0;j<ni;j++) {
			orderloc[ibs+j] = ibsn+j;
		};
	};

// Reorder the matrix

	mtro = this->OrdMtr (orderloc);

// Recompute blocks partitioning

	int nblksloc = nblksbnd+2;

	int *blksobnd;

	blksobnd = new int [nblksloc+1];
	if(!blksobnd) MemoryFail(funcname);

	blksobnd[0] = 0;
	blksobnd[1] = blkso[nblks1];
	blksobnd[2] = blkso[nblks1+nblks2];

	for (i=0;i<nblksbnd;i++) blksobnd[i+3] = blksobnd[i+2]+(blkso[nblks1+nblks2+i+1]-blkso[nblks1+nblks2+i]);

// Compute weights data

	sp2blkloc = new int [n];
	if(!sp2blkloc) MemoryFail(funcname);

	for (i=0;i<nblksloc;i++) {
		for (j=blksobnd[i];j<blksobnd[i+1];j++) {
			sp2blkloc[j] = i;
		};
	};

	aspblk = mtro.BlockSparsityDiag (nblksloc, sp2blkloc);

	int nzjablk = aspblk.GetNzja ();

	int *jaweights;

	jaweights = new int [nzjablk];
	if(!jaweights) MemoryFail(funcname);

	for (i=0;i<nzjablk;i++) jaweights[i] = 0;

	int *nzblkarr;

	nzblkarr = new int [nblksloc];
	if(!nzblkarr) MemoryFail(funcname);

	piablk = aspblk.GetIa ();
	pjablk = aspblk.GetJa ();
	int *piamtr = mtro.GetIa ();
	int *pjamtr = mtro.GetJa ();

	int jj;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (j=piablk[iblk];j<piablk[iblk+1];j++) {
			jblk = pjablk[j];
			nzblkarr[jblk] = 0;
		};
		for (i=blksobnd[iblk];i<blksobnd[iblk+1];i++) {
			for (j=piamtr[i];j<piamtr[i+1];j++) {
				jj = pjamtr[j];
				jblk = sp2blkloc[jj];
				if (iblk != jblk) nzblkarr[jblk]++;
			};
		};
		for (j=piablk[iblk];j<piablk[iblk+1];j++) {
			jblk = pjablk[j];
			jaweights[j] = nzblkarr[jblk];
		};
	};

	int *iablkflt;
	int *jablkflt;
	int *ndweights;
	int *adjweights;
	int *partition;

	iablkflt = new int [nblksloc+1];
	if(!iablkflt) MemoryFail(funcname);
	jablkflt = new int [nzjablk];
	if(!jablkflt) MemoryFail(funcname);
	ndweights = new int [nblksloc];
	if(!ndweights) MemoryFail(funcname);
	adjweights = new int [nzjablk];
	if(!adjweights) MemoryFail(funcname);
	partition = new int [nblksloc];
	if(!partition) MemoryFail(funcname);

	iablkflt[0] = 0;
	nzjablk = 0;

	for (iblk=0;iblk<nblksloc;iblk++) {
		for (j=piablk[iblk];j<piablk[iblk+1];j++) {
			jblk = pjablk[j];
			if (iblk != jblk) {
				jablkflt[nzjablk] = pjablk[j];
				adjweights[nzjablk] = jaweights[j];
				nzjablk++;
			};
		};
		iablkflt[iblk+1] = nzjablk;
		ndweights[iblk] = blksobnd[iblk+1]-blksobnd[iblk];
	};

// Partition

	int pzero = 0;
	int ptwo = 2;
	int pthree = 3;
	int pvol;

	int options[5];

	options[0] = 0;

	METIS_PartGraphRecursive (&nblksloc,iablkflt, jablkflt, 
										ndweights, adjweights, 
										&pthree, &pzero, &ptwo, 
										options, &pvol, partition);

// Store partitioning

	int *partitionext;

	partitionext = new int [_nblks];
	if(!partitionext) MemoryFail(funcname);

	for (i=0;i<nblks1;i++) partitionext[i] = partition[0];
	for (i=0;i<nblks2;i++) partitionext[nblks1+i] = partition[1];
	for (i=0;i<nblksbnd;i++) partitionext[nblks1+nblks2+i] = partition[i+2];

	for (i=0;i<_nblks;i++) {
		iblkold = iorderblk[i];
		markblk[iblkold] = partitionext[i];
	};

// Free current work arrays

	delete [] sp2blkloc;
	delete [] orderblk;
	delete [] iorderblk;
	delete [] blkso;
	delete [] orderloc;
	delete [] blksobnd;
	delete [] jaweights;
	delete [] nzblkarr;
	delete [] iablkflt;
	delete [] jablkflt;
	delete [] ndweights;
	delete [] adjweights;
	delete [] partition;
	delete [] partitionext;

// Determine the two-sided boundaries

	int *imasknd;
	int *listloc;

	imasknd = new int [n];
	if(!imasknd) MemoryFail(funcname);
	listloc = new int [n];
	if(!listloc) MemoryFail(funcname);

	int ipart;

	for (iblk=0;iblk<_nblks;iblk++) {
		ipart = markblk[iblk];
		for (i=_blks[iblk];i<_blks[iblk+1];i++) {
			imasknd[i] = ipart;
		};
	};

	int nlistloc = 0;

	this->OptimalBoundary (imasknd,
									nlistloc, listloc);

	for (i=0;i<nlistloc;i++) {
		jj = listloc[i];
		imasknd[jj] = -1;
	};

// Create the final ordering array

	int *list1loc;
	int *list2loc;
	int *listbnd;

	list1loc = new int [n];
	if(!list1loc) MemoryFail(funcname);
	list2loc = new int [n];
	if(!list2loc) MemoryFail(funcname);
	listbnd = new int [n];
	if(!listbnd) MemoryFail(funcname);

	int nlist1loc = 0;
	int nlist2loc = 0;
	int nlistbnd = 0;

	int itype;

	for (i=0;i<n;i++) {
		itype = imasknd[i];
		if (itype == 0) {
			list1loc[nlist1loc] = i;
			nlist1loc++;
		} else if (itype == 1) {
			list2loc[nlist2loc] = i;
			nlist2loc++;
		} else {
			listbnd[nlistbnd] = i;
			nlistbnd++;
		};
	};

	qsort (list1loc, nlist1loc, sizeof(int), compint);
	qsort (list2loc, nlist2loc, sizeof(int), compint);

	_ind1 = nlist1loc;
	_ind2 = nlist1loc+nlist2loc;
//	cout << " Ind1 = " << _ind1 << " Ind2 = " << _ind2 << endl;

	for (i=0;i<nlist1loc;i++) {
		jj = list1loc[i];
		_order[jj] = i;
	};
	for (i=0;i<nlist2loc;i++) {
		jj = list2loc[i];
		_order[jj] = nlist1loc+i;
	};
	for (i=0;i<nlistbnd;i++) {
		jj = listbnd[i];
		_order[jj] = nlist1loc+nlist2loc+i;
	};

// Create the final block ordering array

	int *listblk1;
	int *listblk2;

	listblk1 = new int [_nblks];
	if(!listblk1) MemoryFail(funcname);
	listblk2 = new int [_nblks];
	if(!listblk2) MemoryFail(funcname);

	nblks1 = 0;
	nblks2 = 0;

	for (i=0;i<_nblks;i++) {
		itype = markblk[i];
		if (itype == 0) {
			listblk1[nblks1] = i;
			nblks1++;
		} else if (itype == 1) {
			listblk2[nblks2] = i;
			nblks2++;
		};
	};

	qsort (listblk1, nblks1, sizeof(int), compint);
	qsort (listblk2, nblks2, sizeof(int), compint);

// Check block ordering and point ordering

	for (i=0;i<nblks1;i++) {
		iblk = listblk1[i];
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (_order[j] >= _ind1 && _order[j] < _ind2) {
				throw " CSMatrix::Compute2LevelOrderingWithBlockSurfaceBoundary: error in point ordering ";
			};
		};
	};

	for (i=0;i<nblks2;i++) {
		iblk = listblk2[i];
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (_order[j] < _ind1) {
				throw " CSMatrix::Compute2LevelOrderingWithBlockSurfaceBoundary: error in point ordering ";
			};
		};
	};

// Modify the numbers of blocks in first and second parts (move empty blocks into the boundary)

	int nblks1new = 0;
	int nblks2new = 0;
	nblksbnd = 0;

	int icheck;

	for (i=0;i<nblks1;i++) {
		iblk = listblk1[i];
		icheck = -1;
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (_order[j] < _ind2) icheck = 0;
		};
		if (icheck == 0) {
			nblks1new++;
		} else {
			nblksbnd++;
		};
	};

	for (i=0;i<nblks2;i++) {
		iblk = listblk2[i];
		icheck = -1;
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (_order[j] < _ind2) icheck = 0;
		};
		if (icheck == 0) {
			nblks2new++;
		} else {
			nblksbnd++;
		};
	};

// Prepare blocks ordering and partitioning arrays

	for (i=0;i<=_nblks;i++) _blkso_nobnd[i] = 0;

	_ind1blk = nblks1new;
	_ind2blk = nblks1new+nblks2new;

	nblks1new = 0;
	nblks2new = _ind1blk;
	nblksbnd = _ind2blk;

	int niloc;

	for (i=0;i<nblks1;i++) {
		iblk = listblk1[i];
		icheck = -1;
		niloc = 0;
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (_order[j] < _ind2) {
				icheck = 0;
				niloc++;
			};
		};
		if (icheck == 0) {
			_blkso_nobnd[nblks1new+1] = niloc;
			_orderblk[iblk] = nblks1new;
			nblks1new++;
		} else {
			_orderblk[iblk] = nblksbnd;
			nblksbnd++;
		};
	};

	for (i=0;i<nblks2;i++) {
		iblk = listblk2[i];
		icheck = -1;
		niloc = 0;
		for (j=_blks[iblk];j<_blks[iblk+1];j++) {
			if (_order[j] < _ind2) {
				icheck = 0;
				niloc++;
			};
		};
		if (icheck == 0) {
			_blkso_nobnd[nblks2new+1] = niloc;
			_orderblk[iblk] = nblks2new;
			nblks2new++;
		} else {
			_orderblk[iblk] = nblksbnd;
			nblksbnd++;
		};
	};

	for (i=0;i<_nblks;i++) _blkso_nobnd[i+1] = _blkso_nobnd[i]+_blkso_nobnd[i+1];

//	cout << " IndBlk1 = " << _blkso_nobnd[_ind1blk] << " IndBlk2 = " << _blkso_nobnd[_ind2blk] << endl;

// Free work arrays

	delete [] iorder;
	delete [] ind2blk;
	delete [] markblk;
	delete [] imasknd;
	delete [] listloc;
	delete [] list1loc;
	delete [] list2loc;
	delete [] listbnd;
	delete [] listblk1;
	delete [] listblk2;

};

// Author: Kharchenko S.A.
// CSMatrix: Partition matrix according to the prescribed binary tree and local supernodes partitioning, internal matrix ordering is assumed to be ND
//========================================================================================
void CSMatrix::PartBinaryTreeNDSups (CTree &_tree, CTree &_treend, // Partition matrix according to the prescribed binary tree and local supernodes partitioning, internal matrix ordering is assumed to be ND
												int _nsups, int *_sprnds, 
												int *_ordersp, int *_ordernd,
												int &_nblkssp, int *&_blkssp, int *&_blk2cpusp, 
												int &_nblksnd, int *&_blksnd, int *&_blk2cpund) const {

	const char *funcname = "PartBinaryTreeNDSups";

// Compute symmetrized matrix

	CSMatrix as;

	as = SymmMtr ();

// Allocate work arrays

	int nlev = _tree.nlev;

	int *inodelv, *ichildlv, *nchildslv;
	int *iorder;
	int *iordersp;
	int *sprndsord;
	int *nzspord;

	inodelv = new int [nlev+1];
	if (!inodelv) MemoryFail (funcname);
	ichildlv = new int [nlev+1];
	if (!ichildlv) MemoryFail (funcname);
	nchildslv = new int [nlev+1];
	if (!nchildslv) MemoryFail (funcname);
	iorder = new int [nsupr];
	if (!iorder) MemoryFail (funcname);
	iordersp = new int [_nsups];
	if (!iordersp) MemoryFail (funcname);
	sprndsord = new int [_nsups+1];
	if (!sprndsord) MemoryFail (funcname);
	nzspord = new int [_nsups+1];
	if (!nzspord) MemoryFail (funcname);

	int i;

	for (i=0;i<_nsups+1;i++) sprndsord[i] = _sprnds[i];
	for (i=0;i<_nsups;i++) nzspord[i] = _sprnds[i+1]-_sprnds[i];

// Prepare root node data

	_treend = _tree;

	int ilevcurr = 0;

	int rootnode = _tree.rootid;

	_treend.nodes[rootnode].indbeg = 0;
	_treend.nodes[rootnode].indend = nsupr-1;
	_treend.nodes[rootnode].indbegtot = 0;
	_treend.nodes[rootnode].indendtot = nsupr-1;

	_tree.nodes[rootnode].indbeg = 0;
	_tree.nodes[rootnode].indend = _nsups-1;
	_tree.nodes[rootnode].indbegtot = 0;
	_tree.nodes[rootnode].indendtot = _nsups-1;

	inodelv[ilevcurr] = rootnode;
	ichildlv[ilevcurr] = 0;
	nchildslv[ilevcurr] = 1;

// Search the tree

	CSMatrix mtrdummy;

	int inodecurr, inodeloc, ichildnode, ifathernode, ichildcurr, nchildscurr;
	int indbeg, indend, nloc, nzloc;
	int isupbeg, isupend, j, jj, nsuploc;
	int ind1blk, ind2blk, ind1, ind2;

	for (i=0;i<nsupr;i++) _ordernd[i] = i;
	for (i=0;i<_nsups;i++) _ordersp[i] = i;

	while (ilevcurr >= 0) {

// Take data of the node

		inodecurr = inodelv[ilevcurr];
		ichildcurr = ichildlv[ilevcurr];

		indbeg = _treend.nodes[inodecurr].indbeg;
		indend = _treend.nodes[inodecurr].indend;

		isupbeg = _tree.nodes[inodecurr].indbeg;
		isupend = _tree.nodes[inodecurr].indend;

		nchildscurr = _tree.nodes[inodecurr].nchilds;

// Take current submatrix

		nloc = indend-indbeg+1;
		nzloc = 0;

		for (i=indbeg;i<=indend;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj >= indbeg && jj <= indend) nzloc++;
			};
		};

		CSMatrix aloc (nloc,nzloc);

		nloc = 0;
		nzloc = 0;
		aloc.ia[0] = 0;
		for (i=indbeg;i<=indend;i++) {
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj >= indbeg && jj <= indend) {
					aloc.ja[nzloc] = jj-indbeg;
					nzloc++;
				};
			};
			aloc.ia[nloc+1] = nzloc;
			nloc++;
		};
		aloc.m = nloc;
		aloc.n = nloc;
		aloc.nsupr = nloc;
		aloc.nsupc = nloc;
		aloc.nzja = nzloc;
		aloc.nlist = nloc;

// Compute new ordering

		nsuploc = isupend-isupbeg+1;

		int *sprndsloc;
		int *orderloc, *iorderloc;
		int *ordersploc, *iordersploc;
		int *sprndso_nobnd;

		sprndsloc = new int [nsuploc+1];
		if (!sprndsloc) MemoryFail (funcname);
		orderloc = new int [nloc];
		if (!orderloc) MemoryFail (funcname);
		iorderloc = new int [nloc];
		if (!iorderloc) MemoryFail (funcname);
		ordersploc = new int [nsuploc];
		if (!ordersploc) MemoryFail (funcname);
		iordersploc = new int [nsuploc];
		if (!iordersploc) MemoryFail (funcname);
		sprndso_nobnd = new int [nsuploc+1];
		if (!sprndso_nobnd) MemoryFail (funcname);

		sprndsloc[0] = 0;

		for (i=0;i<nsuploc;i++) {
//			sprndsloc[i+1] = sprndsord[isupbeg+i+1]-sprndsord[isupbeg+i];
			sprndsloc[i+1] = nzspord[isupbeg+i];
		};

		for (i=0;i<nsuploc;i++) sprndsloc[i+1] = sprndsloc[i]+sprndsloc[i+1];

		if (nchildscurr == 1) {
			int ordtyp = -1;
			aloc.OrderPrfMtr (ordtyp, orderloc);
			for (i=0;i<nsuploc;i++) ordersploc[i] = i;
			for (i=0;i<=nsuploc;i++) sprndso_nobnd[i] = sprndsloc[i];
			ind1blk = 0;
			ind2blk = 0;
			ind1 = 0;
			ind2 = 0;
		} else {
			aloc.Compute2LevelOrderingWithBlockSurfaceBoundary (nsuploc, sprndsloc,
																					ordersploc, sprndso_nobnd, orderloc, 
																					ind1blk, ind2blk,
																					ind1, ind2);
		};

		for (i=0;i<nloc;i++) iorderloc[orderloc[i]] = i;
		for (i=0;i<nsuploc;i++) iordersploc[ordersploc[i]] = i;

// Free submatrix data

		aloc = mtrdummy;

// Update the ordering

		int iold;

		for (i=indbeg;i<=indend;i++) {
			iold = iorderloc[i-indbeg];
			orderloc[i-indbeg] = _ordernd[iold+indbeg];
		};

		for (i=indbeg;i<=indend;i++) {
			_ordernd[i] = orderloc[i-indbeg];
		};

		for (i=0;i<nloc;i++) orderloc[iorderloc[i]] = i;

		for (i=isupbeg;i<=isupend;i++) {
			iold = iordersploc[i-isupbeg];
			ordersploc[i-isupbeg] = _ordersp[iold+isupbeg];
		};

		for (i=isupbeg;i<=isupend;i++) {
			_ordersp[i] = ordersploc[i-isupbeg];
		};

		for (i=0;i<nsuploc;i++) ordersploc[iordersploc[i]] = i;

// Update ordered blocks partitioning

		int iloc;

		for (i=isupbeg;i<=isupend;i++) {
			iloc = i-isupbeg;
			sprndsord[i+1] = sprndsord[i]+(sprndso_nobnd[iloc+1]-sprndso_nobnd[iloc]);
			nzspord[i] = sprndso_nobnd[iloc+1]-sprndso_nobnd[iloc];
		};

// Reorder current part of the matrix in place

		int *ialoc, *jaloc;

		nzloc = 0;
		for (i=indbeg;i<=indend;i++) {
			nzloc += as.ia[i+1]-as.ia[i];
		};

		ialoc = new int [nloc+1];
		if (!ialoc) MemoryFail (funcname);
		jaloc = new int [nzloc];
		if (!jaloc) MemoryFail (funcname);

		nzloc = 0;

		ialoc[0] = 0;

		for (i=indbeg;i<=indend;i++) {
			int nzloc0 = nzloc;
			for (j=as.ia[i];j<as.ia[i+1];j++) {
				jj = as.ja[j];
				if (jj>=indbeg && jj<=indend) {
					int jjloc = jj-indbeg;
					int jjnew = orderloc[jjloc];
					jj = jjnew + indbeg;
				};
				jaloc[nzloc] = jj;
				nzloc++;
			};
			qsort (jaloc+nzloc0, nzloc-nzloc0, sizeof(int), compint);
			ialoc[i-indbeg+1] = nzloc;
		};

		nzloc = as.ia[indbeg];

		for (i=indbeg;i<=indend;i++) {
			int iold = iorderloc[i-indbeg];
			for (j=ialoc[iold];j<ialoc[iold+1];j++) {
				jj = jaloc[j];
				as.ja[nzloc] = jj;
				nzloc++;
			};
			as.ia[i+1] = nzloc;
		};

		delete [] ialoc;
		delete [] jaloc;

// Register the partitioning in the child nodes if necessary

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			_treend.nodes[ichildnode].indbeg = indbeg;
			_treend.nodes[ichildnode].indend = indbeg+ind1-1;
			_treend.nodes[ichildnode].indbegtot = indbeg;
			_treend.nodes[ichildnode].indendtot = indbeg+ind1-1;

			_tree.nodes[ichildnode].indbeg = isupbeg;
			_tree.nodes[ichildnode].indend = isupbeg+ind1blk-1;
			_tree.nodes[ichildnode].indbegtot = isupbeg;
			_tree.nodes[ichildnode].indendtot = isupbeg+ind1blk-1;

			ichildnode = _tree.nodes[inodecurr].childs[1];

			_treend.nodes[ichildnode].indbeg = indbeg+ind1;
			_treend.nodes[ichildnode].indend = indbeg+ind2-1;
			_treend.nodes[ichildnode].indbegtot = indbeg+ind1;
			_treend.nodes[ichildnode].indendtot = indbeg+ind2-1;

			_tree.nodes[ichildnode].indbeg = isupbeg+ind1blk;
			_tree.nodes[ichildnode].indend = isupbeg+ind2blk-1;
			_tree.nodes[ichildnode].indbegtot = isupbeg+ind1blk;
			_tree.nodes[ichildnode].indendtot = isupbeg+ind2blk-1;

			_treend.nodes[inodecurr].indbeg = indbeg+ind2;
			_treend.nodes[inodecurr].indend = indend;

			_tree.nodes[inodecurr].indbeg = isupbeg+ind2blk;
			_tree.nodes[inodecurr].indend = isupend;

		};

// Free local arrays

		delete [] sprndsloc;
		delete [] orderloc;
		delete [] iorderloc;
		delete [] ordersploc;
		delete [] iordersploc;
		delete [] sprndso_nobnd;

// Go down level

		if (nchildscurr > 1) {

			ichildnode = _tree.nodes[inodecurr].childs[0];

			ilevcurr++;

			inodelv[ilevcurr] = ichildnode;
			ichildlv[ilevcurr] = 0;
			nchildslv[ilevcurr] = nchildscurr;

			goto exitcycle;

		};

// Look for the next child

nextchild:;

		if (ichildcurr < nchildslv[ilevcurr]-1) {
			if (ilevcurr == 0) {
				throw " CSMatrix::PartBinaryTreeNDSups: internal error";
			};
			ifathernode = _tree.nodes[inodecurr].fatherid;
			ichildcurr ++;
			inodeloc = _tree.nodes[ifathernode].childs[ichildcurr];
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

		if (ichildcurr == nchildslv[ilevcurr]-1) {
			goto uplevel;
		} else {
			goto nextchild;
		};

exitcycle:;

	};

// Inverse order

	for (i=0;i<nsupr;i++) iorder[_ordernd[i]] = i;
	for (i=0;i<nsupr;i++) _ordernd[i] = iorder[i];

	for (i=0;i<_nsups;i++) iordersp[_ordersp[i]] = i;
	for (i=0;i<_nsups;i++) _ordersp[i] = iordersp[i];

// Return blocks partitionings

	_tree.PartitionTree (_nblkssp, _blkssp, _blk2cpusp);

	_treend.PartitionTree (_nblksnd, _blksnd, _blk2cpund);

// Free work arrays

	delete [] inodelv;
	delete [] ichildlv;
	delete [] nchildslv;
	delete [] iorder;
	delete [] iordersp;
	delete [] sprndsord;
	delete [] nzspord;

};

#ifdef __FV_2004__
} // namespace EqnSolver
} // namespace flowvision
#endif
