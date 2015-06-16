//------------------------------------------------------------------------------------------------
// File: rcm.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
#endif

#include "globals.h"


////////////////////////////////////////////////////////////////////////////////////////////////
//...reordering algorithms (reverse Cuthill-Mckee ordering);
void genrcm(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask);
void subrcm(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n);
/////////////////////////
//...parameters:
//   n          - the dimension of the matrix;
//   xadj, iadj - the matrix structure: xadj[n+1], iadj[*]; information about row i is stored 
//                in xadj[i-1] -- xadj[i]-1 of the the adjacency structure iadj[*];
//                for each row, it contains the column indices of the nonzero entries;
//   perm[n]    - contains the rcm ordering;
//   mask[n]    - marks variables that have been numbered (working array);
//   xls[n+1]   - the index vector for a level structure; the level structure is stored 
//                in the currently unused spaces in the permutation vector perm;
//   nsubg      - the size of the subgraph;
//   subg[n]    - contains the nodes in subgraph (which may be disconnected);
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
//...generates the connected level structure rooted at a given node;
void rootls(int root, int * xadj, int * iadj, int * mask, int & nlvl, int * xls, int * ls, int n)
{
   int i, iccsze, j, jstop, jstrt, lbegin, lvlend, lvsize, nbr, node;

   mask[root-1] = 0;
   ls[0]        = root;
   nlvl         = 0;
   lvlend       = 0;
   iccsze       = 1;

   do {
       lbegin      = lvlend + 1;
       lvlend      = iccsze;
       xls[nlvl++] = lbegin;

/////////////////////////////////////////////////////////////////////////////////////////////
//...Generate the next level by finding all the masked neighbors of nodes in the current level;
       for (i = lbegin; i <= lvlend;  i++) {
            jstrt = xadj[(node = ls[i-1])-1];
            jstop = xadj[ node]-1;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 ls[iccsze++] = nbr;
                 mask[nbr-1]  = 0;
            }
       }
   }
   while ((lvsize = iccsze-lvlend) > 0);

////////////////////////////////////////////////////////////
//...Reset MASK to one for the nodes in the level structure;
   for (xls[nlvl] = lvlend+1, i = 1; i <= iccsze; i++)
       mask[(node = ls[i-1])-1] = 1;
}

///////////////////////////////////
//...finds pseudo-peripheral nodes;
void fnroot(int & root, int * xadj, int * iadj, int * mask, int & nlvl, int * xls, int * ls, int n)
{
   int iccsze, j, jstrt, k, kstop, kstrt, mindeg, nabor, ndeg, node, nunlvl;

   rootls(root, xadj, iadj, mask, nlvl, xls, ls, n);
   if (nlvl == 1 || nlvl == (iccsze = xls[nlvl]-1)) return;

   do {
       mindeg = iccsze;
       root   = ls[(jstrt = xls[nlvl-1])-1];

       if (iccsze > jstrt) {
           for (j = jstrt; j <= iccsze; j++) {
                ndeg  = 0;
                kstrt = xadj[(node = ls[j-1])-1];
                kstop = xadj[ node]-1;

                for (k = kstrt; k <= kstop; k++) 
                if (mask[(nabor = iadj[k-1])-1] > 0 ) ndeg++;

                if (ndeg < mindeg) {
                    root   = node;
                    mindeg = ndeg;
                }
           }
       }

       rootls (root, xadj, iadj, mask, nunlvl, xls, ls, n);
       if (nunlvl <= nlvl) return;
   }
   while ((nlvl = nunlvl) < iccsze);
}

//////////////////////////////////////////////////////////////////
//...computes the degrees of the nodes in the connected component;
void degree (int root, int * xadj, int * iadj, int * mask, int * deg, int & iccsze, int * ls, int n)
{
   int i, ideg, j, jstop, jstrt, lbegin, lvlend, lvsize, nbr, node;

   ls[0]        = root;
   lvlend       = 0;
   iccsze       = 1;
   xadj[root-1] = -xadj[root-1];

   do {
       lbegin = lvlend+1;
       lvlend = iccsze;

       for (i = lbegin; i <= lvlend; i++) {
            jstrt = -xadj[(node = ls[i-1])-1];
//          jstop =  abs(xadj[node])-1;
            jstop = xadj[node];
            if (jstop < 0) jstop = -jstop;
            jstop--;
            ideg  = 0;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 ideg = ideg+1;
                 if (xadj[nbr-1] >= 0) {
                     xadj[nbr-1]  = -xadj[nbr-1];
                     ls[iccsze++] = nbr;
                 }
            }
            deg[node-1] = ideg;
       }
   }
   while ((lvsize = iccsze - lvlend) > 0);

///////////////////////////////////////////////
//...Reset XADJ to its correct sign and return;
   for (i = 1; i <= iccsze; i++) {
        node         = ls[i-1];
        xadj[node-1] = -xadj[node-1];
   }
}


////////////////////////////////////////////////
//...reverses the elements of an integer vector;
#define iSWAP(A, B) { int iSWAP_temp = (A); (A) = (B); (B) = iSWAP_temp; }

void ivec_reverse (int n, int * a)
{
  int  m, i;
  for (m = n/2, i = 1; i <= m; i++)
       iSWAP(a[i-1], a[n-i]);
}

#undef iSWAP

/////////////////////////////////////////////////////////////////////////////
//...numbers a connected component using the reverse Cuthill McKee algorithm;
void rcm(int root, int * xadj, int * iadj, int * mask, int * perm, int & iccsze, int * deg, int n)
{
   int fnbr, i, j, jstop, jstrt, k, l, lbegin, lnbr, lperm, lvlend, nbr, node;

   degree (root, xadj, iadj, mask, deg, iccsze, perm, n);
   mask[root-1] = 0;
   if ( iccsze <= 1) return;

   lvlend = 0;
   lnbr   = 1;

   do {
       lbegin = lvlend+1;
       lvlend = lnbr;

       for (i = lbegin; i <= lvlend; i++) {
            jstrt = xadj[(node = perm[i-1])-1];
            jstop = xadj[ node]-1;
            fnbr  = lnbr+1;

            for (j = jstrt; j <= jstop; j++)
            if  (mask[(nbr = iadj[j-1])-1] != 0) {
                 mask[ nbr-1] = 0;
                 perm[lnbr++] = nbr;
            }
/////////////////////////////////////////////////////////////
//...Sort the neighbors of node in increasing order by degree;
            if (fnbr < lnbr) {
                k = fnbr;
                do {
                    l   = k;
                    nbr = perm[k++];
label40:
                    if (l > fnbr && deg[(lperm = perm[l-1])-1] > deg[nbr-1]) {
                        perm[l--] = lperm;
                        goto label40;
                    }
                    perm[l] = nbr;
                }
                while (k < lnbr);
            }
       }
   }
   while (lnbr > lvlend);

////////////////////////////////////////
//...Reverse the Cuthill-McKee ordering;
   ivec_reverse(iccsze, perm);
}

//////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a general graph;
void genrcm(int n, int * xadj, int * iadj, int * perm, int * xls, int * mask)
{
   int i, iccsze, nlvl, num, root;

   for (         i = 1; i <= n; i++) mask[i-1] = 1;
   for (num = 1, i = 1; i <= n; i++)
   if  (mask[i-1] != 0) {
       fnroot(root = i, xadj, iadj, mask, nlvl, xls,  perm+num-1,  n);
       rcm   (root,     xadj, iadj, mask, perm+num-1, iccsze, xls, n);

       if ((num += iccsze) > n) return;
   }
}

///////////////////////////////////////////////////////////////////
//...finds the reverse Cuthill-Mckee ordering for a given subgraph;
void subrcm(int * xadj, int * iadj, int * mask, int nsubg, int * subg, int * perm, int * xls, int n)
{
   int  i, iccsze, nlvl, node, num;
   for (i = 1; i <= nsubg; i++)
        mask[(node = subg[i-1])-1] = 1;

   for (num = 0, i = 1; i <= nsubg; i++)
   if  (mask[(node = subg[i-1])-1] > 0 ) {
       fnroot(node, xadj, iadj, mask, nlvl, xls,  perm+num,  n);
       rcm   (node, xadj, iadj, mask, perm+num, iccsze, xls, n);

       if ((num += iccsze) >= nsubg ) return;
   }
}

