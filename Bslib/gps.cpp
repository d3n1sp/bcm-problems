//------------------------------------------------------------------------------------------------
// File: gps.cpp
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#if defined(_WINDOWS) && defined(_DEBUG)
   #define _CRTDBG_MAP_ALLOC
   #include <stdlib.h>
#endif

#include "globals.h"

#define   iSWAP(A, B)     { int          iSWAP_temp = (A); (A) = (B); (B) =   iSWAP_temp; }

inline int    min(int    A, int    B) { return A < B ? A : B;};
inline float  min(float  A, float  B) { return A < B ? A : B;};
inline double min(double A, double B) { return A < B ? A : B;};
inline int    max(int    A, int    B) { return A > B ? A : B;};
inline float  max(float  A, float  B) { return A > B ? A : B;};
inline double max(double A, double B) { return A > B ? A : B;};

//inline int    abs(int    A) { return A > 0 ? A : -A;};

//#include "stdafx.h"
//#include "utils.h"
 
void GPSKCB(int N, int * DEGREE, int * RSTART, int * CONNEC, int & AVAIL,  int & NLEFT,  int & STNODE,
                   int & RVNODE, int * WORK,   int & FORWD,  int & BESTBK, int & NNODES, int & DEPTH,
                   int & FWIDTH, int & BWIDTH, int & ERGPS,  int & SPACE);
void GPSKCC(int N, int * DEGREE, int * RSTART, int * CONNEC, int & STNODE, int & AVAIL,  int & NLEFT,
                   int * LIST,   int & ACTIVE, int & DEPTH,  int & WIDTH,  int & ERGPS,  int & SPACE);
void GPSKCD(int N, int * DEGREE, int * RSTART, int * CONNEC, int & STNODE, int & AVAIL,  int & ACTIVE, int & MXDPTH,
                   int * LIST,   int & DEPTH,  int & WIDTH,  int & MAXWID, int & ERGPS,  int & SPACE);
void GPSKCE(int N, int & AVAIL,  int & ACTIVE, int & DEPTH,  int & WRKLEN, int * LVLLST, int * LVLPTR, int * WORK,
                   int & NXTNUM, int & TREE1,  int & TREE2,  int & WIDTH1, int & WIDTH2, int & ONEIS1,
                   int & ERGPS,  int & SPACE);
void GPSKCF(int N, int & ACTIVE, int & DEPTH,  int * LVLLST, int * LVLPTR, int * LVLNUM, int REVERS);
void GPSKCG(int N, int * DEGREE, int * RSTART, int * CONNEC, int & ACTIVE, int & WIDTH1, int & WIDTH2,
                   int * TREE1,  int * TREE2,  int * WORK,   int & WRKLEN, int & DEPTH, int * INC1, int * INC2,
                   int * TOTAL,  int & ONEIS1, int & REVRS1, int & ERGPS,  int & SPACE);
void GPSKCH(int N, int * DEGREE, int * RSTART, int * CONNEC, int * STATUS, int & NREDUC, int * WORK, int & MXCOMP,
                   int * START,  int * SIZE,   int & COMPNS, int & ERGPS,  int & SPACE);
void GPSKCI(int N, int & ACTIVE, int & DEPTH,  int * LSTRUC, int * LVLLST, int * LVLPTR,
                   int * LTOTAL, int & ERGPS,  int & SPACE);
void GPSKCJ(int N, int * DEGREE, int * RSTART, int * CONNEC, int & NCOMPN, int * INVNUM, int & SNODE1, int & SNODE2,
                   int & REVRS1, int & DEPTH,  int * LVLLST, int * LVLPTR, int * LVLNUM, int & ERGPS,  int & SPACE);
void GPSKCK(int N, int * DEGREE, int * RSTART, int * CONNEC, int WRKLEN,   int & NXTNUM, int * WORK,   int & NCOMPN,
                   int & DEPTH,  int * LVLLST, int * LVLPTR, int * LVLNUM, int & ERGPS,  int & SPACE);
void GPSKCL(int N, int * DEGREE, int * RSTART, int * CONNEC, int * INVNUM, int * NEWNUM, int * OLDNUM,
                   int & BANDWD, int & PROFIL, int & ERGPS,  int & SPACE);
void GPSKCM(int N, int * DEGREE, int * RSTART, int * CONNEC, int * INVNUM, int * NEWNUM, int * OLDNUM,
                   int & BANDWD, int & PROFIL, int & ERGPS,  int & SPACE);
void GPSKCN(int N, int * KEY,    int * DATA,   int & ERGPS);
void GPSKCO(int N, int * KEY,                  int & ERGPS);
void GPSKCP(int N, int * INDEX,  int & NVEC,   int * DEGREE, int & ERGPS);
void GPSKCQ(int N, int * INDEX,  int & NVEC,   int * DEGREE, int & ERGPS);
///////////////////////////////////////////////////////////////////////////////////////////////////
//...GIBBS-POOLE-STOCKMEYER AND GIBBS-KING ALGORITHMS
void GPSKCA(int N,        // -- THE DIMENSION OF THE MATRIX
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- DESCRIBE THE STRUCTURE OF THE SPARSE MATRIX.
                          //    DEGREE(I) SPECIFIES THE NUMBER OF NON-ZERO
                          //    OFF-DIAGONAL ENTRIES IN THE  I-TH  ROW OF THE
                          //    SPARSE MATRIX.  THE COLUMN INDICES OF THESE
                          //    ENTRIES ARE GIVEN IN CONSECUTIVE LOCATIONS IN
                          //    CONNEC, STARTING AT LOCATION  RSTART(I).
                          //    IN OTHER WORDS, THE INDICES OF THE NON-ZERO
                          //    OFF-DIAGONAL ELEMENTS OF THE  I-TH  ROW ARE FOUND
                          //    IN:
                          //          CONNEC (RSTART(I)),
                          //          CONNEC (RSTART(I) + 1),
                          //          . . .
                          //          CONNEC (RSTART(I) + DEGREE(I) - 1)
                          //
                          //  DIMENSIONS:
                          //           RSTART IS DIMENSION  N  (OR LONGER).
                          //           DEGREE IS DIMENSION  N  (OR LONGER).
                          //           CONNEC IS DIMENSION ROUGHLY THE NUMBER OF NON-
                          //           ZERO ENTRIES IN THE MATRIX.
            int OPTPRO,   // -- .TRUE. IF REDUCING THE PROFILE OF THE MATRIX
                          //           IS MORE IMPORTANT THAN REDUCING THE
                          //           BANDWIDTH
                          //    .FALSE. IF BANDWIDTH REDUCTION IS MOST IMPORTANT
            int WRKLEN,   // -- THE  ACTUAL  LENGTH OF THE VECTOR  WORK  AS SUPPLIED
                          //    BY THE USER.  SEE THE DISCUSSION OF THE WORKSPACE
                          //    'WORK'  BELOW FOR TYPICAL STORAGE REQUIREMENTS.
                          //    THE VALUE OF  WRKLEN  WILL BE USED TO ENSURE THAT
                          //    THE ROUTINE WILL NOT USE MORE STORAGE THAN IS
                          //    AVAILABLE.  IF NOT ENOUGH SPACE IS GIVEN IN  WORK
                          //    TO PERMIT A SOLUTION TO BE FOUND, THE  ERGPS  FLAG
                          //    WILL BE SET AND FURTHER COMPUTATION STOPPED.
            int * PERMUT, // -- ON INPUT, AN ALTERNATIVE REORDERING FOR THE
                          //    ROWS AND COLUMNS OF THE MATRIX.  PERMUT(I) GIVES
                          //    THE POSITION IN WHICH ROW AND COLUMN  I  SHOULD
                          //    BE PLACED TO REDUCE THE BANDWIDTH OR THE PROFILE.
                          //    IF THE USER HAS NO ALTERNATIVE TO THE NATURAL
                          //    ORDERING IMPLICIT IN  DEGREE,  RSTART  AND  CONNEC,
                          //    HE SHOULD INITIALIZE  PERMUT  TO BE THE IDENTITY
                          //    PERMUTATION  PERMUT(I) = I .
                          //
                          //    ON OUTPUT,  PERMUT  WILL CONTAIN THE PERMUTATION
                          //    FOR REORDERING THE ROWS AND COLUMNS WHICH REDUCES
                          //    THE BANDWIDTH AND/OR PROFILE.  THE RESULT WILL BE
                          //    THE REORDERING FOUND BY 'GPSKCA' OR THE REORDERING
                          //    GIVEN BY THE USER IN 'PERMUT', WHICHEVER DOES THE
                          //    JOB BETTER.
            int * WORK,   // -- A TEMPORARY STORAGE VECTOR, OF LENGTH SOMEWHAT
                          //    GREATER THAN  3N.  THE SPACE BEYOND  3N  REQUIRED
                          //    IS PROBLEM-DEPENDENT.  ANY PROBLEM CAN BE SOLVED
                          //    IN  6N+3  LOCATIONS.
                          //    MOST PROBLEMS CAN BE REORDERED WITH  4N  LOCATIONS
                          //    IN 'WORK'.  IF SPACE IS NOT A CONSTRAINT, PROVIDE
                          //    6N+3  LOCATIONS IN 'WORK'.  OTHERWISE, PROVIDE AS
                          //    MUCH MORE THAN  3N  AS IS CONVENIENT AND CHECK THE
                          //    ERGPS FLAG AND SPACE REQUIRED PARAMETERS (SEE BELOW)
                          //
                          //    ON OUTPUT, THE 1ST  N  LOCATIONS OF WORK WILL BE
                          //    A LISTING OF THE ORIGINAL ROW AND COLUMN INDICES AS
                          //    THEY APPEAR IN THE COMPUTED REORDERING.
                          //    LOCATIONS  N+1, ... , 2N  OF  WORK  WILL CONTAIN
                          //    THE NEW POSITIONS FOR THE EQUATIONS IN THE ORDER
                          //    FOUND BY GPSKCA.  THUS, THE TWO VECTORS ARE INVERSE
                          //    PERMUTATIONS OF EACH OTHER.  IF THE ORDERING
                          //    FOUND BY THIS ALGORITHM IS BETTER THAN THE USER-
                          //    SUPPLIED ORDER, THE SECOND PERMUTATION VECTOR IS
                          //    IDENTICAL TO THE RESULT RETURNED IN  'PERMUT'.
            int & BANDWD, // -- THE BANDWIDTH OF THE MATRIX WHEN ROWS AND COLUMNS
                          //    ARE REORDERED IN THE ORDERING RETURNED IN  PERMUT.
            int & PROFIL, // -- THE PROFILE OF THE MATRIX WHEN ROWS AND COLUMNS ARE
                          //    REORDERED IN THE ORDERING RETURNED IN  PERMUT.
            int & ERGPS,  // -- WILL BE EQUAL TO ZERO IF A NEW NUMBERING COULD BE
                          //    FOUND IN THE SPACE PROVIDED.  OTHERWISE,  ERGPS
                          //    WILL BE SET TO A POSITIVE ERGPS CODE (SEE TABLE
                          //    GIVEN BELOW).  IF THE REORDERING ALGORITHM HAS BEEN
                          //    STOPPED BY LACK OF WORKSPACE, THE SPACE PARAMETER
                          //    WILL BE SET TO THE NUMBER OF ADDITIONAL LOCATIONS
                          //    REQUIRED TO COMPLETE AT LEAST THE NEXT PHASE OF
                          //    THE ALGORITHM.
                          //
                          //    WHENEVER A NON-ZERO VALUE FOR  ERGPS  IS GIVEN
                          //    PERMUT  WILL RETAIN THE VALUES PROVIDED BY THE USER
                          //    AND THE SCALARS  BANDWD  AND  PROFIL  WILL BE SET TO
                          //    OUTRAGEOUS VALUES.  IT IS THE USER'S RESPONSIBILITY
                          //    TO CHECK THE STATUS OF  ERGPS.
            int & SPACE)  // -- WILL INDICATE EITHER HOW MUCH SPACE THE REORDERING
                          //    ACTUALLY REQUIRED OR HOW MUCH SPACE WILL BE
                          //    REQUIRED TO COMPLETE THE NEXT PHASE OF THE
                          //    REORDERING ALGORITHM.  THE POSSIBLE OUTCOMES ARE ..
                          //
                          //         ERGPS = 0          SPACE IS THE MINIMAL VALUE FOR
                          //                            WRKLEN  REQUIRED TO REORDER
                          //                            THIS MATRIX AGAIN.
                          //
                          //         ERGPS <> 0         SPACE IS THE MINIMUM NUMBER
                          //         DUE TO LACK OF     OF EXTRA WORKSPACE REQUIRED
                          //         WORKSPACE          TO CONTINUE THE REORDERING
                          //                            ALGORITHM ON THIS MATRIX.
                          //
                          //         ERGPS <> 0         SPACE = -1
                          //         DUE TO ERGPS
                          //         IN DATA STRUCTURES
{
//INTEGER     DEGREE(N), CONNEC(1), PERMUT(N), WORK(WRKLEN)
  int I, INC1, INC2, AVAIL, NXTNUM, LOWDG, STNODE, NLEFT, TREE1, TREE2, DEPTH, EMPTY, STOTAL,
         REQD, CSPACE, LVLLST, LVLPTR, ACTIVE, RVNODE, WIDTH1, WIDTH2, MXDG, REVRS1, ONEIS1;
///////////////////////////////////////////////////////////////////
//=================================================================
//    << NUMBER ANY DEGREE ZERO NODES >>
//
//    WHILE << SOME NODES YET UNNUMBERED >> DO
//         << FIND A PSEUDO-DIAMETER OF THE MATRIX GRAPH >>
//         << CONVERT FORM OF LEVEL TREES >>
//         << COMBINE LEVEL TREES INTO ONE LEVEL STRUCTURE >>
//         << CONVERT FORM OF LEVEL STRUCTURE >>
//         IF OPTPRO THEN
//             << RENUMBER BY KING ALGORITHM >>
//         ELSE
//             << RENUMBER BY REVERSE CUTHILL-MCKEE ALGORITHM >>
//==================================================================
//////////////////////////////////////////////////////////////////////////////////////////////
// ... INITIALIZE COUNTERS, THEN NUMBER ANY NODES OF DEGREE  0.
//     THE LIST OF NODES, BY NEW NUMBER, WILL BE BUILT IN PLACE AT THE FRONT OF THE WORK AREA;
      NXTNUM = 1; ERGPS = 0; SPACE = 2*N;
      MXDG   = 0;
      for (I = 1; I <= N; I++) {
          if (DEGREE[I-1] < 0) {
              ERGPS = 1; SPACE = -1; goto err;
          }
          else if (DEGREE[I-1] ==   0) WORK[NXTNUM++-1] = I;
          else if (DEGREE[I-1] > MXDG) MXDG = DEGREE[I-1];
      }

/////////////////////////////
//... WHILE  NXTNUM <= N  DO;
      while (NXTNUM <= N) {

////////////////////////////////////////////////
//... FIND AN UNNUMBERED NODE OF MINIMAL DEGREE;
          LOWDG  = MXDG+1;
          STNODE = 0;
          for (I = 1; I <= N; I++)
          if (0 < DEGREE[I-1] && DEGREE[I-1] < LOWDG) LOWDG = DEGREE[(STNODE = I)-1];

          if (STNODE == 0) {
              ERGPS = 2; SPACE = -1; goto err;
          }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//... SET UP POINTERS FOR THREE LISTS IN WORK AREA, THEN LOOK FOR PSEUDO-DIAMETER, BEGINNING WITH STNODE;
          AVAIL = (WRKLEN - NXTNUM + 1) / 3;
          NLEFT = N - NXTNUM + 1;
          SPACE = max(SPACE, NXTNUM + 3*N - 1);
          if (AVAIL < N) {
              ERGPS = 101; SPACE = -1; goto err;
          }
          GPSKCB(N, DEGREE, RSTART, CONNEC, AVAIL, NLEFT, STNODE, RVNODE, WORK+NXTNUM-1, TREE1, TREE2,
                    ACTIVE, DEPTH, WIDTH1, WIDTH2, ERGPS, SPACE);
          if (ERGPS != 0) goto err;

          SPACE = max(SPACE, NXTNUM + 3*(ACTIVE+DEPTH+1) - 1);

/////////////////////////////////////////////////////////////
//... DYNAMIC SPACE CHECK FOR MOST OF REMAINDER OF ALGORITHM;
          REQD  = max(NXTNUM + 2*N + 3*DEPTH - 1, 3*N + 2*DEPTH + 1);
          SPACE = max(SPACE, REQD);
          if (WRKLEN < REQD) {
              ERGPS = 102; SPACE = -1; goto err;
          }

//////////////////////////////////////////////////////////////////////
//... OUTPUT FROM GPSKCB IS A PAIR OF LEVEL TREES, IN THE FORM OF LISTS OF NODES BY LEVEL.
//    CONVERT THIS TO TWO LISTS OF LEVEL NUMBER BY NODE. AT THE SAME TIME PACK
//    STORAGE SO THAT ONE OF THE LEVEL TREE VECTORS IS AT THE BACK END OF THE WORK AREA.
          LVLPTR = NXTNUM + AVAIL - DEPTH;
          GPSKCE (N, AVAIL, ACTIVE, DEPTH, WRKLEN, WORK+NXTNUM-1, WORK+LVLPTR-1, WORK, NXTNUM,
                     TREE1, TREE2, WIDTH1, WIDTH2, ONEIS1, ERGPS, SPACE);
          if (ERGPS != 0) goto err;
          if (TREE1 != WRKLEN-N+1 || TREE2 != NXTNUM) {
              ERGPS = 3; SPACE = -1; goto err;
          }

//////////////////////////////////////////////////////////////////////
//... COMBINE THE TWO LEVEL TREES INTO A MORE GENERAL LEVEL STRUCTURE;
          AVAIL  = WRKLEN - NXTNUM + 1 - 2*N - 3*DEPTH;
          STOTAL = N + NXTNUM;
          EMPTY  = STOTAL + DEPTH;
          INC1   = TREE1 - DEPTH;
          INC2   = INC1 - DEPTH;
          GPSKCG (N, DEGREE, RSTART, CONNEC, ACTIVE, WIDTH1, WIDTH2, WORK+TREE1-1, WORK+TREE2-1, WORK+EMPTY-1,
                     AVAIL,  DEPTH,  WORK+INC1-1, WORK+INC2-1, WORK+STOTAL-1, ONEIS1, REVRS1, ERGPS, CSPACE);
          if (ERGPS != 0) goto err;
          SPACE = max(SPACE, NXTNUM + CSPACE - 1);

/////////////////////////////////////////////////////////////////////////////////////
//... COMBINED LEVEL STRUCTURE IS REPRESENTED BY GPSKCG AS A VECTOR OF LEVEL NUMBERS.
//    FOR RENUMBERING PHASE, CONVERT THIS ALSO TO THE INVERSE PERMUTATION;
          LVLPTR = TREE1 - (DEPTH + 1);
          LVLLST = LVLPTR - ACTIVE;
          if (STOTAL + DEPTH > LVLPTR) {
              ERGPS = 4; SPACE = -1; goto err;
          }
          GPSKCI (N, ACTIVE, DEPTH, WORK+TREE1-1, WORK+LVLLST-1, WORK+LVLPTR-1, WORK+STOTAL-1, ERGPS, SPACE);
          if (ERGPS != 0) goto err;

////////////////////////////////////////////////////////////////////////////////////////
//... NOW RENUMBER ALL MEMBERS OF THIS COMPONENT USING EITHER A REVERSE CUTHILL-MCKEE OR
//    A KING STRATEGY, AS PROFILE OR BANDWIDTH REDUCTION IS MORE IMPORTANT.
          if (OPTPRO) {
              GPSKCK (N, DEGREE, RSTART, CONNEC, LVLLST-1, NXTNUM, WORK, ACTIVE, DEPTH,
                         WORK+LVLLST-1, WORK+LVLPTR-1, WORK+TREE1-1, ERGPS, SPACE);
              if (ERGPS != 0) goto err;
          }
          else {
              GPSKCJ (N, DEGREE, RSTART, CONNEC, ACTIVE, WORK+NXTNUM-1, STNODE, RVNODE,
                         REVRS1, DEPTH,  WORK+LVLLST-1, WORK+LVLPTR-1, WORK+TREE1-1, ERGPS, SPACE);
              if (ERGPS != 0) goto err;
              NXTNUM += ACTIVE;
          }
      }

/////////////////////////////////////////////////////////////////////////////////
//... CHECK WHETHER INITIAL NUMBERING OR FINAL NUMBERING PROVIDES BETTER RESULTS;
      if (WRKLEN < 2*N) {
          ERGPS = 10; SPACE = 2*N-WRKLEN; goto err;
      }
      if (OPTPRO) GPSKCM (N, DEGREE, RSTART, CONNEC, WORK, WORK+N, PERMUT, BANDWD,
                             PROFIL, ERGPS, SPACE);
      else        GPSKCL (N, DEGREE, RSTART, CONNEC, WORK, WORK+N, PERMUT, BANDWD,
                             PROFIL, ERGPS, SPACE);
      return;
err:
      for (I = 1; I  <= N; I++) if (DEGREE[I-1] < 0) DEGREE[I-1] = -DEGREE[I-1];
      BANDWD = PROFIL = -1;

      return;
}

/////////////////////////////////////////////////////////////////////////////////////
//...FIND A PSEUDO-DIAMETER OF THE MATRIX GRAPH
void GPSKCB(int N,        // -- THE DIMENSION OF THE MATRIX.
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- DESCRIBE THE MATRIX STRUCTURE.
            int & AVAIL,
            int & NLEFT,
            int & STNODE, // -- IS INITIALLY THE NUMBER OF A NODE TO BE USED TO START THE PROCESS,
                          //    TO BE THE ROOT OF THE FIRST TREE. ON OUTPUT, STNODE IS THE END
                          //    OF THE PSEUDO-DIAMETER WHOSE LEVEL TREE IS NARROWEST.
            int & RVNODE, // -- WILL BE THE OTHER END OF THE PSEUDO-DIAMETER.
            int * WORK,   // -- WORKING SPACE, OF LENGTH 3*AVAIL, USED TO STORE THREE LEVEL TREES.
            int & FORWD,
            int & BESTBK,
            int & NNODES, // -- WILL BE THE NUMBER OF NODES IN THIS CONNECTED COMPONNENT OF THE MATRIX GRAPH,
                          // I.E., THE LENGTH OF THE LEVEL TREES.
            int & DEPTH,  // -- THE DEPTH OF THE LEVEL TREES BEING RETURNED, I.E., THE LENGTH OF THE PSEUDO-DIAMETER.
            int & FWIDTH,
            int & BWIDTH,
            int & ERGPS,
            int & SPACE)
{
/////////////////////////////////////////////////////////////////////////
//==================================================================
//         << BUILD A LEVEL TREE FROM STNODE >>
//         REPEAT
//             << BUILD A LEVEL TREE FROM EACH NODE 'BKNODE' IN THE
//                DEEPEST LEVEL OF  STNODE'S TREE >>
//             << REPLACE 'STNODE' WITH 'BKNODE' IF A DEEPER AND
//                NARROWER TREE WAS FOUND. >>
//         UNTIL
//             << NO FURTHER IMPROVEMENT MADE >>
//==================================================================
// STRUCTURE OF WORKSPACE ...
// ---------------------------------------------------------------
// : NUMBERED :  TLIST1  PTR1  :  TLIST2  PTR2  :  TLIST3  PTR3  :
// ---------------------------------------------------------------
// TLISTI IS A LIST OF NODES OF LENGTH  'ACTIVE'
// PTRI   IS A LIST OF POINTERS INTO TLISTI, OF LENGTH  'DEPTH+1'
//==================================================================
//    INTEGER     WORK(AVAIL,3)
      int BACKWD, MXDPTH, WIDTH, FDEPTH, LSTLVL, NLAST, I, BKNODE, LSTLVI, IMPROV, 
          id_FORWD  = -1, 
          id_BACKWD = AVAIL+id_FORWD, 
          id_BESTBK = AVAIL+id_BACKWD;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//... BUILD INITIAL LEVEL TREE FROM 'STNODE'. FIND OUT HOW MANY NODES LIE IN THE CURRENT CONNECTED COMPONENT;
      FORWD  = 1;
      BACKWD = 2;
      BESTBK = 3;
      GPSKCC (N, DEGREE, RSTART, CONNEC, STNODE, AVAIL, NLEFT,
                 WORK,   NNODES, DEPTH, WIDTH, ERGPS, SPACE);
      if (ERGPS != 0) return;

////////////////////////////////////////////
//...REPEAT UNTIL NO DEEPER TREES ARE FOUND;
      MXDPTH = AVAIL - NNODES - 1;
      do {
          FWIDTH = WIDTH;
          FDEPTH = DEPTH;
          LSTLVL = AVAIL - DEPTH + 1;
          NLAST  = WORK[LSTLVL-1+id_FORWD] - WORK[LSTLVL+id_FORWD];
          LSTLVL = WORK[LSTLVL  +id_FORWD];
          BWIDTH = N+1;

//////////////////////////////////////////////////////////////////////////////////
//... SORT THE DEEPEST LEVEL OF 'FORWD' TREE INTO INCREASING ORDER OF NODE DEGREE;
          GPSKCQ (NLAST, WORK+LSTLVL+id_FORWD, N, DEGREE, ERGPS);
          if (ERGPS != 0) goto err;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//... BUILD LEVEL TREE FROM NODES IN 'LSTLVL' UNTIL A DEEPER AND NARROWER TREE IS FOUND OR THE LIST IS EXHAUSTED;
          for (IMPROV = 0, I = 1; ! IMPROV && I <= NLAST; I++) {
               LSTLVI = LSTLVL + I - 1;
               BKNODE = WORK[LSTLVI+id_FORWD];
               GPSKCD (N, DEGREE, RSTART, CONNEC, BKNODE, AVAIL, NNODES, MXDPTH,
                          WORK+1+id_BACKWD, DEPTH, WIDTH, BWIDTH, ERGPS, SPACE);
               if (ERGPS != 0) return;

//////////////////////////////////////////////////////////////////////////////
//... NEW DEEPER TREE ... MAKE IT NEW 'FORWD' TREE AND BREAK OUT OF 'DO' LOOP;
              if (DEPTH > FDEPTH) {
                  iSWAP(FORWD, BACKWD); iSWAP(id_FORWD, id_BACKWD);
                  IMPROV = 1;
                  STNODE = BKNODE;
              }
              else
///////////////////////////////////
//... ELSE CHECK FOR NARROWER TREE;
              if (WIDTH < BWIDTH ) {
                  iSWAP(BESTBK, BACKWD); iSWAP(id_BESTBK, id_BACKWD);
                  BWIDTH = WIDTH;
                  RVNODE = BKNODE;
              }
          }
      }
      while (IMPROV);

/////////////////////
//... END OF PROGRAM;
      DEPTH = FDEPTH; return;
err:
      ERGPS = 11; SPACE = -1; return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...BUILD THE LEVEL TREE ROOTED AT 'STNODE' IN THE SPACE PROVIDED IN LIST.
//   CHECK FOR OVERRUN OF SPACE ALLOCATION.
void GPSKCC(int N,        // -- THE DIMENSION OF THE MATRIX.
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- DESCRIBE THE MATRIX STRUCTURE.
            int & STNODE, // -- THE ROOT OF THE LEVEL TREE.
            int & AVAIL,  // -- THE LENGTH OF THE WORKING SPACE AVAILABLE.
            int & NLEFT,  // -- THE NUMBER OF NODES YET TO BE NUMBERED.
            int * LIST,   // -- THE WORKING SPACE.
            int & ACTIVE, // -- THE NUMBER OF NODES IN THE COMPONENT.
            int & DEPTH,  // -- THE DEPTH OF THE LEVEL TREE ROOTED AT  STNODE.
            int & WIDTH,  // -- THE WIDTH OF THE LEVEL TREE ROOTED AT  STNODE.
            int & ERGPS,  // -- ZERO UNLESS STORAGE WAS INSUFFICIENT.
            int & SPACE)  //
{
//    INTEGER     LIST(AVAIL)
      int LSTART, NLEVEL, FRONT, J, NEWNOD, PTR, CDGREE, LFRONT, LISTJ = N;

////////////////////////////////////////////////////////////////////////////////
//... BUILD THE LEVEL TREE USING  LIST  AS A QUEUE AND LEAVING THE NODES IN PLACE.
//    THIS GENERATES THE NODES ORDERED BY LEVEL PUT POINTERS TO THE BEGINNING OF EACH LEVEL,
//    BUILDING FROM THE BACK OF THE WORK AREA.
      ACTIVE = 1;
      DEPTH  = 0;
      WIDTH  = 0;
      ERGPS  = 0;
      LSTART = 1;
      FRONT  = 1;
      LIST  [ACTIVE-1] =  STNODE;
      DEGREE[STNODE-1] = -DEGREE[STNODE-1];
      LIST  [AVAIL-1]  =  1;
      NLEVEL           =  AVAIL;

//////////////////////////////////////////////////////////////
//... REPEAT UNTIL QUEUE BECOMES EMPTY OR WE RUN OUT OF SPACE;
      do {

///////////////////////////////////////////
//... FIRST NODE OF LEVEL. UPDATE POINTERS;
          if (FRONT >= LSTART ) {
              LSTART = ACTIVE + 1;
              WIDTH  = max(WIDTH, LSTART - LIST[NLEVEL-1]);
              DEPTH++;
              if (--NLEVEL <= ACTIVE ) {
                  SPACE = 3*((NLEFT+1-ACTIVE)*DEPTH/NLEFT+(NLEFT+1-ACTIVE));
                  ERGPS = 110; return;
              }
              LIST[NLEVEL-1] = LSTART;
          }

////////////////////////////////////////////////////////////
//... FIND ALL NEIGHBORS OF CURRENT NODE, ADD THEM TO QUEUE;
          LFRONT =  LIST  [ FRONT-1];
          PTR    =  RSTART[LFRONT-1];
          CDGREE = -DEGREE[LFRONT-1];
          if (CDGREE <= 0) {
              ERGPS = 12; SPACE = -1; return;
          }

          for (J = 1; J <= CDGREE; J++) {
               NEWNOD = CONNEC[PTR++-1];

///////////////////////////////////////////////////
//... ADD TO QUEUE ONLY NODES NOT ALREADY IN QUEUE;
               if (DEGREE[NEWNOD-1] > 0) {
                   DEGREE[NEWNOD-1] = -DEGREE[NEWNOD-1];
                   if (NLEVEL <= ++ACTIVE) {
                       SPACE = 3*((NLEFT+1-ACTIVE)*DEPTH/NLEFT+(NLEFT+1-ACTIVE));
                       ERGPS = 110; return;
                   }
                   if (ACTIVE > NLEFT) {
                       ERGPS = 12; SPACE = -1; return;
                   }
                   LIST[ACTIVE-1] = NEWNOD;
               }
          }
      }
      while (++FRONT <= ACTIVE);

////////////////////////////////////////////////
//... YES, THE TREE IS BUILT. UNDO OUR MARKINGS;
      for (J = 1; J <= ACTIVE; J++) {
          LISTJ           =  LIST[J-1];
          DEGREE[LISTJ-1] = -DEGREE[LISTJ-1];
      }
      return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...BUILD THE LEVEL TREE ROOTED AT 'STNODE' IN THE SPACE PROVIDED IN LIST.  OVERFLOW CHECK NEEDED
//   ONLY ON DEPTH OF TREE. BUILD THE LEVEL TREE TO COMPLETION ONLY IF THE WIDTH OF ALL LEVELS IS
//   SMALLER THAN 'MAXWID'.  IF A WIDER LEVEL IS FOUND TERMINATE THE CONSTRUCTION.
void GPSKCD(int N,        //
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- DESCRIBE THE MATRIX STRUCTURE.
            int & STNODE, // -- THE ROOT OF THE LEVEL TREE.
            int & AVAIL,  // -- THE LENGTH OF THE WORKING SPACE AVAILABLE.
            int & ACTIVE, // -- THE NUMBER OF NODES IN THE COMPONENT.
            int & MXDPTH, // -- MAXIMUM DEPTH OF LEVEL TREE POSSIBLE IN ALLOTTED WORKING SPACE.
            int * LIST,   // -- THE WORKING SPACE.
            int & DEPTH,  // -- THE DEPTH OF THE LEVEL TREE ROOTED AT  STNODE.
            int & WIDTH,  // -- THE WIDTH OF THE LEVEL TREE ROOTED AT  STNODE.
            int & MAXWID, // -- LIMIT ON WIDTH OF THE TREE.  TREE WILL NOT BE USED IF WIDTH OF ANY LEVEL IS AS GREAT AS
                          //    MAXWID, SO CONSTRUCTION OF TREE NEED NOT CONTINUE IF ANY LEVEL THAT WIDE IS FOUND.
            int & ERGPS,  // -- ZERO UNLESS STORAGE WAS INSUFFICIENT.
            int & SPACE)  //
{
//    INTEGER     LIST(AVAIL)
      int LSTART, NLEVEL, FRONT, J, NEWNOD, PTR, BACK, SPTR, FPTR, LFRONT, LISTJ;

/////////////////////////////////////////////////////////////////////////////////////////////
//... BUILD THE LEVEL TREE USING  LIST  AS A QUEUE AND LEAVING THE NODES IN PLACE.
//    THIS GENERATES THE NODES ORDERED BY LEVEL PUT POINTERS TO THE BEGINNING OF EACH LEVEL,
//    BUILDING FROM THE BACK OF THE WORK AREA.
      BACK   = 1;
      DEPTH  = 0;
      WIDTH  = 0;
      ERGPS  = 0;
      LSTART = 1;
      FRONT  = 1;
      LIST  [BACK  -1] =  STNODE;
      DEGREE[STNODE-1] = -DEGREE[STNODE-1];
      LIST  [AVAIL -1] =  1;
      NLEVEL           =  AVAIL;

//////////////////////////////////////////////////////////////
//... REPEAT UNTIL QUEUE BECOMES EMPTY OR WE RUN OUT OF SPACE.
      do {

///////////////////////////////////////////
//... FIRST NODE OF LEVEL. UPDATE POINTERS.
          if (FRONT >= LSTART) {
              LSTART = BACK + 1;
              WIDTH  = max(WIDTH, LSTART - LIST[NLEVEL-1]);

//////////////////////////////////////////////////////////////
//... ABORT GENERATION OF TREE BECAUSE IT IS ALREADY TOO WIDE;
              if (WIDTH >= MAXWID) {
                  WIDTH  = N+1; DEPTH = 0; goto err;
              }
              NLEVEL--;
              if (++DEPTH > MXDPTH) {
                  SPACE = 3*((ACTIVE+1-BACK)*DEPTH/ACTIVE+(ACTIVE+1-BACK)); ERGPS = 111;
                  return;
              }
              LIST[NLEVEL-1] = LSTART;
          }

////////////////////////////////////////////////////////////
//... FIND ALL NEIGHBORS OF CURRENT NODE, ADD THEM TO QUEUE;
          LFRONT   = LIST  [ FRONT-1];
          SPTR     = RSTART[LFRONT-1];
          FPTR     = SPTR - DEGREE[LFRONT-1] - 1;
///////////////////////////////////////////////////
//... ADD TO QUEUE ONLY NODES NOT ALREADY IN QUEUE;
          for (PTR = SPTR; PTR <= FPTR; PTR++) {
              if (DEGREE[(NEWNOD   =  CONNEC[PTR-1])-1] > 0 ) {
                  DEGREE[NEWNOD-1] = -DEGREE[NEWNOD -1];
                  LIST[++BACK-1]   =  NEWNOD;
              }
          }
      }
      while (++FRONT <= BACK);

////////////////////////////////////////////////
//...YES, THE TREE IS BUILT.  UNDO OUR MARKINGS;
     if (BACK != ACTIVE) {
         ERGPS = 13; SPACE = -1; return;
      }
err:
      for (J = 1; J <= BACK; J++) {
          LISTJ           =  LIST[J-1];
          DEGREE[LISTJ-1] = -DEGREE[LISTJ-1];
      }
      return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...TRANSITION BETWEEN ALGORITHM I AND ALGORITHM II OF THE GIBBS-POOLE-STOCKMEYER PAPER
void GPSKCE(int N,
            int & AVAIL,
            int & ACTIVE,
            int & DEPTH,
            int & WRKLEN,
            int * LVLLST,
            int * LVLPTR,
            int * WORK,
            int & NXTNUM,
            int & TREE1,
            int & TREE2,
            int & WIDTH1,
            int & WIDTH2,
            int & ONEIS1,
            int & ERGPS,
            int & SPACE)
{
////////////////////////////////////////////////////////////////////
//==================================================================
//...IN THIS IMPLEMENTATION ALGORITHM I REPRESENTS LEVEL TREES AS
//   LISTS OF NODES ORDERED BY LEVEL.  ALGORITHM II APPEARS TO REQUIRE
//   LEVEL NUMBERS INDEXED BY NODE -- VECTORS FOR EFFICIENCY.
//   THIS SUBROUTINE CHANGES THE LEVEL TREE REPRESENTATION TO THAT
//   REQUIRED BY ALGORITHM II.  NOTE THAT THE FIRST ALGORITHM CAN BE
//   CARRIED OUT WITH THE LEVEL NUMBER VECTOR FORMAT, PROBABLY REQURING
//   MORE COMPUTATION TIME, BUT PERHAPS LESS STORAGE.
//
//   INPUT:  TWO LEVEL TREES, AS LEVEL LISTS AND LEVEL POINTERS,
//           FOUND IN TWO OF THE THREE COLUMNS OF THE ARRAYS 'LVLLST'
//           AND 'LVLPTR'
//
//   OUTPUT: TWO LEVEL TREES, AS VECTORS OF LEVEL NUMBERS,
//           ONE PACKED TO THE FRONT, ONE TO THE REAR OF THE WORKING
//           AREA 'WORK'.  NOTE THAT 'WORK', 'LVLLST' AND 'LVLPTR'
//           SHARE COMMON LOCATIONS.
//================================================================
//... STRUCTURE OF WORKSPACE
//
//    INPUT .. (OUTPUT FROM GPSKCB)
//--------------------------------------------------------------
//: NUMBERED : TLIST1  PTR1  :  TLIST2  PTR2  :  TLIST3  PTR3  :
//--------------------------------------------------------------
//
//    OUTPUT .. (GOES TO COMBIN)
//--------------------------------------------------------------
//: NUMBERED :  TREE2  :           ...               :  TREE1  :
//--------------------------------------------------------------
//==================================================================
//    INTEGER  LVLLST(AVAIL,3), LVLPTR(AVAIL,3), WORK(WRKLEN)
      int   I, BTREE, FTREE, FWIDTH, BWIDTH,
            id_AVAIL1 = AVAIL-1,
            id_AVAIL2 = id_AVAIL1+AVAIL;

///////////////////////////////////////////////////////////////////
//... CHECK THAT WE HAVE ENOUGH ROOM TO DO THE NECESSARY UNPACKING;
      if (3*AVAIL > WRKLEN) {
          ERGPS = 20; SPACE = -1; return;
      }
      if (AVAIL < N) {
          SPACE = 3*(N-AVAIL); ERGPS = 120; return;
      }

/////////////////////////////////////
//... INPUT HAS THREE POSSIBLE CASES:
//    LVLLST(*,1) IS EMPTY;
//    LVLLST(*,2) IS EMPTY;
//    LVLLST(*,3) IS EMPTY;
      FTREE  = TREE1;
      BTREE  = TREE2;
      FWIDTH = WIDTH1;
      BWIDTH = WIDTH2;
      TREE1  = WRKLEN-N+1;
      TREE2  = NXTNUM;

///////////////////////////////////////////////////////////
//... CASE 1: 1ST SLOT IS EMPTY. UNPACK 3 INTO 1, 2 INTO 3;
      if (FTREE != 1 && BTREE != 1) {
          if (FTREE == 2) {
              ONEIS1 = 1;
              WIDTH2 = BWIDTH;
              WIDTH1 = FWIDTH;
          }
          else {
              ONEIS1 = 0;
              WIDTH1 = BWIDTH;
              WIDTH2 = FWIDTH;
          }
          GPSKCF (N, ACTIVE, DEPTH, LVLLST+1+id_AVAIL2, LVLPTR+1+id_AVAIL2, WORK+TREE2-1,   ONEIS1);
          GPSKCF (N, ACTIVE, DEPTH, LVLLST+1+id_AVAIL1, LVLPTR+1+id_AVAIL1, WORK+TREE1-1, ! ONEIS1);
          return;
      }

//////////////////////////////////////////////////////////////////////////////////////////////////////
//... CASE 2: 2ND SLOT IS EMPTY. TO ENABLE COMPLETE REPACKING, MOVE 3 INTO 2, THEN FALL INTO NEXT CASE
      if (FTREE != 2 && BTREE != 2) {
          for (I = 1; I <= ACTIVE; I++) LVLLST[I+id_AVAIL1] = LVLLST[I+id_AVAIL2];
          for (I = 1; I <= DEPTH ; I++) LVLPTR[I+id_AVAIL1] = LVLPTR[I+id_AVAIL2];
      }

//////////////////////////////////////////////////////////////
//... CASE 3:  SLOT 3 IS EMPTY.  MOVE 1 INTO 3, THEN 2 INTO 1;
      if (FTREE != 1) {
          ONEIS1 = 0;
          WIDTH1 = BWIDTH;
          WIDTH2 = FWIDTH;
      }
      else {
          ONEIS1 = 1;
          WIDTH1 = FWIDTH;
          WIDTH2 = BWIDTH;
      }
      GPSKCF (N, ACTIVE, DEPTH, LVLLST,             LVLPTR,             WORK+TREE1-1, ! ONEIS1);
      GPSKCF (N, ACTIVE, DEPTH, LVLLST+1+id_AVAIL1, LVLPTR+1+id_AVAIL1, WORK+TREE2-1,   ONEIS1);
      return;
}

////////////////////////////////////////////////////////////////////////////////
//...CONVERT LEVEL STRUCTURE REPRESENTATION FROM A LIST OF NODES GROUPED BY LEVEL
//   TO A VECTOR GIVING LEVEL NUMBER FOR EACH NODE.
void GPSKCF(int N,        //
            int & ACTIVE, //
            int & DEPTH,  //
            int * LVLLST, //
            int * LVLPTR, // -- LIST OF LISTS.
            int * LVLNUM, // -- OUTPUT VECTOR OF LEVEL NUMBERS.
            int REVERS)   // -- IF == 1, NUMBER LEVEL STRUCTURE FROM BACK END INSTEAD OF FROM FRONT.
{
//    INTEGER     LVLLST(ACTIVE), LVLPTR(DEPTH), LVLNUM(N)
      int I, LEVEL, LSTART, LEND, XLEVEL, PLSTRT;

////////////////////////////////////////////////////////////////////////////////////
//... IF NOT ALL NODES OF GRAPH ARE ACTIVE, MASK OUT THE NODES WHICH ARE NOT ACTIVE;
      if (ACTIVE != N)
      for (I = 1; I <= N; I++) LVLNUM[I-1] = 0;
      for (LEVEL = 1; LEVEL <= DEPTH; LEVEL++) {
          PLSTRT = DEPTH-LEVEL+1;
          XLEVEL = REVERS ? PLSTRT : LEVEL;
          LSTART = LVLPTR[PLSTRT-1];
          LEND   = LVLPTR[PLSTRT-2]-1;
          for (I = LSTART; I <= LEND; I++) LVLNUM[LVLLST[I-1]-1] = XLEVEL;
      }
      return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...COMBINE THE TWO ROOTED LEVEL TREES INTO A SINGLE LEVEL STRUCTURE WHICH MAY HAVE SMALLER WIDTH
//   THAN EITHER OF THE TREES. THE NEW STRUCTURE IS NOT NECESSARILY A ROOTED STRUCTURE.
void GPSKCG(int N,        //
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- GIVE THE DIMENSION AND STRUCTURE OF THE SPARSE SYMMETRIC MATRIX.
            int & ACTIVE, // -- THE NUMBER OF NODES IN THIS CONNECTED COMPONENT OF THE MATRIX GRAPH.
            int & WIDTH1, // -- THE MAXIMUM WIDTH OF A LEVEL IN TREE1.
            int & WIDTH2, // -- THE MAXIMUM WIDTH OF A LEVEL IN TREE2.
            int * TREE1,  // -- ON INPUT, ONE OF THE INPUT LEVEL TREES. ON OUTPUT, THE COMBINED LEVEL STRUCTURE.
            int * TREE2,  // -- THE SECOND INPUT LEVEL TREE.
            int * WORK,   // -- A WORKING AREA OF LENGTH 'WRKLEN'.
            int & WRKLEN, //
            int & DEPTH,  //
            int * INC1,   //
            int * INC2,   //
            int * TOTAL,  // -- VECTORS OF LENGTH 'DEPTH'.
            int & ONEIS1, // -- INDICATES WHETHER TREE1 OR TREE2 REPRESENTS THE FORWARD TREE OR
                          //    THE BACKWARDS TREE OF PHASE 1. USED TO MIMIC ARBITRARY TIE-BREAKING
                          //    PROCEDURE OF ORIGINAL GIBBS-POOLE-STOCKMEYER CODE.
            int & REVRS1, // -- OUTPUT PARAMETER INDICATING WHETHER A BACKWARDS ORDERING WAS USED FOR
                          //    THE LARGEST COMPONENT OF THE REDUCED GRAPH
            int & ERGPS,  // -- NON-ZERO ONLY IF FAILURE OF SPACE ALLOCATION OR DATA STRUCTURE ERGPS FOUND.
            int & SPACE)  // -- MINIMUM SPACE REQUIRED TO RERUN OR COMPLETE PHASE.
{
////////////////////////////////////////////////////////////////////
//==================================================================
//     << REMOVE ALL NODES OF PSEUDO-DIAMETERS >>
//     << FIND CONNECTED COMPONENTS OF REDUCED GRAPH >>
//     << COMBINE LEVEL TREES, COMPONENT BY COMPONENT >>
//==================================================================
//  STRUCTURE OF WORKSPACE ...
//------------------------------------------------------------------
//: NUMBERED : TREE2 : TOTAL : NODES : START : SIZE : INC1 : INC2 :
//------------------------------------------------------------------
//
//--------
//  TREE1 : NUMBERED  IS THE SET OF  NUMBERED NODES (PROBABLY EMPTY)
//--------
//  TREE1 AND TREE1 ARE LEVEL TREES (LENGTH N)
//  TOTAL, INC1 AND INC2  ARE VECTORS OF NODE COUNTS PER LEVEL (LENGTH 'DEPTH')
//  NODES IS THE SET OF NODES IN THE REDUCED GRAPH (THE NODES NOT ON ANY SHORTEST PATH FROM ONE END OF THE
//         PSEUDODIAMETER TO THE OTHER)
//  START, SIZE ARE POINTERS INTO 'NODES', ONE OF EACH FOR EACH CONNECTED COMPONENT OF THE REDUCED GRAPH.
//         THE SIZES OF NODES, START AND SIZE ARE NOT KNOWN APRIORI.
//==================================================================
//    INTEGER     TREE1(N), TREE2(N), INC1(DEPTH), INC2(DEPTH), TOTAL(DEPTH)
      int I, SIZE, AVAIL, CSTOP, START, COMPON, TREE1I, CSTART, MXINC1, MXINC2, COMPNS,
             MXCOMP, OFFDIA, CSIZE, WORKI;

//////////////////////////////////////////////////////////////////////////////////////////////////////
//... FIND ALL SHORTEST PATHS FROM START TO FINISH. REMOVE NODES ON THESE PATHS AND IN OTHER CONNECTED
//    COMPONENTS OF FULL GRAPH FROM FURTHER CONSIDERATION. SIGN OF ENTRIES IN TREE1 IS USED AS A MASK.
      OFFDIA = ACTIVE;
      for (I = 1; I <= DEPTH; I++) TOTAL[I-1] = 0;
      for (I = 1; I <= N;     I++) {
          TREE1I = TREE1[I-1];
          if (TREE1I == TREE2[I-1] && TREE1I != 0) {
              TOTAL[TREE1I-1]++;
              TREE1[I-1] = -TREE1[I-1];
              OFFDIA--;
          }
      }
      if (OFFDIA == 0) {//... DEFAULT WHEN THE REDUCED GRAPH IS EMPTY;
          REVRS1  = 1; SPACE = 2*N; return;
      }
      if (OFFDIA  < 0) {
          ERGPS = 30; SPACE = -1; return;
      }

///////////////////////////////////////////////////////////////////////////////////////////////
//... FIND CONNECTED COMPONENTS OF GRAPH INDUCED BY THE NODES NOT REMOVED.
//   'MXCOMP' IS THE LARGEST NUMBER OF COMPONENTS REPRESENTABLE IN THE WORKING SPACE AVAILABLE.
      AVAIL  = WRKLEN - OFFDIA;
      MXCOMP = AVAIL/2;
      START  = OFFDIA + 1;
      SIZE  = START + MXCOMP;
      if (MXCOMP <= 0) {
          ERGPS = 131; SPACE = 2-AVAIL; return;
      }
      GPSKCH(N, DEGREE, RSTART, CONNEC, TREE1, OFFDIA, WORK, MXCOMP,
                WORK+START-1, WORK+SIZE-1, COMPNS, ERGPS, SPACE);
      if (ERGPS != 0) {
          SPACE = -1; return;
      }

/////////////////////////////////////////////////////////
//... RECORD SPACE ACTUALLY USED (NOT INCLUDING NUMBERED);
      SPACE = 2*N+3*DEPTH+2*COMPNS+OFFDIA;

///////////////////////////////////////////////////////////////////////////////////
//... SORT THE COMPONENT START POINTERS INTO INCREASING ORDER OF SIZE OF COMPONENT;
      if (COMPNS > 1) GPSKCN(COMPNS, WORK+SIZE-1, WORK+START-1, ERGPS);
      if (ERGPS != 0) {
          ERGPS = 32; SPACE = -1; return;
      }

/////////////////////////////////////////////////////////////////////////////////////////////////
//... FOR EACH COMPONENT IN TURN, CHOOSE TO USE THE ORDERING OF THE 'FORWARD' TREE1 OR OF THE
//   'BACKWARD' TREE2 TO NUMBER THE NODES IN THIS COMPONENT. THE NUMBERING IS CHOSEN TO MINIMIZE
//    THE MAXIMUM INCREMENT TO ANY LEVEL.
      for (COMPON = 1; COMPON <= COMPNS; COMPON++) {
           CSTART = WORK[(START+COMPON-1)-1];
           CSIZE  = WORK[(SIZE +COMPON-1)-1];
           CSTOP  = CSTART + CSIZE - 1;
           if (CSIZE < 0 ||  CSIZE > OFFDIA) {
               ERGPS = 31; SPACE = -1; return;
           }
           for (I = 1; I <= DEPTH; I++) INC1[I-1] = INC2[I-1] = 0; MXINC1 = MXINC2 = 0;
           for (I = CSTART; I <= CSTOP; I++) {
                INC1[(-TREE1[(WORKI = WORK[I-1])-1])-1]++;
                INC2[( TREE2[ WORKI-1])-1]++;
           }

/////////////////////////////////////////////////////////////////////////
//... BAROQUE TESTS BELOW DUPLICATE THE GIBBS-POOLE-STOCKMEYER-CRANE PROGRAM,
//*** NOT *** THE PUBLISHED ALGORITHM.
          for (I = 1; I <= DEPTH; I++)
           if (INC1[I-1] != 0 || INC2[I-1] != 0) {
               if (MXINC1 < TOTAL[I-1]+INC1[I-1]) MXINC1 = TOTAL[I-1]+INC1[I-1];
               if (MXINC2 < TOTAL[I-1]+INC2[I-1]) MXINC2 = TOTAL[I-1]+INC2[I-1];
          }

///////////////////////////////////////////////////////////////////////////////////
//... USE ORDERING OF NARROWER TREE UNLESS IT INCREASES WIDTH MORE THAN WIDER TREE.
//    IN CASE OF TIE, USE TREE 2!
          if (MXINC1  > MXINC2 ||
             (MXINC1 == MXINC2 && (WIDTH1  > WIDTH2 ||
                                  (WIDTH1 == WIDTH2 && ONEIS1)))) {
              if  (COMPON == 1) REVRS1 = ONEIS1;
              for (I = CSTART; I <= CSTOP; I++) {
                   WORKI = WORK[I-1];
                   TREE1[WORKI-1] = -TREE2[WORKI-1];
              }
              for (I = 1; I <= DEPTH; I++) TOTAL[I-1] += INC2[I-1];
          }
          else {
              if  (COMPON == 1)  REVRS1 = ! ONEIS1;
              for (I = 1; I <= DEPTH; I++) TOTAL[I-1] += INC1[I-1];
          }
      }
      return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
//...FIND THE CONNECTED COMPONENTS OF THE GRAPH INDUCED BY THE SET OF NODES WITH POSITIVE 'STATUS'.
//   WE SHALL BUILD THE LIST OF CONNECTED COMPONENTS IN 'WORK', WITH A LIST OF POINTERS TO
//   THE BEGINNING NODES OF COMPONENTS LOCATED IN 'START'
void GPSKCH(int N,        // -- DIMENSION OF THE ORIGINAL MATRIX.
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- THE STRUCTURE OF THE ORIGINAL MATRIX.
            int * STATUS, // -- DERIVED FROM A LEVEL TREE. POSITIVE ENTRIES INDICATE
                          //    ACTIVE NODES.  NODES WITH STATUS <= 0 ARE IGNORED.
            int & NREDUC, // -- THE NUMBER OF ACTIVE NODES.
            int * WORK,   // -- WORK SPACE, USED AS A QUEUE TO BUILD CONNECTED COMPONENTS IN PLACE.
            int & MXCOMP, // -- MAXIMUM NUMBER OF COMPONENTS ALLOWED BY CURRENT SPACE ALLOCATION. MUST NOT BE VIOLATED.
            int * START,  // -- POINTER TO BEGINNING OF  I-TH  CONNECTED COMPONENT.
            int * SIZE,   // -- SIZE OF EACH COMPONENT.
            int & COMPNS, // -- NUMBER OF COMPONENTS ACTUALLY FOUND.
            int & ERGPS,  // -- SHOULD BE ZERO ON RETURN UNLESS WE HAVE TOO LITTLE SPACE
                          //    OR WE ENCOUNTER AN ERGPS IN THE DATA STRUCTURE.
            int & SPACE)  // -- MAXIMUM AMOUNT OF WORKSPACE USED / NEEDED.
{
////////////////////////////////////////////////////////////////////
//==================================================================
//     REPEAT
//         << FIND AN UNASSIGNED NODE AND START A NEW COMPONENT >>
//         REPEAT
//             << ADD ALL NEW NEIGHBORS OF FRONT NODE TO QUEUE, >>
//             << REMOVE FRONT NODE.                            >>
//         UNTIL <<QUEUE EMPTY>>
//     UNTIL << ALL NODES ASSIGNED >>
//==================================================================
//    INTEGER   WORK(NREDUC), START(MXCOMP), SIZE(MXCOMP)
      int I, J, FREE, JPTR, NODE, JNODE, FRONT, CDGREE, ROOT;

//////////////////////////////////////////////////////////
//... START OF OUTER REPEAT LOOP, FIND AN UNASSIGNED NODE;
      FREE   = 1;
      COMPNS = 0;
      ROOT   = 1;
      do {
          for (I = ROOT; I <= N; I++) if (STATUS[I-1] > 0) {
               NODE = I; break;
          }
          if (I > N) {
              ERGPS = 34; SPACE = -1; return;
          }

//////////////////////////
//... START NEW COMPONENT;
          COMPNS++;
          ROOT   = NODE   + 1;
          if (COMPNS > MXCOMP) {
              ERGPS = 130; SPACE = NREDUC - FREE + 1; return;
          }
          START [COMPNS-1] =  FREE;
          WORK  [FREE  -1] =  NODE;
          STATUS[NODE  -1] = -STATUS[NODE-1];
          FRONT = FREE++;

/////////////////////////////////////////////
//... INNER REPEAT UNTIL QUEUE BECOMES EMPTY;
          do {
              NODE = WORK[FRONT++ -1];
              JPTR = RSTART  [NODE-1];
              CDGREE = DEGREE[NODE-1];
              for (J = 1; J <= CDGREE; J++) {
                   JNODE = CONNEC[JPTR++-1];
                   if (STATUS[JNODE-1] > 0) {
                       STATUS[JNODE-1] = -STATUS[JNODE-1];
                       WORK[FREE++ -1] = JNODE;
                   }
                   if (STATUS[JNODE-1] == 0) {
                       ERGPS = 33; SPACE = -1; return;
                   }
              }
          }
          while (FRONT < FREE);

////////////////////////////////////////////////////////////////////////////////////////////////////
//... END OF INNER REPEAT. COMPUTE SIZE OF COMPONENT AND SEE IF THERE ARE MORE NODES TO BE ASSIGNED
          SIZE[COMPNS-1] = FREE-START[COMPNS-1];
      }
      while (FREE <= NREDUC);

//////////////////
//... ERROR FOUND.
      if (FREE != NREDUC+1) {
          ERGPS = 35; SPACE = -1;
      }
      return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...TRANSITIONAL SUBROUTINE, ALGORITHM II TO IIIA OR IIIB. CONVERT LEVEL STRUCTURE GIVEN AS VECTOR
//   OF LEVEL NUMBERS FOR NODES TO STRUCTURE AS LIST OF NODES BY LEVEL
void GPSKCI(int N,        //
            int & ACTIVE, //
            int & DEPTH,  // -- PROBLEM SIZES.
            int * LSTRUC, // -- INPUT LEVEL STRUCTURE.
            int * LVLLST, //
            int * LVLPTR, // -- OUTPUT LEVEL STRUCTURE.
            int * LTOTAL, // -- NUMBER OF NODES AT EACH LEVEL (PRECOMPUTED).
            int & ERGPS,  //
            int & SPACE)  //
{
////////////////////////////////////////////////////////////////////////
//==================================================================
//     STRUCTURE OF WORKSPACE ..
//
// INPUT (FROM COMBIN) ..
//------------------------------------------------------------------
//:  NUMBERED  :  ..(N)..  :  TOTAL  :         ...        :  TREE  :
//------------------------------------------------------------------
//
// OUTPUT (TO GPSKCJ OR GPSKCK) ..
//------------------------------------------------------------------
//:  NUMBERED  :       ...             :  TLIST  :  TPTR  :  TREE  :
//------------------------------------------------------------------
// HERE, NUMBERED IS THE SET OF NODES IN NUMBERED COMPONENTS TOTAL IS A VECTOR OF
// LENGTH 'DEPTH' GIVING THE NUMBER OF NODES IN EACH LEVEL OF THE 'TREE'.
// TLIST, TPTR ARE LISTS OF NODES OF THE TREE, ARRANGED BY LEVEL.
// TLIST IS OF LENGTH 'ACTIVE', TPTR 'DEPTH+1'.
//=================================================================
//    INTEGER LSTRUC(N), LVLLST(ACTIVE), LVLPTR(1), LTOTAL(DEPTH)
      int  I, ACOUNT, START, LEVEL; ACTIVE = ACTIVE;

////////////////////////////////////////////////////////////
//... ESTABLISH STARTING AND ENDING POINTERS FOR EACH LEVEL;
      for (START = 1, I = 1; I <= DEPTH; I++) {
           LVLPTR[I-1] =  START;
           LTOTAL[I-1] = (START += LTOTAL[I-1]);
      }
      LVLPTR[DEPTH] = START;
      for (ACOUNT = 0, I = 1; I <= N; I++) if (LSTRUC[I-1] > 0) {
           ERGPS = 40; SPACE = -1; return;
      }
      else
      if (LSTRUC[I-1] < 0) {
          LSTRUC[I-1] = LEVEL       = -LSTRUC[I-1];
          LVLLST[LVLPTR[LEVEL-1]-1] =  I;
          LVLPTR[LEVEL-1]++;
          ACOUNT         ++;
          if (LVLPTR[LEVEL-1] > LTOTAL[LEVEL-1]) {
              ERGPS = 41; SPACE = -1; return;
          }
      }

//////////////////////////////
//... RESET STARTING POINTERS;
      for (LVLPTR[0] = 1, I = 1; I <= DEPTH; I++) LVLPTR[I] = LTOTAL[I-1];
      return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...NUMBER THE NODES IN A GENERALIZED LEVEL STRUCTURE ACCORDING TO A GENERALIZATION OF
//   THE CUTHILL MCKEE STRATEGY.
void GPSKCJ(int N,        // -- DIMENSION OF ORIGINAL PROBLEM.
            int * DEGREE, //
            int * RSTART, //
            int * CONNEC, // -- GIVE STRUCTURE OF SPARSE AND SYMMETRIC MATRIX.
            int & NCOMPN, // -- NUMBER OF NODES IN THIS COMPONENT OF MATRIX GRAPH.
            int * INVNUM, // -- WILL BECOME A LIST OF THE ORIGINAL NODES IN THE ORDER
                          //    WHICH REDUCES THE BANDWIDTH OF THE MATRIX.
            int & SNODE1, //
            int & SNODE2, //
            int & REVRS1, // -- IF == 1, FIRST COMPONENT OF REDUCED GRAPH WAS NUMBERED BACKWARDS.
            int & DEPTH,  //
            int * LVLLST, // -- LIST OF NODES IN LEVEL TREE ORDERED BY LEVEL.
            int * LVLPTR, // -- POSITION OF INITIAL NODE IN EACH LEVEL OF LVLLST.
            int * LVLNUM, // -- LEVEL NUMBER OF EACH NODE IN COMPONENT.
            int & ERGPS,  //
            int & SPACE)  //
{
// NXTNUM -- THE NEXT INDEX TO BE ASSIGNED (1 FOR FIRST COMPONENT)
/////////////////////////////////////////////////////////////////////////
//=======================================================================
// NUMBERING REQUIRES TWO QUEUES, WHICH CAN BE BUILD IN PLACE IN INVNUM.
//
//     << SET QUEUE1 TO BE THE SET CONTAINING ONLY THE START NODE. >>
//     FOR LEVEL = 1 TO DEPTH DO
//         BEGIN
//         LOOP
//             REPEAT
//                 BEGIN
//                 << CNODE <- FRONT OF QUEUE1                        >>
//                 << ADD UNNUMBERED NEIGHBORS OF CNODE TO THE BACK   >>
//                 << OF QUEUE1 OR QUEUE2 (USE QUEUE1 IF NEIGHBOR     >>
//                 << AT SAME LEVEL, QUEUE2 IF AT NEXT LEVEL).  SORT  >>
//                 << THE NEWLY QUEUED NODES INTO INCREASING ORDER OF >>
//                 << DEGREE.  NUMBER CNODE, DELETE IT FROM QUEUE1.   >>
//                 END
//             UNTIL
//                 << QUEUE1 IS EMPTY >>
//         EXIT IF << ALL NODES AT THIS LEVEL NUMBERED >>
//             BEGIN
//             << FIND THE UNNUMBERED NODE OF MINIMAL DEGREE AT THIS >>
//             << LEVEL, RESTART QUEUE1 WITH THIS NODE.              >>
//             END
//         END << LOOP LOOP >>
//         << PROMOTE QUEUE2 TO BE INITIAL QUEUE1 FOR NEXT ITERATION >>
//         << OF  FOR  LOOP.                                         >>
//         END <<FOR LOOP>>
//=======================================================================
//  STRUCTURE OF WORKSPACE ..
//--------------------------------------------------------------
//: NUMBERED :  QUEUE1  :  QUEUE2  : ... : TLIST : TPTR : TREE :
//--------------------------------------------------------------
//  ON COMPLETION, WE HAVE ONLY A NEW, LONGER NUMBERED SET.
//=======================================================================
//    INTEGER     INVNUM(NCOMPN), LVLLST(NCOMPN), LVLPTR(DEPTH), LVLNUM(N)
      int I, BQ1, BQ2, FQ1, INC, CPTR, CNODE, INODE, LEVEL, NLEFT, LSTART, LWIDTH, QUEUE1, QUEUE2,
             CDGREE, XLEVEL, STNODE, ILEVEL, SQ1, SQ2, NSORT, LOWDG, BPTR, LVLLSC, LVLLSB, INVNMI,
             FORWRD, RLEVEL;

///////////////////////////////////////////////////////
//... GIBBS-POOLE-STOCKMEYER HEURISTIC CHOICE OF ORDER;
      if (DEGREE[SNODE1-1] <= DEGREE[SNODE2-1]) {
          FORWRD = REVRS1;
          STNODE = SNODE1;
      }
      else {
          FORWRD = ! REVRS1;
          STNODE = SNODE2;
      }

//////////////////////////////////////////////////////////////////////////////////////////////
//... SET UP INITIAL QUEUES AT FRONT OF 'INVNUM' FOR FORWRD ORDER, AT BACK FOR REVERSED ORDER.
      if (! FORWRD) {
          INC    = -1;
          QUEUE1 =  NCOMPN;
      }
      else {
          INC    = 1;
          QUEUE1 = 1;
      }
      INVNUM[QUEUE1-1] = STNODE;
      RLEVEL           = LVLNUM[STNODE-1] == DEPTH;
      LVLNUM[STNODE-1] = 0;
      FQ1 = QUEUE1;
      BQ1 = QUEUE1 + INC;

/////////////////////////////////
//... NUMBER NODES LEVEL BY LEVEL;
      for (XLEVEL = 1; XLEVEL <= DEPTH; XLEVEL++) {
           LEVEL  = RLEVEL ? DEPTH-XLEVEL+1 : XLEVEL;
           LSTART = LVLPTR[LEVEL-1];
           LWIDTH = LVLPTR[LEVEL  ] - LSTART;
           NLEFT  = LWIDTH;
           QUEUE2 = QUEUE1 + INC*LWIDTH;
           BQ2    = QUEUE2;

///////////////////////////////////////////////////////////////////////////////////////////
//...'LOOP' CONSTRUCT BEGINS, THE INNER 'REPEAT' WILL BE DONE AS MANY TIMES AS IS NECESSARY
//    TO NUMBER ALL THE NODES AT THIS LEVEL.
rep:
           do {
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//... REPEAT ..., UNTIL QUEUE1 BECOMES EMPTY TAKE NODE FROM FRONT OF QUEUE1, FIND EACH OF ITS NEIGHBORS
//    WHICH HAVE NOT YET BEEN NUMBERED, AND ADD THE NEIGHBORS TO QUEUE1 OR QUEUE2 ACCORDING TO THEIR LEVELS.
               CNODE  = INVNUM[FQ1-1];
               FQ1   += INC;
               SQ1    = BQ1;
               SQ2    = BQ2;
               NLEFT--;
               CPTR   = RSTART[CNODE-1];
               CDGREE = DEGREE[CNODE-1];
               for (I = 1; I <= CDGREE; I++) {
                    INODE  = CONNEC[CPTR++-1];
                    ILEVEL = LVLNUM[INODE -1];
                    if (ILEVEL != 0) {
                        LVLNUM[INODE-1] = 0;
                        if (ILEVEL != LEVEL) {
                            if (std::abs(LEVEL-ILEVEL) != 1) {
                                ERGPS = 54; SPACE = -1; return;
                            }
                            INVNUM[BQ2-1] = INODE;
                            BQ2 += INC;
                        }
                        else {
                            INVNUM[BQ1-1] = INODE;
                            BQ1 += INC;
                        }
                    }
               }

////////////////////////////////////////////////////////////////////////////////////////////////
//... SORT THE NODES JUST ADDED TO QUEUE1 AND QUEUE2 SEPARATELY INTO INCREASING ORDER OF DEGREE.
               if ((NSORT = std::abs(BQ1-SQ1)) > 1) {
                   if (! FORWRD) {
                       GPSKCP (NSORT, INVNUM+BQ1, N, DEGREE, ERGPS);
                       if (ERGPS != 0) goto err;
                   }
                   else {
                       GPSKCQ (NSORT, INVNUM+SQ1-1, N, DEGREE, ERGPS);
                       if (ERGPS != 0) goto err;
                   }
               }
               if ((NSORT = std::abs(BQ2-SQ2)) > 1) {
                   if (! FORWRD) {
                       GPSKCP (NSORT, INVNUM+BQ2, N, DEGREE, ERGPS);
                       if (ERGPS != 0) goto err;
                   }
                   else {
                       GPSKCQ (NSORT, INVNUM+SQ2-1, N, DEGREE, ERGPS);
                       if (ERGPS != 0) goto err;
                   }
               }
/////////////////////////
//... END OF REPEAT LOOP;
           }
           while (FQ1 != BQ1);

/////////////////////////////////////////////////////////////////////////////////////////////////
//... QUEUE1 IS NOW EMPTY ..., IF THERE ARE ANY UNNUMBERED NODES LEFT AT THIS LEVEL, FIND THE ONE
//    OF MINIMAL DEGREE AND RETURN TO THE REPEAT LOOP ABOVE.
           if (BQ1 != QUEUE2 || NLEFT != 0) {
              if (NLEFT <= 0 || NLEFT != INC*(QUEUE2-BQ1)) {
                  ERGPS = 52; SPACE = -1; return;
              }
              for (LOWDG = BPTR = N+1, CPTR = LSTART-1, I = 1; I <= NLEFT; I++) {
                   while (LVLNUM[(LVLLSC = LVLLST[++CPTR-1])-1] != LEVEL)
                      if (LVLNUM[LVLLSC-1] != 0) {
                          ERGPS = 53; SPACE = -1; return;
                   }
                   if (DEGREE[LVLLSC-1] < LOWDG) {
                       LOWDG = DEGREE[LVLLSC-1];
                       BPTR  = CPTR;
                   }
              }

///////////////////////////////////////////
//... MINIMAL DEGREE UNNUMBERED NODE FOUND;
              if (BPTR > N) {
                  ERGPS = 55; SPACE = -1; return;
              }
              LVLLSB = LVLLST[BPTR-1];
              INVNUM[BQ1-1]    = LVLLSB;
              LVLNUM[LVLLSB-1] = 0;
              BQ1 += INC;

              goto rep;
           }

//////////////////////////////////////////////////////////////////////////////////
//... ADVANCE QUEUE POINTERS TO MAKE QUEUE2 THE NEW QUEUE1 FOR THE NEXT ITERATION.
           else {
               QUEUE1 = QUEUE2;
               FQ1    = QUEUE1;
               BQ1    = BQ2;
               if (BQ1 == FQ1 && XLEVEL < DEPTH) {
                   ERGPS = 51; SPACE = -1; return;
               }
           }
      }

/////////////////////////////////////////////////////////////
//... CHANGE SIGN OF DEGREE TO MARK THESE NODES AS 'NUMBERED';
      for (I = 1; I <= NCOMPN; I++) {
           INVNMI = INVNUM[I-1];
           DEGREE[INVNMI-1] = -DEGREE[INVNMI-1];
      }
      return;
err:
      ERGPS = 56; SPACE = -1; return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...NUMBER NODES IN A GENERALIZED LEVEL STRUCTURE ACCORDING TO A GENERALIZATION OF THE KING
//   ALGORITHM, WHICH REDUCES THE PROFILE OF THE SPARSE SYMMETRIC MATRIX.
void GPSKCK(int N,
            int * DEGREE,
            int * RSTART,
            int * CONNEC,
            int WRKLEN,
            int & NXTNUM,
            int * WORK,
            int & NCOMPN,
            int & DEPTH,
            int * LVLLST,
            int * LVLPTR,
            int * LVLNUM,
            int & ERGPS,
            int & SPACE)
{
////////////////////////////////////////////////////////////////////////
//======================================================================
// CODE USES A PRIORITY QUEUE TO CHOOSE THE NEXT NODE TO BE NUMBERED
// THE PRIORITY QUEUE IS REPRESENTED BY A SIMPLE LINEAR-LINKED LIST
// TO SAVE SPACE.  THIS WILL REQUIRE MORE SEARCHING THAN A FULLY
// LINKED REPRESENTATION, BUT THE DATA MANIPULATION IS SIMPLER.
//
//     << ESTABLISH PRIORITY QUEUE 'ACTIVE' FOR LEVEL 1 NODES >>
//     FOR I = 1 TO DEPTH DO
//         << SET QUEUE 'QUEUED' TO BE EMPTY, LIST 'NEXT' TO BE >>
//         << SET OF NODES AT NEXT LEVEL.                       >>
//         FOR J = 1 TO 'NODES AT THIS LEVEL' DO
//             << FIND FIRST NODE IN ACTIVE WITH MINIMAL CONNECTIONS >>
//             << TO 'NEXT'.  NUMBER THIS NODE AND REMOVE HIM FROM   >>
//             << 'ACTIVE'.  FOR EACH NODE IN 'NEXT' WHICH CONNECTED >>
//             << TO THIS NODE, MOVE IT TO 'QUEUED' AND REMOVE IT    >>
//             << FROM 'NEXT'.                                       >>
//         << SET NEW QUEUE 'ACTIVE' TO BE 'QUEUED' FOLLOWED BY ANY >>
//         << NODES STILL IN 'NEXT'.                                >>
//======================================================================
// DATA STRUCTURE ASSUMPTIONS:
// THE FIRST 'NXTNUM-1' ELEMENTS OF  WORK  ARE ALREADY IN USE. THE LEVEL
// STRUCTURE 'LVLLST' IS CONTIGUOUS WITH  WORK, THAT IS, IT RESIDES IN
// ELEMENTS  WRKLEN+1, ...  OF  WORK.  'LVLPTR' AND 'LVLNUM' ARE ALSO EMBEDDED
// IN WORK, BEHIND 'LVLLST'.  THE THREE VECTORS ARE PASSED SEPARATELY TO CLARIFY
// THE INDEXING, BUT THE QUEUES DEVELOPED WILL BE ALLOWED TO OVERRUN 'LVLLST' AS NEEDED.
//
//... BUILD THE FIRST 'ACTIVE' QUEUE STARTING W1 LOCATIONS FROM THE FRONT OF THE CURRENT
//    WORKING AREA  (W1 IS THE WIDTH OF THE FIRST LEVEL).  BUILD THE FIRST 'QUEUED' QUEUE
//    STARTING FROM THE BACK OF WORK SPACE.  THE LIST 'NEXT' WILL BE REALIZED
//    IMPLICITLY IN 'LVLNUM' AS:
//       LVLNUM(I) > 0   <== LEVEL NUMBER OF NODE.  'NEXT' IS SET WITH LVLNUM(I) = LEVEL+1
//       LVLNUM(I) = 0   <== I-TH NODE IS IN 'QUEUED' OR IS NOT IN THIS COMPONENT OF GRAPH, OR HAS JUST BEEN NUMBERED.
//       LVLNUM(I) < 0   <== I-TH NODE IS IN 'ACTIVE' AND IS CONNECTED TO -LVLNUM(I) NODES IN 'NEXT'.
//======================================================================
// STRUCTURE OF WORKSPACE ..
//--------------------------------------------------------------
//: NUMBERED : DONE : ACTIVE : ALEVEL : ... : QUEUED : LVLLST :
//--------------------------------------------------------------
//-------------------
// LVLPTR : LVLNUM :
//-------------------
// IN THE ABOVE, NUMBERED IS THE SET OF NODES ALREADY NUMBERED FROM PREVIOUS COMPONENTS AND EARLIER LEVELS OF
// THIS COMPONENT. DONE, ACTIVE, ALEVEL  ARE VECTORS OF LENGTH THE WIDTH OF THE CURRENT LEVEL.
// ACTIVE IS A SET OF INDICES INTO ALEVEL.  AS THE NODES IN ALEVEL ARE NUMBERED, THEY  ARE PLACED INTO 'DONE'.
// QUEUED IS A QUEUE OF NODES IN THE 'NEXT' LEVEL, WHICH GROWS FROM THE START OF THE 'NEXT' LEVEL IN LVLLST
// FORWARDS TOWARD 'ALEVEL'.  QUEUED IS OF LENGTH NO MORE THAN THE WIDTH OF THE NEXT LEVEL.
// LVLLST IS THE LIST OF UNNUMBERED NODES IN THE TREE, ARRANGED BY LEVEL.
//======================================================================
//    INTEGER   LVLLST(N), LVLPTR(DEPTH), LVLNUM(N)
      int I, J, K, PTR, JPTR, KPTR, lptr, MPTR, PPTR, RPTR, MPPTR, JNODE, KNODE, CNODE, LEVEL,
          LOWDG, UNUSED, MXQUE, NNEXT, ASTART, MINDG, LSTART, LWIDTH, ACTIVE, QUEUEB, QUEUED,
          QCOUNT, NCONNC, NACTIV, CDGREE, LDGREE, NFINAL, JDGREE, STRTIC, ADDED, TWRKLN, LVLLSL,
          CONNER, ACTIVI, BREAK;

      TWRKLN = WRKLEN + NCOMPN + N + DEPTH + 1;
      UNUSED = TWRKLN;

      ASTART = LVLPTR[0];
      LWIDTH = LVLPTR[1] - ASTART;
      ASTART = WRKLEN  + 1;
      ACTIVE = NXTNUM + LWIDTH + 1;
      NACTIV = LWIDTH;
      NFINAL = NXTNUM + NCOMPN;

      NNEXT  = LVLPTR[2] - LVLPTR[1];
      QUEUED = WRKLEN;
      QUEUEB = QUEUED;
      MXQUE  = ACTIVE + LWIDTH;

//////////////////////////////////////////
//... BUILD FIRST PRIORITY QUEUE 'ACTIVE';
      for (LOWDG  = -(N+1), lptr = LVLPTR[0], I = 1; I <= LWIDTH; I++) {
           LVLLSL =  LVLLST[lptr  -1];
           JPTR   =  RSTART[LVLLSL-1];
           LDGREE =  DEGREE[LVLLSL-1];

           for (NCONNC = 0, J = 1; J <= LDGREE; J++)
            if (LVLNUM[CONNEC[JPTR++ -1]-1] == 2) NCONNC--;

           ACTIVI = ACTIVE+I-1;
           WORK  [ACTIVI-1] = I;
           LVLNUM[LVLLSL-1] = NCONNC;
           LOWDG = max(LOWDG, NCONNC);
           lptr++;
      }
      WORK[ACTIVE-2] = 0;

/////////////////////////////////////////
//... NOW NUMBER NODES LEVEL BY LEVEL...;
      for (LEVEL = 1; LEVEL <= DEPTH; LEVEL++) {

/////////////////////////////////////
//... NUMBER ALL NODES IN THIS LEVEL;
          for (I = 1; I <= LWIDTH; I++) {
              PPTR = -1;
              PTR  = WORK[ACTIVE-2];

//////////////////////////////////////////////////////////////////////////////////
//... IF NODES REMAIN IN NEXT, FIND THE EARLIEST NODE IN ACTIVE OF MINIMAL DEGREE.
              if (NNEXT != 0) {
                  for (MINDG = -(N+1), BREAK = 0, J = 1; J <= NACTIV; J++) {
                       if (LVLNUM[(CNODE = WORK[(ASTART+PTR)-1])-1] == LOWDG) {
                           BREAK = 1; break;
                       }
                       if (LVLNUM[CNODE-1] > MINDG) {
                           MPPTR = PPTR;
                           MPTR  = PTR;
                           MINDG = LVLNUM[CNODE-1];
                       }
                       PPTR = PTR;
                       PTR  = WORK[(ACTIVE+PTR)-1];
                  }
/////////////////////////////////////////////////////////////////////////
//... ESTABLISH PTR AS FIRST MIN DEGREE NODE PPTR AS PREDECESSOR IN LIST.
                  if (! BREAK) {
                      PTR  = MPTR;
                      PPTR = MPPTR;
                  }
                  CNODE = WORK  [(ASTART+PTR)-1];
                  LOWDG = LVLNUM[CNODE-1];
                  LVLNUM[CNODE-1] = 0;
                  JPTR  = RSTART[CNODE-1];

///////////////////////////////////////////////////////////////////////////////////////////
//... UPDATE CONNECTION COUNTS FOR ALL NODES WHICH CONNECT TO  CNODE'S  NEIGHBORS IN  NEXT.
                  for (CDGREE = DEGREE[CNODE-1], STRTIC = QUEUEB, J = 1; J <= CDGREE; J++)
                   if (LVLNUM[(JNODE = CONNEC[JPTR++-1])-1] == LEVEL+1 ) {
                       if (QUEUEB < MXQUE) goto err;

                       WORK[QUEUEB-- -1] = JNODE;
                       NNEXT--;
                       LVLNUM[JNODE-1]   = 0;

                       if (NACTIV != 1)
                       for (KPTR  = RSTART[JNODE-1], JDGREE = DEGREE[JNODE-1], K = 1; K <= JDGREE; K++) {
                            KNODE = CONNEC[KPTR++-1];
                            if (LVLNUM[KNODE-1] < 0)
                            if (LOWDG < ++LVLNUM[KNODE-1]) LOWDG = LVLNUM[KNODE-1];
                       }
                  }

////////////////////////////////////////////////////////////////////////////////////////////////////
//... TO MIMIC THE ALGORITHM AS IMPLEMENTED BY GIBBS, SORT THE NODES JUST ADDED TO THE QUEUE INTO
//    INCREASING ORDER OF ORIGINAL INDEX. (BUT, BECAUSE THE QUEUE IS STORED BACKWARDS IN MEMORY, THE
//    SORT ROUTINE IS CALLED FOR DECREASING INDEX.) TREAT  0, 1 OR 2  NODES ADDED AS SPECIAL CASES
                  if ((ADDED = STRTIC - QUEUEB) > 2) {
                     GPSKCO(ADDED, WORK+QUEUEB, ERGPS);
                     if (ERGPS != 0) {
                         ERGPS = 64; return;
                     }
                  }
                  else if (ADDED == 2   && WORK[STRTIC-2] <= WORK[STRTIC-1])
                     iSWAP(WORK[STRTIC-1], WORK[STRTIC-2]);

              }

///////////////////////////////////////////////////////////////////////////////////////////////////
//... NUMBER THIS NODE AND DELETE IT FROM 'ACTIVE'. MARK IT UNAVAILABLE BY CHANGING SIGN OF DEGREE;
              NACTIV--;
              CNODE = WORK[(ASTART+PTR)-1];
              WORK[NXTNUM++-1] =  CNODE;
              DEGREE[CNODE -1] = -DEGREE[CNODE-1];

/////////////////////////////////////////
//... DELETE LINK TO THIS NODE FROM LIST;
              WORK[(ACTIVE+PPTR)-1] = WORK[(ACTIVE+PTR)-1];
          }

///////////////////////////////////////////////////////////////////////////////////////////////
//... NOW MOVE THE QUEUE 'QUEUED' FORWARD, AT THE SAME TIME COMPUTING CONNECTION COUNTS FOR ITS
//    ELEMENTS. THEN DO THE SAME FOR THE REMAINING NODES IN 'NEXT'.
          UNUSED = min(UNUSED, QUEUEB-MXQUE);
          if (NXTNUM != ACTIVE-1) {
              ERGPS = 61; return;
          }
          if (LEVEL != DEPTH) {
              LSTART = LVLPTR[LEVEL];
              LWIDTH = LVLPTR[LEVEL+1] - LSTART;
              ACTIVE = NXTNUM + LWIDTH + 1;
              ASTART = ACTIVE + LWIDTH;
              NACTIV = LWIDTH;
              if ((MXQUE = ASTART+LWIDTH) > QUEUEB+1) goto err;
              UNUSED = min(UNUSED, QUEUEB-MXQUE+1);

              QCOUNT = QUEUED - QUEUEB;
              LOWDG  = -N-1;
              WORK[ACTIVE-2] = 0;

//////////////////////////////////////////////////////
//... CHOOSE NEXT NODE FROM EITHER 'QUEUED' OR 'NEXT';
              for (PTR = LSTART, I = 1; I <= LWIDTH; I++) {
                  if (I <= QCOUNT) CNODE = WORK[(QUEUED+1-I)-1];
                  else do {
                      CNODE = LVLLST[PTR++-1];
                      if (PTR > LVLPTR[LEVEL+1]) {
                          ERGPS = 62; return;
                      }
                  } while (LVLNUM[CNODE-1] <= 0);

                  if (LEVEL+1 != DEPTH) {
                      for (RPTR  = RSTART[CNODE-1], NCONNC = 0, JDGREE = DEGREE[CNODE-1], J = 1; J <= JDGREE; J++) {
                          CONNER = CONNEC[RPTR++-1];
                          if (LVLNUM[CONNER-1] == LEVEL+2) NCONNC--;
                      }
                      LVLNUM[CNODE-1]  = NCONNC;
                      LOWDG = max(LOWDG, NCONNC);
                  }
//////////////////////////////////////
//... ADD CNODE TO NEW 'ACTIVE' QUEUE;
                  WORK[(ACTIVI = ACTIVE+(I-1))-1] = I;
                  WORK[(         ASTART+(I-1))-1] = CNODE;
              }
              if (DEPTH == LEVEL+1) NNEXT = 0; else {
                  NNEXT  = LVLPTR[LEVEL+2] - LVLPTR[LEVEL+1];
                  QUEUED = LSTART - 1 + LWIDTH + WRKLEN;
                  QUEUEB = QUEUED;
              }
          }
      }
      if (NXTNUM != NFINAL) {
          ERGPS = 63; return;
      }
      SPACE = max(SPACE, TWRKLN-UNUSED); return;
err:
      ERGPS = 160; SPACE = NACTIV+NNEXT; return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...COMPUTE THE BANDWIDTH AND PROFILE FOR THE RENUMBERING GIVEN BY 'INVNUM' AND ALSO FOR THE RENUMBERING
//   GIVEN BY 'OLDNUM'. 'NEWNUM' WILL BE A PERMUTATION VECTOR COPY OF THE NODE LIST 'INVNUM'.
void GPSKCL(int N,
            int * DEGREE,
            int * RSTART,
            int * CONNEC,
            int * INVNUM,
            int * NEWNUM,
            int * OLDNUM,
            int & BANDWD,
            int & PROFIL,
            int & ERGPS,
            int & SPACE)
{
//    INTEGER   INVNUM(N), NEWNUM(N), OLDNUM(N)
      int I, J, JPTR,  IDGREE, OLDBND, OLDPRO, NEWBND, NEWPRO, OLDRWD, NEWRWD, OLDORG, NEWORG,
                JNODE;

////////////////////////////////////////////
//... CREATE NEWNUM AS A PERMUTATION VECTOR;
      for (I = 1; I <= N; I++) NEWNUM[(INVNUM[I-1])-1] = I;

///////////////////////////////////////////////////////////////////////////
//... COMPUTE PROFILE AND BANDWIDTH FOR BOTH THE OLD AND THE NEW ORDERINGS.
      for (OLDBND = OLDPRO = NEWBND = NEWPRO = 0, I = 1; I <= N; I++)
      if  (DEGREE[I-1] != 0) {
           if (DEGREE[I-1] > 0) goto err;

           IDGREE = -DEGREE[I-1];
           DEGREE[I-1] = IDGREE;
           NEWORG =  NEWNUM[I-1];
           OLDORG =  OLDNUM[I-1];
           NEWRWD =  0;
           OLDRWD =  0;
           JPTR   =  RSTART[I-1];
/////////////////////////////////////////////////////////////////////////
//... FIND NEIGHBOR WHICH IS NUMBERED FARTHEST AHEAD OF THE CURRENT NODE.
              for (J = 1; J <= IDGREE; J++) {
                   JNODE  = CONNEC[JPTR++-1];
                   NEWRWD = max(NEWRWD, NEWORG - NEWNUM[JNODE-1]);
                   OLDRWD = max(OLDRWD, OLDORG - OLDNUM[JNODE-1]);
              }
              NEWPRO += NEWRWD;
              NEWBND  = max(NEWBND, NEWRWD);
              OLDPRO += OLDRWD;
              OLDBND  = max(OLDBND, OLDRWD);
      }

///////////////////////////////////////////////////////////////////////////////////////////////////
//... IF NEW ORDERING HAS BETTER BANDWIDTH THAN OLD ORDERING, REPLACE OLD ORDERING BY NEW ORDERING;
      if (NEWBND <= OLDBND) {
          BANDWD = NEWBND;
          PROFIL = NEWPRO;
          for (I = 1; I <= N; I++) OLDNUM[I-1] = NEWNUM[I-1];
      }

//////////////////////////
//... RETAIN OLD ORDERING;
      else {
          BANDWD = OLDBND;
          PROFIL = OLDPRO;
      }
      return;
err:
      ERGPS = 70; SPACE = -1; return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...COMPUTE THE BANDWIDTH AND PROFILE FOR THE RENUMBERING GIVEN BY 'INVNUM', BY THE REVERSE OF NUMBERING
//  'INVNUM', AND ALSO BY THE RENUMBERING GIVEN IN 'OLDNUM'. 'NEWNUM' WILL BE A PERMUTATION VECTOR COPY OF
//   THE NODE LIST 'INVNUM'.
void GPSKCM(int N,
            int * DEGREE,
            int * RSTART,
            int * CONNEC,
            int * INVNUM,
            int * NEWNUM,
            int * OLDNUM,
            int & BANDWD,
            int & PROFIL,
            int & ERGPS,
            int & SPACE)
{
//    INTEGER     DEGREE(N), CONNEC(N), INVNUM(N), NEWNUM(N), OLDNUM(N)
      int   I, J, JPTR, IDGREE, OLDBND, OLDPRO, NEWBND, NEWPRO, OLDRWD, NEWRWD, OLDORG, NEWORG,
               JNODE, NRVBND, NRVPRO, NRVORG, NRVRWD, NMIP1;

////////////////////////////////////////////
//... CREATE NEWNUM AS A PERMUTATION VECTOR;
      for (I = 1; I <= N; I++) NEWNUM[(INVNUM[I-1])-1] = I;

///////////////////////////////////////////////////////////////////////////
//... COMPUTE PROFILE AND BANDWIDTH FOR BOTH THE OLD AND THE NEW ORDERINGS.
      for (OLDBND = OLDPRO = NEWBND = NEWPRO = NRVBND = NRVPRO = 0, I = 1; I <= N; I++)
      if  (DEGREE[I-1]) {
           if (DEGREE[I-1] > 0)  goto err;

           IDGREE = -DEGREE[I-1];
           DEGREE[I-1] = IDGREE;
           NEWRWD =
           OLDRWD =
           NRVRWD =  0;
           NEWORG =  NEWNUM[I-1];
           OLDORG =  OLDNUM[I-1];
           NRVORG =  N - NEWNUM[I-1] + 1;
           JPTR   =  RSTART[I-1];

/////////////////////////////////////////////////////////////////////////
//... FIND NEIGHBOR WHICH IS NUMBERED FARTHEST AHEAD OF THE CURRENT NODE.
              for (J = 1; J <= IDGREE; J++) {
                   JNODE  = CONNEC[JPTR++-1];
                   NEWRWD = max(NEWRWD, NEWORG - NEWNUM[JNODE-1]);
                   OLDRWD = max(OLDRWD, OLDORG - OLDNUM[JNODE-1]);
                   NRVRWD = max(NRVRWD, NRVORG - N + NEWNUM[JNODE-1] - 1);
              }
              NEWPRO += NEWRWD;
              NEWBND  = max(NEWBND, NEWRWD);
              NRVPRO += NRVRWD;
              NRVBND  = max(NRVBND, NRVRWD);
              OLDPRO += OLDRWD;
              OLDBND  = max(OLDBND, OLDRWD);
      }

///////////////////////////////////////////////////////////////////////////////////////////////////
//... IF NEW ORDERING HAS BETTER BANDWIDTH THAN OLD ORDERING, REPLACE OLD ORDERING BY NEW ORDERING;
      if (NEWPRO <= OLDPRO && NEWPRO <= NRVPRO) {
          BANDWD = NEWBND;
          PROFIL = NEWPRO;
          for (I = 1; I <= N; I++) OLDNUM[I-1] = NEWNUM[I-1];
      }
      else
///////////////////////////////////////////////////
//... CHECK NEW REVERSED ORDERING FOR BEST PROFILE;
      if (NRVPRO <= OLDPRO) {
          BANDWD = NRVBND;
          PROFIL = NRVPRO;
          for (I = 1;  I <= N; I++) {
              OLDNUM[I-1] = N - NEWNUM[I-1] + 1;
              if (I <= N/2) {
                  J     = INVNUM[I-1];
                  NMIP1 = (N + 1) - I;
                  INVNUM[I-1]     = INVNUM[NMIP1-1];
                  INVNUM[NMIP1-1] = J;
              }
          }
      }

//////////////////////////
//... RETAIN OLD ORDERING;
      else {
          BANDWD = OLDBND;
          PROFIL = OLDPRO;
      }
      return;
err:
      ERGPS = 71; SPACE = -1; return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//... I N S E R T I O N    S O R T
void GPSKCN(int N,       // -- NUMBER OF ELEMENTS TO BE SORTED.
            int * KEY,   // -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES WHICH ARE TO BE SORTED.
            int * DATA,  // -- A SECOND ARRAY OF LENGTH  N  CONTAINING DATA ASSOCIATED WITH THE INDIVIDUAL KEYS.
            int & ERGPS) // -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING, IN WHICH CASE IT WILL BE EQUAL TO 1.
{
////////////////////////////////////////////////////////////////////
//==================================================================
// KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING ORDER
// DATA -- REARRANGED TO CORRESPOND TO REARRANGED KEYS
//==================================================================
//    INTEGER  KEY(N), DATA(N)
      int   I, J, D, K, IP1, JM1;
      if (N == 1) return;
      if (N <= 0) goto err;

//////////////////////////////////////////////////////////
//... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...;
      ERGPS = 0;
      I     = N - 1;
      IP1   = N;
      do {
////////////////////////////////////////////////
//... OUT OF ORDER ... MOVE UP TO CORRECT PLACE;
          if (KEY[I-1] < KEY[IP1-1]) {
              K   = KEY [I-1];
              D   = DATA[I-1];
              J   = IP1;
              JM1 = I;

///////////////////////////////////////////////////
//... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND';
              do {
                  KEY [JM1-1] = KEY [J-1];
                  DATA[JM1-1] = DATA[J-1];
                  JM1         = J++;
              }
              while (J <= N && KEY[J-1] >= K);

              KEY [JM1-1] = K;
              DATA[JM1-1] = D;
          }
          IP1 = I--;
      }
      while (I > 0);

      return;
err:
      ERGPS = 1; return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//... I N S E R T I O N    S O R T
void GPSKCO(int N,       // -- NUMBER OF ELEMENTS TO BE SORTED.
            int * KEY,   // -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES WHICH ARE TO BE SORTED.
            int & ERGPS) //
{
////////////////////////////////////////////////////////////////////
//==================================================================
//  KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING ORDER
//==================================================================
//    INTEGER    KEY(N)
      int     I, J, K, IP1, JM1;
      if (N == 1) return;
      if (N <= 0) goto err;

//////////////////////////////////////////////////////////
//... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...;
      ERGPS = 0;
      I     = N - 1;
      IP1   = N;
      do {
////////////////////////////////////////////////
//... OUT OF ORDER ... MOVE UP TO CORRECT PLACE;
          if (KEY[I-1] < KEY[IP1-1]) {
              K   = KEY[I-1];
              J   = IP1;
              JM1 = I;

///////////////////////////////////////////////////
//... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND';
              do {
                  KEY[JM1-1] = KEY[J-1];
                  JM1        = J++;
              }
              while (J <= N && KEY[J-1] > K);

              KEY[JM1-1] = K;
          }
          IP1 = I--;
      }
      while (I > 0);

      return;
err:
      ERGPS = 1; return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//... I N S E R T I O N    S O R T
void GPSKCP(int N,        // -- NUMBER OF ELEMENTS TO BE SORTED.
            int * INDEX,  // -- AN ARRAY OF LENGTH  N  CONTAINING THE INDICES WHOSE DEGREES ARE TO BE SORTED.
            int & NVEC,   //
            int * DEGREE, // -- AN  NVEC  VECTOR, GIVING THE DEGREES OF NODES WHICH ARE TO BE SORTED.
            int & ERGPS)  // -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING, IN WHICH CASE IT WILL BE EQUAL TO 1.
{
////////////////////////////////////////////////////////////////////
//==================================================================
// INDEX  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING ORDER
//==================================================================
//    INTEGER     INDEX(N), DEGREE(NVEC)
      int   I, J, V, IP1, JM1, INDEXI, INDXI1, INDEXJ; NVEC = NVEC;
      if (N == 1) return;
      if (N <= 0) goto err;

//////////////////////////////////////////////////////////
//... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...;
      ERGPS = 0;
      I     = N - 1;
      IP1   = N;
      do {
          INDEXI = INDEX[I-1];
          INDXI1 = INDEX[IP1-1];

////////////////////////////////////////////////
//... OUT OF ORDER ... MOVE UP TO CORRECT PLACE;
          if (DEGREE[INDEXI-1] < DEGREE[INDXI1-1]) {
              V      = DEGREE[INDEXI-1];
              J      = IP1;
              JM1    = I;
              INDEXJ = INDEX[J-1];

///////////////////////////////////////////////////
//... REPEAT ... UNTIL 'CORRECT PLACE FOR V FOUND';
              do {
                  INDEX[JM1-1] = INDEXJ;
                  JM1 = J++;
              }
              while (J <= N && DEGREE[(INDEXJ = INDEX[J-1])-1] > V);

              INDEX [JM1-1] = INDEXI;
          }
          IP1 = I--;
      }
      while (I > 0);

      return;
err:
      ERGPS = 1; return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//... I N S E R T I O N    S O R T
void GPSKCQ(int N,        // -- NUMBER OF ELEMENTS TO BE SORTED.
            int * INDEX,  // -- AN ARRAY OF LENGTH  N  CONTAINING THE INDICES WHOSE DEGREES ARE TO BE SORTED.
            int & NVEC,   //
            int * DEGREE, // -- AN  NVEC  VECTOR, GIVING THE DEGREES OF NODES WHICH ARE TO BE SORTED.
            int & ERGPS)  // -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING, IN WHICH CASE IT WILL BE EQUAL TO 1.
{
////////////////////////////////////////////////////////////////////
//==================================================================
// INDEX  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING ORDER
//==================================================================
//    INTEGER     INDEX(N), DEGREE(NVEC)
      int   I, J, V, INDEXI, INDXI1, INDEXJ, IP1, JM1; NVEC = NVEC;
      if (N == 1) return;
      if (N <= 0) goto err;

//////////////////////////////////////////////////////////
//... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...;
      ERGPS = 0;
      I     = N - 1;
      IP1   = N;
      do {
          INDEXI = INDEX[I-1];
          INDXI1 = INDEX[IP1-1];
////////////////////////////////////////////////
//... OUT OF ORDER ... MOVE UP TO CORRECT PLACE;
          if (DEGREE[INDEXI-1] > DEGREE[INDXI1-1]) {
              V      = DEGREE[INDEXI-1];
              J      = IP1;
              JM1    = I;
              INDEXJ = INDEX[J-1];

///////////////////////////////////////////////////
//... REPEAT ... UNTIL 'CORRECT PLACE FOR V FOUND';
              do {
                  INDEX[JM1-1] = INDEXJ;
                  JM1 = J++;
              }
              while (J <= N && DEGREE[(INDEXJ = INDEX[J-1])-1] < V);

              INDEX [JM1-1] = INDEXI;
          }
          IP1 = I--;
      }
      while (I > 0);

      return;
err:
      ERGPS = 1; return;
}
