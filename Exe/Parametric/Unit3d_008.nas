ID FEMAP,FEMAP
SOL SESTATICS
TIME 10000
CEND
  ECHO = NONE
  DISPLACEMENT(PLOT) = ALL
  SPCFORCE(PLOT) = ALL
  OLOAD(PLOT) = ALL
  FORCE(PLOT,CORNER) = ALL
  STRESS(PLOT,CORNER) = ALL
BEGIN BULK
$ ***************************************************************************
$   Written by : FEMAP
$   Version    : 7.00
$   Translator : MSC/NASTRAN
$   From Model : 
$   Date       : Fri Apr 09 10:26:14 2010
$ ***************************************************************************
$
PARAM,POST,-1
PARAM,OGEOM,NO
PARAM,AUTOSPC,YES
PARAM,GRDPNT,0
CORD2C         1       0      0.      0.      0.      0.      0.      1.+FEMAPC1
+FEMAPC1      1.      0.      1.
CORD2S         2       0      0.      0.      0.      0.      0.      1.+FEMAPC2
+FEMAPC2      1.      0.      1.
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.      1.       0        
GRID           2       0      0.     0.5     0.5       0        
GRID           3       0      0.     0.5      1.       0        
GRID           4       0      0.      0.     0.5       0        
GRID           5       0      0.      0.      0.       0        
GRID           6       0     0.5      0.      1.       0        
GRID           7       0      1.      0.      1.       0        
GRID           8       0      1.      0.     0.5       0        
GRID           9       0     0.5     0.5     0.5       0        
GRID          10       0     0.5      0.     0.5       0        
GRID          11       0      1.      0.      0.       0        
GRID          12       0     0.5      0.      0.       0        
GRID          13       0      0.      1.      1.       0        
GRID          14       0     0.5      1.      1.       0        
GRID          15       0     0.5     0.5      1.       0        
GRID          16       0     0.5      1.     0.5       0        
GRID          17       0      0.      1.     0.5       0        
GRID          18       0      0.      1.      0.       0        
GRID          19       0      0.     0.5      0.       0        
GRID          20       0     0.5     0.5      0.       0        
GRID          21       0      1.     0.5     0.5       0        
GRID          22       0      1.      1.      1.       0        
GRID          23       0      1.     0.5      1.       0        
GRID          24       0      1.      1.      0.       0        
GRID          25       0      1.      1.     0.5       0        
GRID          26       0      1.     0.5      0.       0        
GRID          27       0     0.5      1.      0.       0        
CHEXA          1       1      10       4       1       6       9       2+EL    1
+EL    1       3      15                                                        
CHEXA          2       1      12       5       4      10      20      19+EL    2
+EL    2       2       9                                                        
CHEXA          3       1      23       7       6      15      21       8+EL    3
+EL    3      10       9                                                        
CHEXA          4       1       9      10      12      20      21       8+EL    4
+EL    4      11      26                                                        
CHEXA          5       1      17      13       3       2      16      14+EL    5
+EL    5      15       9                                                        
CHEXA          6       1      17      16       9       2      18      27+EL    6
+EL    6      20      19                                                        
CHEXA          7       1      21       9      15      23      25      16+EL    7
+EL    7      14      22                                                        
CHEXA          8       1      20       9      21      26      27      16+EL    8
+EL    8      25      24                                                        
ENDDATA
