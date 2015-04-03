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
$   Date       : Fri Apr 09 10:24:13 2010
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
GRID           2       0     0.5      1.      1.       0        
GRID           3       0      0.      1.      1.       0        
GRID           4       0      0.      1.      0.       0        
GRID           5       0     0.5      1.     0.5       0        
GRID           6       0      0.      1.     0.5       0        
GRID           7       0      0.      0.     0.5       0        
GRID           8       0     0.5      1.      0.       0        
GRID           9       0      0.      0.      0.       0        
GRID          10       0     0.5      0.     0.5       0        
GRID          11       0      1.      1.     0.5       0        
GRID          12       0      1.      0.      1.       0        
GRID          13       0     0.5      0.      1.       0        
GRID          14       0      1.      1.      1.       0        
GRID          15       0      1.      0.     0.5       0        
GRID          16       0      1.      1.      0.       0        
GRID          17       0     0.5      0.      0.       0        
GRID          18       0      1.      0.      0.       0        
CHEXA          1       1       6       3       1       7       5       2+EL    1
+EL    1      13      10                                                        
CHEXA          2       1      10       7       6       5      17       9+EL    2
+EL    2       4       8                                                        
CHEXA          3       1      14      12      13       2      11      15+EL    3
+EL    3      10       5                                                        
CHEXA          4       1       5      10      17       8      11      15+EL    4
+EL    4      18      16                                                        
ENDDATA
