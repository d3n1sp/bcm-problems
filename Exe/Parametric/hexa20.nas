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
$   Date       : Sun Feb 25 09:47:50 2007
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
GRID           2       0      1.      0.      1.       0        
GRID           3       0      1.      0.      0.       0        
GRID           4       0      0.      1.      0.       0        
GRID           5       0      0.      1.      1.       0        
GRID           6       0      1.      1.      1.       0        
GRID           7       0      1.      1.      0.       0        
GRID           8       0      0.      0.      0.       0        
GRID           9       0     0.5      0.      1.       0        
GRID          10       0      0.     0.5      1.       0        
GRID          11       0      1.     0.5      1.       0        
GRID          12       0     0.5      1.      1.       0        
GRID          13       0     0.5      0.      0.       0        
GRID          14       0      0.     0.5      0.       0        
GRID          15       0      1.     0.5      0.       0        
GRID          16       0     0.5      1.      0.       0        
GRID          17       0      0.      0.     0.5       0        
GRID          18       0      1.      0.     0.5       0        
GRID          19       0      0.      1.     0.5       0        
GRID          20       0      1.      1.     0.5       0        
CHEXA          7       1       1       5       6       2       8       4+EL    7
+EL    7       7       3      10      12      11       9      17      19+EA    7
+EA    7      20      18      14      16      15      13
ENDDATA
