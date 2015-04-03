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
$   Date       : Fri Apr 09 10:21:10 2010
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
GRID           1       0      1.      0.      1.       0        
GRID           2       0      1.      0.      0.       0        
GRID           3       0      1.      1.      1.       0        
GRID           4       0      1.      1.      0.       0        
GRID           5       0     0.5      1.      0.       0        
GRID           6       0      0.      1.      0.       0        
GRID           7       0     0.5      0.      0.       0        
GRID           8       0     0.5      1.      1.       0        
GRID           9       0      0.      0.      1.       0        
GRID          10       0      0.      1.      1.       0        
GRID          11       0      0.      0.      0.       0        
GRID          12       0     0.5      0.      1.       0 
CHEXA          1       1       1       3       4       2      12       8+EL    1
+EL    1       5       7                                                        
CHEXA          2       1      12       8       5       7       9      10+EL    2
+EL    2       6      11                                                        
ENDDATA
