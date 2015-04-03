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
$   Date       : Mon Feb 26 18:50:21 2007
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
PSHELL         1       1      0.       1      1.       1 0.83333      0.
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
CQUAD4         1       1       1       2       6       5                
CQUAD4         2       1       3       2       1       8                
CQUAD4         3       1       8       1       5       4                
CQUAD4         4       1       4       5       6       7                
CQUAD4         5       1       7       3       8       4                
CQUAD4         6       1       3       7       6       2                
ENDDATA
