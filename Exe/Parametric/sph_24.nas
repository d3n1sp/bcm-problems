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
$   Date       : Sun Jun 11 08:12:45 2006
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
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -1.      0.      0.       0        
GRID           2       0-0.70711      0. 0.70711       0        
GRID           3       0      0. 0.70711 0.70711       0        
GRID           4       0-0.70711 0.70711      0.       0        
GRID           5       0    -0.5 0.61237 0.61237       0        
GRID           6       0      0.      0.     -1.       0        
GRID           7       0      0.-0.70711-0.70711       0        
GRID           8       0    -0.5-0.61237-0.61237       0        
GRID           9       0-0.70711-0.70711      0.       0        
GRID          10       0      0.     -1.      0.       0        
GRID          11       0      0.-0.70711 0.70711       0        
GRID          12       0      0.      0.      1.       0        
GRID          13       0    -0.5-0.61237 0.61237       0        
GRID          14       0      0.      1.      0.       0        
GRID          15       0      0. 0.70711-0.70711       0        
GRID          16       0-0.70711      0.-0.70711       0        
GRID          17       0    -0.5 0.61237-0.61237       0        
GRID          18       0 0.70711 0.70711      0.       0        
GRID          19       0     0.5 0.61237-0.61237       0        
GRID          20       0 0.70711-0.70711      0.       0        
GRID          21       0      1.      0.      0.       0        
GRID          22       0 0.70711      0. 0.70711       0        
GRID          23       0     0.5-0.61237 0.61237       0        
GRID          24       0 0.70711      0.-0.70711       0        
GRID          25       0     0.5-0.61237-0.61237       0        
GRID          26       0     0.5 0.61237 0.61237       0        
CQUAD4         1       1       3      14       4       5                
CQUAD4         2       1       2      12       3       5                
CQUAD4         3       1       1       2       5       4                
CQUAD4         4       1       7      10       9       8                
CQUAD4         5       1      16       6       7       8                
CQUAD4         6       1       1      16       8       9                
CQUAD4         7       1      11      12       2      13                
CQUAD4         8       1       9      10      11      13                
CQUAD4         9       1       1       9      13       2                
CQUAD4        10       1      15       6      16      17                
CQUAD4        11       1       4      14      15      17                
CQUAD4        12       1       1       4      17      16                
CQUAD4        13       1      18      21      24      19                
CQUAD4        14       1      15      14      18      19                
CQUAD4        15       1       6      15      19      24                
CQUAD4        16       1      20      21      22      23                
CQUAD4        17       1      11      10      20      23                
CQUAD4        18       1      12      11      23      22                
CQUAD4        19       1      24      21      20      25                
CQUAD4        20       1       7       6      24      25                
CQUAD4        21       1      10       7      25      20                
CQUAD4        22       1      22      21      18      26                
CQUAD4        23       1       3      12      22      26                
CQUAD4        24       1      14       3      26      18                
ENDDATA
