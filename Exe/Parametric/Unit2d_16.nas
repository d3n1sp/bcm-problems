ID C:\Users\Dima\Lame3d2\Mpl,FEMAP
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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Parametric\Solid\UnitBox.MOD
$   Date       : Mon Oct 10 15:28:47 2005
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
$ FEMAP Property 1 : UnitPlane
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : MaterialNull
MAT1           1                              0.      0.      0.        
GRID           5       0      0.      1.      0.       0        
GRID           6       0      0.    0.75      0.       0        
GRID           7       0      0.     0.5      0.       0        
GRID           8       0      0.    0.25      0.       0        
GRID           9       0      0.      0.      0.       0        
GRID          10       0    0.25      0.      0.       0        
GRID          11       0     0.5      0.      0.       0        
GRID          12       0    0.75      0.      0.       0        
GRID          13       0      1.      0.      0.       0        
GRID          14       0      1.    0.25      0.       0        
GRID          15       0      1.     0.5      0.       0        
GRID          16       0      1.    0.75      0.       0        
GRID          17       0      1.      1.      0.       0        
GRID          18       0    0.75      1.      0.       0        
GRID          19       0     0.5      1.      0.       0        
GRID          20       0    0.25      1.      0.       0        
GRID          21       0    0.25    0.75      0.       0        
GRID          22       0    0.25     0.5      0.       0        
GRID          23       0    0.25    0.25      0.       0        
GRID          24       0     0.5    0.75      0.       0        
GRID          25       0     0.5     0.5      0.       0        
GRID          26       0     0.5    0.25      0.       0        
GRID          27       0    0.75    0.75      0.       0        
GRID          28       0    0.75     0.5      0.       0        
GRID          29       0    0.75    0.25      0.       0        
CQUAD4         2       1       5       6      21      20                
CQUAD4         3       1       6       7      22      21                
CQUAD4         4       1       7       8      23      22                
CQUAD4         5       1       8       9      10      23                
CQUAD4         6       1      20      21      24      19                
CQUAD4         7       1      21      22      25      24                
CQUAD4         8       1      22      23      26      25                
CQUAD4         9       1      23      10      11      26                
CQUAD4        10       1      19      24      27      18                
CQUAD4        11       1      24      25      28      27                
CQUAD4        12       1      25      26      29      28                
CQUAD4        13       1      26      11      12      29                
CQUAD4        14       1      18      27      16      17                
CQUAD4        15       1      27      28      15      16                
CQUAD4        16       1      28      29      14      15                
CQUAD4        17       1      29      12      13      14                
ENDDATA
