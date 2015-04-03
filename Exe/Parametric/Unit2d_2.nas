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
$   Date       : Mon Oct 10 15:14:44 2005
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
GRID           6       0      0.     0.5      0.       0        
GRID           7       0      0.      0.      0.       0        
GRID           9       0      1.      0.      0.       0        
GRID          10       0      1.     0.5      0.       0        
GRID          11       0      1.      1.      0.       0        
CQUAD4         2       1       5       6      10      11                
CQUAD4         3       1       6       7       9      10                
ENDDATA
