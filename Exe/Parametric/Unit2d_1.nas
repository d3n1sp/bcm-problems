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
$   Date       : Mon Oct 10 14:58:03 2005
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
GRID           1       0      0.      1.      0.       0        
GRID           2       0      0.      0.      0.       0        
GRID           3       0      1.      0.      0.       0        
GRID           4       0      1.      1.      0.       0        
CQUAD4         1       1       1       2       3       4                
ENDDATA
