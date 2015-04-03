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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Parametric\Solid\sphere1.MOD
$   Date       : Sat Nov 19 16:27:51 2005
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
$ FEMAP Property 1 : Sphere surface
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : Mat1
MAT1           1                              0.      0.      0.        
GRID          57       0     -1.      0.      0.       0        
GRID          58       0      0.      0.      1.       0        
GRID          59       0      0.      1.      0.       0        
GRID          61       0      0.      0.     -1.       0        
GRID          64       0      0.     -1.      0.       0        
GRID          74       0      1.      0.      0.       0        
CTRIA3        25       1      57      58      59                        
CTRIA3        26       1      57      61      64                        
CTRIA3        27       1      57      64      58                        
CTRIA3        28       1      57      59      61                        
CTRIA3        29       1      61      59      74                        
CTRIA3        30       1      58      64      74                        
CTRIA3        31       1      64      61      74                        
CTRIA3        32       1      59      58      74                        
ENDDATA
