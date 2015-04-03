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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Parametric\Solid\Hexa.MOD
$   Date       : Mon Mar 01 23:21:04 2010
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
$ FEMAP Property 1 : Matrix
PSOLID         1       1       0        
$ FEMAP Property 2 : Untitled
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     0.5      0.      0.       0        
GRID           2       0      0.      0.      0.       0        
GRID           3       0     0.5      0.     0.5       0        
GRID           4       0      1.      0.      0.       0        
GRID           5       0      1.     0.5     0.5       0        
GRID           6       0      0.     0.5     0.5       0        
GRID           7       0     0.5     0.5      0.       0        
GRID           8       0      0.      1.      0.       0        
GRID           9       0      0.      1.     0.5       0        
GRID          10       0      0.     0.5      0.       0        
GRID          11       0      1.     0.5      0.       0        
GRID          12       0     0.5      1.      0.       0        
GRID          13       0      1.      1.      0.       0        
GRID          14       0      0.      0.      1.       0        
GRID          15       0      0.     0.5      1.       0        
GRID          16       0      0.      0.     0.5       0        
GRID          17       0     0.5      0.      1.       0        
GRID          18       0      1.      0.      1.       0        
GRID          19       0      1.      0.     0.5       0        
GRID          20       0     0.5     0.5     0.5       0        
GRID          21       0     0.5     0.5      1.       0        
GRID          22       0      0.      1.      1.       0        
GRID          23       0      1.     0.5      1.       0        
GRID          24       0     0.5      1.     0.5       0        
GRID          25       0     0.5      1.      1.       0        
GRID          26       0      1.      1.      1.       0        
GRID          27       0      1.      1.     0.5       0        
CHEXA          1       1      16       6      20       3       2      10+EL    1
+EL    1       7       1                                                        
CHEXA          2       1       3      20       5      19       1       7+EL    2
+EL    2      11       4                                                        
CHEXA          3       1       6       9      24      20      10       8+EL    3
+EL    3      12       7                                                        
CHEXA          4       1      20      24      27       5       7      12+EL    4
+EL    4      13      11                                                        
CHEXA          5       1      14      15      21      17      16       6+EL    5
+EL    5      20       3                                                        
CHEXA          6       1      17      21      23      18       3      20+EL    6
+EL    6       5      19                                                        
CHEXA          7       1      15      22      25      21       6       9+EL    7
+EL    7      24      20                                                        
CHEXA          8       1      21      25      26      23      20      24+EL    8
+EL    8      27       5                                                        
ENDDATA
