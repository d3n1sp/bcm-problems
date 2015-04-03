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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box3d_homog\Solid\model4.MOD
$   Date       : Sun Feb 12 21:15:39 2006
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
$ FEMAP Property 2 : Fulleren
PSOLID         2       1       0        
$ FEMAP Material 1 : Mat
MAT1           1                              0.      0.      0.        
GRID           1       0 0.24504     0.5 0.24504       0        
GRID           2       0 0.13944     0.5     0.5       0        
GRID           3       0     0.5     0.5 0.13944       0        
GRID           4       0     0.5 0.24504 0.24504       0        
GRID           5       0 0.29183 0.29183 0.29183       0        
GRID           6       0     0.5     0.5     0.5       0        
GRID           7       0     0.5 0.13944     0.5       0        
GRID           8       0 0.24504 0.24504     0.5       0        
GRID           9       0      0.      0.      0.       0        
GRID          10       0      0.     0.5     0.5       0        
GRID          11       0     0.5     0.5      0.       0        
GRID          12       0      0.     0.5      0.       0        
GRID          13       0     0.5      0.     0.5       0        
GRID          14       0      0.      0.     0.5       0        
GRID          15       0     0.5      0.      0.       0        
CTETRA         1       2       5       2       6       1                        
CTETRA         2       2       5       6       2       8                        
CTETRA         3       2       5       3       6       4                        
CTETRA         4       2       5       6       3       1                        
CTETRA         5       2       5       7       6       8                        
CTETRA         6       2       5       6       7       4                        
CPENTA         7       1       5       2       1       9      10      12        
CPENTA         8       1       9      10      14       5       2       8        
CPENTA         9       1       5       3       4       9      11      15        
CPENTA        10       1       9      11      12       5       3       1        
CPENTA        11       1       5       7       8       9      13      14        
CPENTA        12       1       9      13      15       5       7       4        
ENDDATA
