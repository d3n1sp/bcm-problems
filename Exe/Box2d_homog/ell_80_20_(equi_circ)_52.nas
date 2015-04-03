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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box2d_homog\Solid\model.mod
$   Date       : Fri Jan 26 18:45:44 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.544 0.544 0.344 0.344 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           1       0      0.   1.088      0.       0        
GRID           2       0      0.  0.8704      0.       0        
GRID           3       0      0.  0.6528      0.       0        
GRID           4       0      0.  0.4352      0.       0        
GRID           5       0      0.  0.2176      0.       0        
GRID           6       0      0.      0.      0.       0        
GRID           7       0  0.2176      0.      0.       0        
GRID           8       0  0.4352      0.      0.       0        
GRID           9       0  0.6528      0.      0.       0        
GRID          10       0  0.8704      0.      0.       0        
GRID          11       0   1.088      0.      0.       0        
GRID          12       0   1.088  0.2176      0.       0        
GRID          13       0   1.088  0.4352      0.       0        
GRID          14       0   1.088  0.6528      0.       0        
GRID          15       0   1.088  0.8704      0.       0        
GRID          16       0   1.088   1.088      0.       0        
GRID          17       0  0.8704   1.088      0.       0        
GRID          18       0  0.6528   1.088      0.       0        
GRID          19       0  0.4352   1.088      0.       0        
GRID          20       0  0.2176   1.088      0.       0        
GRID          22       0  0.2657  0.7462      0.       0        
GRID          23       0  0.4377 0.87116      0.       0        
GRID          26       0   0.888   0.544      0.       0        
GRID          27       0  0.8223  0.3418      0.       0        
GRID          29       0  0.4377 0.21684      0.       0        
GRID          31       0 0.86197 0.12937      0.       0        
GRID          32       0 0.86197 0.95863      0.       0        
GRID          33       0 0.22603 0.12937      0.       0        
GRID          35       0  0.8223  0.7462      0.       0        
GRID          36       0  0.6503 0.87116      0.       0        
GRID          39       0     0.2   0.544      0.       0        
GRID          40       0  0.2657  0.3418      0.       0        
GRID          42       0  0.6503 0.21684      0.       0        
GRID          44       0 0.64798 0.46845      0.       0        
GRID          45       0  0.4397 0.61977      0.       0        
GRID          46       0 0.44022 0.40129      0.       0        
GRID          47       0 0.64765  0.6868      0.       0        
CTRIA3         1       1       9      10      31                        
CTRIA3         2       1      42       9      31                        
CTRIA3         3       1      27      42      31                        
CTRIA3         4       1      31      10      11                        
CTRIA3         5       1      31      11      12                        
CTRIA3         6       1      27      31      12                        
CTRIA3         7       1      27      12      13                        
CTRIA3         8       1      26      27      13                        
CTRIA3         9       1      26      13      14                        
CTRIA3        10       1      35      26      14                        
CTRIA3        11       1      35      14      15                        
CTRIA3        12       1      16      17      32                        
CTRIA3        13       1      15      16      32                        
CTRIA3        14       1      35      15      32                        
CTRIA3        15       1      36      35      32                        
CTRIA3        16       1      32      17      18                        
CTRIA3        17       1      36      32      18                        
CTRIA3        18       1      36      18      19                        
CTRIA3        19       1      23      36      19                        
CTRIA3        20       1      23      19      20                        
CTRIA3        21       1       8       9      42                        
CTRIA3        22       1       8      42      29                        
CTRIA3        23       1      29      40      33                        
CTRIA3        24       1       8      29      33                        
CTRIA3        25       1       7       8      33                        
CTRIA3        26       1       6       7      33                        
CTRIA3        27       1       5       6      33                        
CTRIA3        28       1       5      33      40                        
CTRIA3        29       1       4       5      40                        
CTRIA3        30       1       4      40      39                        
CTRIA3        31       1       3       4      39                        
CTRIA3        32       1       3      39      22                        
CTRIA3        33       1       2       3      22                        
CTRIA3        34       1      22      23      20                        
CTRIA3        35       1       2      22      20                        
CTRIA3        36       1       1       2      20                        
CTRIA3        37       2      44      45      46                        
CTRIA3        38       2      42      27      44                        
CTRIA3        39       2      42      44      46                        
CTRIA3        40       2      29      42      46                        
CTRIA3        41       2      40      29      46                        
CTRIA3        42       2      39      40      46                        
CTRIA3        43       2      39      46      45                        
CTRIA3        44       2      22      39      45                        
CTRIA3        45       2      23      22      45                        
CTRIA3        46       2      23      45      47                        
CTRIA3        47       2      36      23      47                        
CTRIA3        48       2      35      36      47                        
CTRIA3        49       2      47      45      44                        
CTRIA3        50       2      26      35      47                        
CTRIA3        51       2      26      47      44                        
CTRIA3        52       2      26      44      27                        
ENDDATA
