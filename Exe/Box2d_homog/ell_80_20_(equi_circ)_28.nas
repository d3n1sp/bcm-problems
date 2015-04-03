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
$   Date       : Fri Jan 26 18:55:26 2007
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
GRID          29       0      0.   1.088      0.       0        
GRID          30       0      0.   0.816      0.       0        
GRID          31       0      0.   0.544      0.       0        
GRID          32       0      0.   0.272      0.       0        
GRID          33       0      0.      0.      0.       0        
GRID          34       0   0.272      0.      0.       0        
GRID          35       0   0.544      0.      0.       0        
GRID          36       0   0.816      0.      0.       0        
GRID          37       0   1.088      0.      0.       0        
GRID          38       0   1.088   0.272      0.       0        
GRID          39       0   1.088   0.544      0.       0        
GRID          40       0   1.088   0.816      0.       0        
GRID          41       0   1.088   1.088      0.       0        
GRID          42       0   0.816   1.088      0.       0        
GRID          43       0   0.544   1.088      0.       0        
GRID          44       0   0.272   1.088      0.       0        
GRID          46       0   0.372 0.84191      0.       0        
GRID          47       0   0.716 0.84191      0.       0        
GRID          48       0   0.888   0.544      0.       0        
GRID          50       0   0.372 0.24609      0.       0        
GRID          54       0     0.2   0.544      0.       0        
GRID          56       0   0.716 0.24609      0.       0        
GRID          57       0   0.544   0.544      0.       0        
CTRIA3        33       1      54      31      32                        
CTRIA3        34       1      50      54      32                        
CTRIA3        35       1      32      33      34                        
CTRIA3        36       1      50      32      34                        
CTRIA3        37       1      50      34      35                        
CTRIA3        38       1      56      50      35                        
CTRIA3        39       1      56      35      36                        
CTRIA3        40       1      36      37      38                        
CTRIA3        41       1      56      36      38                        
CTRIA3        42       1      48      56      38                        
CTRIA3        43       1      48      38      39                        
CTRIA3        44       1      48      39      40                        
CTRIA3        45       1      47      48      40                        
CTRIA3        46       1      40      41      42                        
CTRIA3        47       1      47      40      42                        
CTRIA3        48       1      47      42      43                        
CTRIA3        49       1      46      47      43                        
CTRIA3        50       1      46      43      44                        
CTRIA3        51       1      30      31      54                        
CTRIA3        52       1      30      54      46                        
CTRIA3        53       1      30      46      44                        
CTRIA3        54       1      29      30      44                        
CTRIA3        55       2      54      50      57                        
CTRIA3        56       2      46      54      57                        
CTRIA3        57       2      47      46      57                        
CTRIA3        58       2      57      50      56                        
CTRIA3        59       2      48      47      57                        
CTRIA3        60       2      48      57      56                        
ENDDATA
