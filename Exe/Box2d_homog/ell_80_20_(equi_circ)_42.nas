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
$   Date       : Fri Jan 26 18:57:26 2007
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
GRID          45       0     0.2   0.544      0.       0        
GRID          46       0 0.30076 0.78724      0.       0        
GRID          48       0 0.78724 0.78724      0.       0        
GRID          49       0   0.888   0.544      0.       0        
GRID          50       0 0.78724 0.30076      0.       0        
GRID          51       0   0.544     0.2      0.       0        
GRID          52       0 0.30076 0.30076      0.       0        
GRID          53       0 0.14319 0.14319      0.       0        
GRID          54       0 0.94481 0.14319      0.       0        
GRID          55       0 0.94481 0.94481      0.       0        
GRID          56       0 0.14319 0.94481      0.       0        
GRID          59       0   0.544   0.888      0.       0        
GRID          65       0 0.62781 0.46019      0.       0        
GRID          66       0 0.46019 0.62781      0.       0        
CTRIA3        33       1      32      33      53                        
CTRIA3        34       1      45      31      32                        
CTRIA3        35       1      52      45      32                        
CTRIA3        36       1      52      32      53                        
CTRIA3        37       1      53      33      34                        
CTRIA3        38       1      52      53      34                        
CTRIA3        39       1      51      52      34                        
CTRIA3        40       1      51      34      35                        
CTRIA3        41       1      36      37      54                        
CTRIA3        42       1      51      35      36                        
CTRIA3        43       1      50      51      36                        
CTRIA3        44       1      50      36      54                        
CTRIA3        45       1      54      37      38                        
CTRIA3        46       1      50      54      38                        
CTRIA3        47       1      49      50      38                        
CTRIA3        48       1      49      38      39                        
CTRIA3        49       1      40      41      55                        
CTRIA3        50       1      49      39      40                        
CTRIA3        51       1      48      49      40                        
CTRIA3        52       1      48      40      55                        
CTRIA3        53       1      55      41      42                        
CTRIA3        54       1      48      55      42                        
CTRIA3        55       1      59      48      42                        
CTRIA3        56       1      59      42      43                        
CTRIA3        57       1      30      31      45                        
CTRIA3        58       1      30      45      46                        
CTRIA3        59       1      30      46      56                        
CTRIA3        60       1      29      30      56                        
CTRIA3        61       1      59      43      44                        
CTRIA3        62       1      46      59      44                        
CTRIA3        63       1      56      46      44                        
CTRIA3        64       1      29      56      44                        
CTRIA3        65       2      51      50      65                        
CTRIA3        66       2      52      51      65                        
CTRIA3        67       2      52      65      66                        
CTRIA3        68       2      45      52      66                        
CTRIA3        69       2      46      45      66                        
CTRIA3        70       2      59      46      66                        
CTRIA3        71       2      48      59      66                        
CTRIA3        72       2      48      66      65                        
CTRIA3        73       2      49      48      65                        
CTRIA3        74       2      49      65      50                        
ENDDATA
