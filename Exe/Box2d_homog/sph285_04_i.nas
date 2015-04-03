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
$   Date       : Thu Jun 14 16:42:11 2007
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
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           2       00.040106 0.42554      0.       0        
GRID           3       0 0.10405 0.52125      0.       0        
GRID           4       0 0.19976 0.58519      0.       0        
GRID           7       0 0.20843  0.6253      0.       0        
GRID           8       0 0.10422  0.6253      0.       0        
GRID           9       0      0.  0.6253      0.       0        
GRID          10       0      0. 0.52108      0.       0        
GRID          11       0      0. 0.41687      0.       0        
GRID          19       0 0.02765 0.31265      0.       0        
GRID          21       0 0.11112 0.51418      0.       0        
GRID          22       0 0.20359 0.57596      0.       0        
GRID          23       0 0.31265 0.31265      0.       0        
GRID          24       0 0.31265 0.45515      0.       0        
GRID          25       0 0.31265 0.59765      0.       0        
GRID          28       00.049344 0.42171      0.       0        
GRID          30       0 0.17015 0.31265      0.       0        
GRID          31       0 0.18393 0.44137      0.       0        
GRID          32       00.040106 0.19976      0.       0        
GRID          33       0 0.10405 0.10405      0.       0        
GRID          34       0 0.199760.040106      0.       0        
GRID          35       0 0.31265      0.      0.       0        
GRID          36       0 0.20843      0.      0.       0        
GRID          37       0 0.10422      0.      0.       0        
GRID          38       0      0.      0.      0.       0        
GRID          39       0      0. 0.10422      0.       0        
GRID          40       0      0. 0.20843      0.       0        
GRID          41       0      0. 0.31265      0.       0        
GRID          43       0 0.31265 0.01765      0.       0        
GRID          44       0 0.01765 0.31265      0.       0        
GRID          46       00.049344 0.20359      0.       0        
GRID          48       0 0.203590.049344      0.       0        
GRID          50       0 0.31265 0.17015      0.       0        
GRID          53       0 0.11112 0.11112      0.       0        
GRID          57       0 0.18393 0.18393      0.       0        
GRID          58       0 0.58519 0.42554      0.       0        
GRID          59       0 0.52125 0.52125      0.       0        
GRID          60       0 0.42554 0.58519      0.       0        
GRID          61       0 0.31265  0.6253      0.       0        
GRID          62       0 0.41687  0.6253      0.       0        
GRID          63       0 0.52108  0.6253      0.       0        
GRID          64       0  0.6253  0.6253      0.       0        
GRID          65       0  0.6253 0.52108      0.       0        
GRID          66       0  0.6253 0.41687      0.       0        
GRID          69       0 0.31265 0.60765      0.       0        
GRID          72       0 0.57596 0.42171      0.       0        
GRID          74       0 0.42171 0.57596      0.       0        
GRID          79       0 0.51418 0.51418      0.       0        
GRID          83       0 0.44137 0.44137      0.       0        
GRID          84       0 0.58519 0.19976      0.       0        
GRID          85       0 0.52125 0.10405      0.       0        
GRID          86       0 0.425540.040106      0.       0        
GRID          88       0 0.41687      0.      0.       0        
GRID          89       0 0.52108      0.      0.       0        
GRID          90       0  0.6253      0.      0.       0        
GRID          91       0  0.6253 0.10422      0.       0        
GRID          92       0  0.6253 0.20843      0.       0        
GRID          93       0  0.6253 0.31265      0.       0        
GRID          94       0 0.31265 0.02765      0.       0        
GRID          96       0 0.60765 0.31265      0.       0        
GRID          97       0 0.59765 0.31265      0.       0        
GRID         104       0 0.421710.049344      0.       0        
GRID         105       0 0.51418 0.11112      0.       0        
GRID         106       0 0.57596 0.20359      0.       0        
GRID         108       0 0.45515 0.31265      0.       0        
GRID         109       0 0.44137 0.18393      0.       0        
CQUAD4         1       1       4      69      61       7                
CQUAD4         2       1       3       4       7       8                
CQUAD4         4       1       2       3      10      11                
CQUAD4         5       1      44       2      11      41                
CQUAD4         6       1      25      69       4      22                
CQUAD4         7       1      22       4       3      21                
CQUAD4         8       1      21       3       2      28                
CQUAD4         9       1      28       2      44      19                
CQUAD4        10       2      28      19      30      31                
CTRIA3        11       2      22      21      31                        
CQUAD4        12       2      24      25      22      31                
CQUAD4        13       2      23      24      31      30                
CTRIA3        14       2      21      28      31                        
CQUAD4        15       1      43      34      36      35                
CQUAD4        16       1      34      33      37      36                
CQUAD4        18       1      33      32      40      39                
CQUAD4        19       1      32      44      41      40                
CQUAD4        20       1      43      94      48      34                
CQUAD4        21       1      34      48      53      33                
CQUAD4        22       1      33      53      46      32                
CQUAD4        23       1      32      46      19      44                
CQUAD4        24       2      19      46      57      30                
CTRIA3        25       2      53      48      57                        
CQUAD4        26       2      94      50      57      48                
CQUAD4        27       2      50      23      30      57                
CTRIA3        28       2      46      53      57                        
CQUAD4        29       1      69      60      62      61                
CQUAD4        30       1      60      59      63      62                
CQUAD4        32       1      59      58      66      65                
CQUAD4        33       1      58      96      93      66                
CQUAD4        34       1      69      25      74      60                
CQUAD4        35       1      60      74      79      59                
CQUAD4        36       1      59      79      72      58                
CQUAD4        37       1      58      72      97      96                
CQUAD4        38       2      97      72      83     108                
CTRIA3        39       2      79      74      83                        
CQUAD4        40       2      25      24      83      74                
CQUAD4        41       2      24      23     108      83                
CTRIA3        42       2      72      79      83                        
CQUAD4        43       1      86      43      35      88                
CQUAD4        44       1      85      86      88      89                
CQUAD4        46       1      84      85      91      92                
CQUAD4        47       1      96      84      92      93                
CQUAD4        48       1      94      43      86     104                
CQUAD4        49       1     104      86      85     105                
CQUAD4        50       1     105      85      84     106                
CQUAD4        51       1     106      84      96      97                
CQUAD4        52       2     106      97     108     109                
CTRIA3        53       2     104     105     109                        
CQUAD4        54       2      50      94     104     109                
CQUAD4        55       2      23      50     109     108                
CTRIA3        56       2     105     106     109                        
CTRIA3        57       1       9       3       8                        
CTRIA3        58       1       3       9      10                        
CTRIA3        59       1      59      65      64                        
CTRIA3        60       1      59      64      63                        
CTRIA3        61       1      85      90      91                        
CTRIA3        62       1      90      85      89                        
CTRIA3        63       1      37      33      38                        
CTRIA3        64       1      38      33      39                        
ENDDATA
