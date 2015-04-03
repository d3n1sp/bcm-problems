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
$   Date       : Thu Jun 14 14:20:28 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.34355 0.34355 0.285 0.285 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           2       00.052528  0.4641      0.       0        
GRID           4       0   0.223 0.63457      0.       0        
GRID           7       0 0.22903  0.6871      0.       0        
GRID           8       0 0.11452  0.6871      0.       0        
GRID           9       0      0.  0.6871      0.       0        
GRID          10       0      0. 0.57258      0.       0        
GRID          11       0      0. 0.45807      0.       0        
GRID          12       0      0. 0.34355      0.       0        
GRID          14       0 0.34355 0.65855      0.       0        
GRID          16       0 0.12081 0.56629      0.       0        
GRID          20       00.080244 0.45261      0.       0        
GRID          21       0 0.14202 0.54508      0.       0        
GRID          23       0 0.34355 0.34355      0.       0        
GRID          25       0 0.34355 0.62855      0.       0        
GRID          26       0 0.23449 0.60686      0.       0        
GRID          30       0 0.20105 0.34355      0.       0        
GRID          31       0 0.21483 0.47227      0.       0        
GRID          32       0 0.02855 0.34355      0.       0        
GRID          33       00.052528   0.223      0.       0        
GRID          34       0 0.12081 0.12081      0.       0        
GRID          38       0 0.22903      0.      0.       0        
GRID          39       0 0.11452      0.      0.       0        
GRID          40       0      0.      0.      0.       0        
GRID          41       0      0. 0.11452      0.       0        
GRID          42       0      0. 0.22903      0.       0        
GRID          45       0 0.34355 0.02855      0.       0        
GRID          46       0   0.2230.052528      0.       0        
GRID          51       00.080244 0.23449      0.       0        
GRID          57       0 0.234490.080244      0.       0        
GRID          58       0 0.14202 0.14202      0.       0        
GRID          60       0 0.05855 0.34355      0.       0        
GRID          62       0 0.21483 0.21483      0.       0        
GRID          65       0 0.56629 0.56629      0.       0        
GRID          68       0 0.34355  0.6871      0.       0        
GRID          69       0 0.45807  0.6871      0.       0        
GRID          70       0 0.57258  0.6871      0.       0        
GRID          71       0  0.6871  0.6871      0.       0        
GRID          72       0  0.6871 0.57258      0.       0        
GRID          73       0  0.6871 0.45807      0.       0        
GRID          74       0  0.6871 0.34355      0.       0        
GRID          77       0  0.4641 0.63457      0.       0        
GRID          79       0 0.63457  0.4641      0.       0        
GRID          83       0 0.54508 0.54508      0.       0        
GRID          84       0 0.45261 0.60686      0.       0        
GRID          86       0 0.34355 0.48605      0.       0        
GRID          90       0 0.60686 0.45261      0.       0        
GRID          91       0 0.62855 0.34355      0.       0        
GRID          92       0 0.48605 0.34355      0.       0        
GRID          93       0 0.47227 0.47227      0.       0        
GRID          94       0 0.65855 0.34355      0.       0        
GRID          95       0 0.63457   0.223      0.       0        
GRID          99       0 0.34355      0.      0.       0        
GRID         100       0 0.45807      0.      0.       0        
GRID         101       0 0.57258      0.      0.       0        
GRID         102       0  0.6871      0.      0.       0        
GRID         103       0  0.6871 0.11452      0.       0        
GRID         104       0  0.6871 0.22903      0.       0        
GRID         106       0 0.34355 0.05855      0.       0        
GRID         108       0  0.46410.052528      0.       0        
GRID         109       0 0.56629 0.12081      0.       0        
GRID         117       0 0.34355 0.20105      0.       0        
GRID         119       0 0.452610.080244      0.       0        
GRID         120       0 0.54508 0.14202      0.       0        
GRID         121       0 0.60686 0.23449      0.       0        
GRID         124       0 0.47227 0.21483      0.       0        
CQUAD4         1       1       4      14      68       7                
CQUAD4         2       1      16       4       7       8                
CQUAD4         3       1      16       8       9      10                
CQUAD4         4       1       2      16      10      11                
CQUAD4         5       1      32       2      11      12                
CQUAD4         6       1      25      14       4      26                
CQUAD4         7       1      26       4      16      21                
CQUAD4         8       1      21      16       2      20                
CQUAD4         9       1      20       2      32      60                
CTRIA3        10       2      26      21      31                        
CQUAD4        11       2      86      25      26      31                
CQUAD4        12       2      31      20      60      30                
CQUAD4        13       2      23      86      31      30                
CTRIA3        14       2      21      20      31                        
CQUAD4        15       1      45      46      38      99                
CQUAD4        16       1      46      34      39      38                
CQUAD4        17       1      39      34      41      40                
CQUAD4        18       1      34      33      42      41                
CQUAD4        19       1      33      32      12      42                
CQUAD4        20       1      45     106      57      46                
CQUAD4        21       1      46      57      58      34                
CQUAD4        22       1      34      58      51      33                
CQUAD4        23       1      33      51      60      32                
CTRIA3        24       2      58      57      62                        
CQUAD4        25       2     106     117      62      57                
CQUAD4        26       2      51      62      30      60                
CQUAD4        27       2     117      23      30      62                
CTRIA3        28       2      51      58      62                        
CQUAD4        29       1      14      77      69      68                
CQUAD4        30       1      77      65      70      69                
CQUAD4        31       1      70      65      72      71                
CQUAD4        32       1      65      79      73      72                
CQUAD4        33       1      79      94      74      73                
CQUAD4        34       1      14      25      84      77                
CQUAD4        35       1      77      84      83      65                
CQUAD4        36       1      65      83      90      79                
CQUAD4        37       1      79      90      91      94                
CTRIA3        38       2      83      84      93                        
CQUAD4        39       2      25      86      93      84                
CQUAD4        40       2      90      93      92      91                
CQUAD4        41       2      86      23      92      93                
CTRIA3        42       2      90      83      93                        
CQUAD4        43       1     108      45      99     100                
CQUAD4        44       1     109     108     100     101                
CQUAD4        45       1     109     101     102     103                
CQUAD4        46       1      95     109     103     104                
CQUAD4        47       1      94      95     104      74                
CQUAD4        48       1     106      45     108     119                
CQUAD4        49       1     119     108     109     120                
CQUAD4        50       1     120     109      95     121                
CQUAD4        51       1     121      95      94      91                
CTRIA3        52       2     119     120     124                        
CQUAD4        53       2     117     106     119     124                
CQUAD4        54       2     124     121      91      92                
CQUAD4        55       2      23     117     124      92                
CTRIA3        56       2     120     121     124                        
ENDDATA
