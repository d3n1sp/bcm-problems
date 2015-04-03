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
$   Date       : Thu Jun 14 19:35:02 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.689 0.689 0.623 0.623 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           2       00.094945 0.93507      0.       0        
GRID           3       0 0.23433 1.14367      0.       0        
GRID           4       0 0.44293 1.28305      0.       0        
GRID           7       0 0.45933   1.378      0.       0        
GRID           8       0 0.22967   1.378      0.       0        
GRID           9       0      0.   1.378      0.       0        
GRID          10       0      0. 1.14833      0.       0        
GRID          11       0      0. 0.91867      0.       0        
GRID          18       0   0.046   0.689      0.       0        
GRID          19       0   0.066   0.689      0.       0        
GRID          21       0 0.24847 1.12953      0.       0        
GRID          22       0 0.45059 1.26458      0.       0        
GRID          23       0   0.689   0.689      0.       0        
GRID          24       0   0.689  1.0005      0.       0        
GRID          25       0   0.689   1.312      0.       0        
GRID          28       0 0.11342 0.92741      0.       0        
GRID          30       0  0.3775   0.689      0.       0        
GRID          31       0 0.40763 0.97037      0.       0        
GRID          32       00.094945 0.44293      0.       0        
GRID          33       0 0.23433 0.23433      0.       0        
GRID          34       0 0.442930.094945      0.       0        
GRID          36       0 0.45933      0.      0.       0        
GRID          37       0 0.22967      0.      0.       0        
GRID          38       0      0.      0.      0.       0        
GRID          39       0      0. 0.22967      0.       0        
GRID          40       0      0. 0.45933      0.       0        
GRID          41       0      0.   0.689      0.       0        
GRID          43       0   0.689   0.046      0.       0        
GRID          46       0 0.11342 0.45059      0.       0        
GRID          48       0 0.45059 0.11342      0.       0        
GRID          50       0   0.689  0.3775      0.       0        
GRID          53       0 0.24847 0.24847      0.       0        
GRID          57       0 0.40763 0.40763      0.       0        
GRID          58       0 1.28305 0.93507      0.       0        
GRID          59       0 1.14367 1.14367      0.       0        
GRID          60       0 0.93507 1.28305      0.       0        
GRID          61       0   0.689   1.378      0.       0        
GRID          62       0 0.91867   1.378      0.       0        
GRID          63       0 1.14833   1.378      0.       0        
GRID          64       0   1.378   1.378      0.       0        
GRID          65       0   1.378 1.14833      0.       0        
GRID          66       0   1.378 0.91867      0.       0        
GRID          69       0   0.689   1.332      0.       0        
GRID          74       0 0.92741 1.26458      0.       0        
GRID          79       0 1.12953 1.12953      0.       0        
GRID          80       0 1.26458 0.92741      0.       0        
GRID          82       0  1.0005   0.689      0.       0        
GRID          83       0 0.97037 0.97037      0.       0        
GRID          84       0 1.28305 0.44293      0.       0        
GRID          85       0 1.14367 0.23433      0.       0        
GRID          86       0 0.935070.094945      0.       0        
GRID          87       0   0.689      0.      0.       0        
GRID          88       0 0.91867      0.      0.       0        
GRID          89       0 1.14833      0.      0.       0        
GRID          90       0   1.378      0.      0.       0        
GRID          91       0   1.378 0.22967      0.       0        
GRID          92       0   1.378 0.45933      0.       0        
GRID          93       0   1.378   0.689      0.       0        
GRID          94       0   0.689   0.066      0.       0        
GRID          96       0   1.332   0.689      0.       0        
GRID          97       0   1.312   0.689      0.       0        
GRID         104       0 0.92741 0.11342      0.       0        
GRID         105       0 1.12953 0.24847      0.       0        
GRID         106       0 1.26458 0.45059      0.       0        
GRID         109       0 0.97037 0.40763      0.       0        
CQUAD4         1       1       4      69      61       7                
CQUAD4         2       1       3       4       7       8                
CQUAD4         3       1       3       8       9      10                
CQUAD4         4       1       2       3      10      11                
CQUAD4         5       1      18       2      11      41                
CQUAD4         6       1      25      69       4      22                
CQUAD4         7       1      22       4       3      21                
CQUAD4         8       1      21       3       2      28                
CQUAD4         9       1      28       2      18      19                
CTRIA3        10       2      22      21      31                        
CQUAD4        11       2      24      25      22      31                
CQUAD4        12       2      31      28      19      30                
CQUAD4        13       2      23      24      31      30                
CTRIA3        14       2      21      28      31                        
CQUAD4        15       1      43      34      36      87                
CQUAD4        16       1      34      33      37      36                
CQUAD4        17       1      37      33      39      38                
CQUAD4        18       1      33      32      40      39                
CQUAD4        19       1      32      18      41      40                
CQUAD4        20       1      43      94      48      34                
CQUAD4        21       1      34      48      53      33                
CQUAD4        22       1      33      53      46      32                
CQUAD4        23       1      32      46      19      18                
CTRIA3        24       2      53      48      57                        
CQUAD4        25       2      94      50      57      48                
CQUAD4        26       2      46      57      30      19                
CQUAD4        27       2      50      23      30      57                
CTRIA3        28       2      46      53      57                        
CQUAD4        29       1      69      60      62      61                
CQUAD4        30       1      60      59      63      62                
CQUAD4        31       1      63      59      65      64                
CQUAD4        32       1      59      58      66      65                
CQUAD4        33       1      58      96      93      66                
CQUAD4        34       1      69      25      74      60                
CQUAD4        35       1      60      74      79      59                
CQUAD4        36       1      59      79      80      58                
CQUAD4        37       1      58      80      97      96                
CTRIA3        38       2      79      74      83                        
CQUAD4        39       2      25      24      83      74                
CQUAD4        40       2      80      83      82      97                
CQUAD4        41       2      24      23      82      83                
CTRIA3        42       2      80      79      83                        
CQUAD4        43       1      86      43      87      88                
CQUAD4        44       1      85      86      88      89                
CQUAD4        45       1      85      89      90      91                
CQUAD4        46       1      84      85      91      92                
CQUAD4        47       1      96      84      92      93                
CQUAD4        48       1      94      43      86     104                
CQUAD4        49       1     104      86      85     105                
CQUAD4        50       1     105      85      84     106                
CQUAD4        51       1     106      84      96      97                
CTRIA3        52       2     104     105     109                        
CQUAD4        53       2      50      94     104     109                
CQUAD4        54       2     109     106      97      82                
CQUAD4        55       2      23      50     109      82                
CTRIA3        56       2     105     106     109                        
ENDDATA
