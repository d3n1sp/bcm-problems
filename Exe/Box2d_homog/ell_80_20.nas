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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box2d_fulleren\Solid\model.MOD
$   Date       : Mon Apr 17 13:27:32 2006
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 1. 0.4 0.8 0.2 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID          28       0 0.83333     0.8      0.       0        
GRID          29       0 0.66667     0.8      0.       0        
GRID          30       0     0.5     0.8      0.       0        
GRID          31       0 0.33333     0.8      0.       0        
GRID          32       0 0.16667     0.8      0.       0        
GRID          33       0      0.     0.8      0.       0        
GRID          34       0      0.     0.6      0.       0        
GRID          35       0      0.     0.4      0.       0        
GRID          36       0      0.     0.2      0.       0        
GRID          37       0      0.      0.      0.       0        
GRID          38       0 0.16667      0.      0.       0        
GRID          39       0 0.33333      0.      0.       0        
GRID          40       0     0.5      0.      0.       0        
GRID          41       0 0.66667      0.      0.       0        
GRID          42       0 0.83333      0.      0.       0        
GRID          43       0      1.      0.      0.       0        
GRID          45       0 0.75279 0.20979      0.       0        
GRID          48       0 0.23915  0.3382      0.       0        
GRID          51       0 0.35279 0.51756      0.       0        
GRID          54       0 0.38942 0.14019      0.       0        
GRID          55       0 0.23072 0.18057      0.       0        
GRID          56       0 0.13408  0.3037      0.       0        
GRID          57       0 0.38942 0.65981      0.       0        
GRID          58       0 0.23072 0.61943      0.       0        
GRID          59       0 0.13408  0.4963      0.       0        
GRID          61       0      1.     0.4      0.       0        
GRID          62       0      1.     0.6      0.       0        
GRID          63       0 0.75279 0.59021      0.       0        
GRID          64       0 0.52977  0.5618      0.       0        
GRID          66       0 0.23915  0.4618      0.       0        
GRID          67       0     0.2     0.4      0.       0        
GRID          69       0 0.35279 0.28244      0.       0        
GRID          70       0 0.52977  0.2382      0.       0        
GRID          72       0 0.64597 0.39999      0.       0        
GRID          73       0 0.45267     0.4      0.       0        
GRID          74       0 0.30605     0.4      0.       0        
GRID          75       0 0.85848     0.4      0.       0        
GRID          76       0      1.     0.8      0.       0        
GRID          77       0 1.16667     0.8      0.       0        
GRID          78       0 1.33333     0.8      0.       0        
GRID          79       0     1.5     0.8      0.       0        
GRID          80       0 1.66667     0.8      0.       0        
GRID          81       0 1.83333     0.8      0.       0        
GRID          82       0      2.     0.8      0.       0        
GRID          83       0      2.     0.6      0.       0        
GRID          84       0      2.     0.4      0.       0        
GRID          85       0      2.     0.2      0.       0        
GRID          86       0      2.      0.      0.       0        
GRID          87       0 1.83333      0.      0.       0        
GRID          88       0 1.66667      0.      0.       0        
GRID          89       0     1.5      0.      0.       0        
GRID          90       0 1.33333      0.      0.       0        
GRID          91       0 1.16667      0.      0.       0        
GRID          93       0      1.     0.2      0.       0        
GRID          94       0 1.24721 0.20979      0.       0        
GRID          95       0 1.76085  0.3382      0.       0        
GRID          96       0 1.64721 0.51756      0.       0        
GRID          97       0 1.61058 0.14019      0.       0        
GRID          98       0 1.76928 0.18057      0.       0        
GRID          99       0 1.86592  0.3037      0.       0        
GRID         100       0 1.61058 0.65981      0.       0        
GRID         101       0 1.76928 0.61943      0.       0        
GRID         102       0 1.86592  0.4963      0.       0        
GRID         105       0 1.24721 0.59021      0.       0        
GRID         106       0 1.47023  0.5618      0.       0        
GRID         107       0 1.76085  0.4618      0.       0        
GRID         108       0     1.8     0.4      0.       0        
GRID         109       0 1.64721 0.28244      0.       0        
GRID         110       0 1.47023  0.2382      0.       0        
GRID         111       0 1.35403 0.39999      0.       0        
GRID         112       0 1.54733     0.4      0.       0        
GRID         113       0 1.69395     0.4      0.       0        
GRID         114       0 1.14152     0.4      0.       0        
CTRIA3        31       1      42      43      93                        
CTRIA3        32       1      42      93      45                        
CTRIA3        33       1      41      42      45                        
CTRIA3        34       1      41      45      70                        
CTRIA3        35       1      40      41      70                        
CTRIA3        36       1      40      70      54                        
CTRIA3        37       1      39      40      54                        
CTRIA3        38       1      39      54      55                        
CTRIA3        39       1      38      39      55                        
CTRIA3        40       1      36      37      38                        
CTRIA3        41       1      36      38      55                        
CTRIA3        42       1      54      70      69                        
CTRIA3        43       1      55      54      69                        
CTRIA3        44       1      55      69      48                        
CTRIA3        45       1      55      48      56                        
CTRIA3        46       1      36      55      56                        
CTRIA3        47       1      56      48      67                        
CTRIA3        48       1      35      36      56                        
CTRIA3        49       1      35      56      67                        
CTRIA3        50       1      51      64      57                        
CTRIA3        51       1      51      57      58                        
CTRIA3        52       1      67      66      59                        
CTRIA3        53       1      35      67      59                        
CTRIA3        54       1      34      35      59                        
CTRIA3        55       1      66      51      58                        
CTRIA3        56       1      59      66      58                        
CTRIA3        57       1      34      59      58                        
CTRIA3        58       1      32      33      34                        
CTRIA3        59       1      32      34      58                        
CTRIA3        60       1      31      32      58                        
CTRIA3        61       1      31      58      57                        
CTRIA3        62       1      30      31      57                        
CTRIA3        63       1      30      57      64                        
CTRIA3        64       1      29      30      64                        
CTRIA3        65       1      29      64      63                        
CTRIA3        66       1      28      29      63                        
CTRIA3        67       1      62      76      28                        
CTRIA3        68       1      62      28      63                        
CTRIA3        69       2      70      45      72                        
CTRIA3        70       2      70      72      73                        
CTRIA3        71       2      69      70      73                        
CTRIA3        72       2      69      73      74                        
CTRIA3        73       2      48      69      74                        
CTRIA3        74       2      67      48      74                        
CTRIA3        75       2      66      67      74                        
CTRIA3        76       2      51      66      74                        
CTRIA3        77       2      51      74      73                        
CTRIA3        78       2      64      51      73                        
CTRIA3        79       2      64      73      72                        
CTRIA3        80       2      63      64      72                        
CTRIA3        81       2      63      72      75                        
CTRIA3        82       2      62      63      75                        
CTRIA3        83       2      61      62      75                        
CTRIA3        84       2      75      72      45                        
CTRIA3        85       2      93      61      75                        
CTRIA3        86       2      93      75      45                        
CTRIA3        87       1      43      91      93                        
CTRIA3        88       1      93      91      94                        
CTRIA3        89       1      91      90      94                        
CTRIA3        90       1      94      90     110                        
CTRIA3        91       1      90      89     110                        
CTRIA3        92       1     110      89      97                        
CTRIA3        93       1      89      88      97                        
CTRIA3        94       1      97      88      98                        
CTRIA3        95       1      88      87      98                        
CTRIA3        96       1      86      85      87                        
CTRIA3        97       1      87      85      98                        
CTRIA3        98       1     110      97     109                        
CTRIA3        99       1      97      98     109                        
CTRIA3       100       1     109      98      95                        
CTRIA3       101       1      95      98      99                        
CTRIA3       102       1      98      85      99                        
CTRIA3       103       1      95      99     108                        
CTRIA3       104       1      85      84      99                        
CTRIA3       105       1      99      84     108                        
CTRIA3       106       1     106      96     100                        
CTRIA3       107       1     100      96     101                        
CTRIA3       108       1     107     108     102                        
CTRIA3       109       1     108      84     102                        
CTRIA3       110       1      84      83     102                        
CTRIA3       111       1      96     107     101                        
CTRIA3       112       1     107     102     101                        
CTRIA3       113       1     102      83     101                        
CTRIA3       114       1      82      81      83                        
CTRIA3       115       1      83      81     101                        
CTRIA3       116       1      81      80     101                        
CTRIA3       117       1     101      80     100                        
CTRIA3       118       1      80      79     100                        
CTRIA3       119       1     100      79     106                        
CTRIA3       120       1      79      78     106                        
CTRIA3       121       1     106      78     105                        
CTRIA3       122       1      78      77     105                        
CTRIA3       123       1      76      62      77                        
CTRIA3       124       1      77      62     105                        
CTRIA3       125       2      94     110     111                        
CTRIA3       126       2     111     110     112                        
CTRIA3       127       2     110     109     112                        
CTRIA3       128       2     112     109     113                        
CTRIA3       129       2     109      95     113                        
CTRIA3       130       2      95     108     113                        
CTRIA3       131       2     108     107     113                        
CTRIA3       132       2     107      96     113                        
CTRIA3       133       2     113      96     112                        
CTRIA3       134       2      96     106     112                        
CTRIA3       135       2     112     106     111                        
CTRIA3       136       2     106     105     111                        
CTRIA3       137       2     111     105     114                        
CTRIA3       138       2     105      62     114                        
CTRIA3       139       2      62      61     114                        
CTRIA3       140       2     111     114      94                        
CTRIA3       141       2      61      93     114                        
CTRIA3       142       2     114      93      94                        
ENDDATA
