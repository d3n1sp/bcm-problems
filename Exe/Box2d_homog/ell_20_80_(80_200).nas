ID FEMAP,FEMAP
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
$   From Model : 
$   Date       : Wed Apr 19 00:00:11 2006
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
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID          28       0     0.8 0.83333      0.       0        
GRID          29       0     0.8 0.66667      0.       0        
GRID          30       0     0.8     0.5      0.       0        
GRID          31       0     0.8 0.33333      0.       0        
GRID          32       0     0.8 0.16667      0.       0        
GRID          33       0     0.8      0.      0.       0        
GRID          34       0     0.6      0.      0.       0        
GRID          35       0     0.4      0.      0.       0        
GRID          36       0     0.2      0.      0.       0        
GRID          37       0      0.      0.      0.       0        
GRID          38       0      0. 0.16667      0.       0        
GRID          39       0      0. 0.33333      0.       0        
GRID          40       0      0.     0.5      0.       0        
GRID          41       0      0. 0.66667      0.       0        
GRID          42       0      0. 0.83333      0.       0        
GRID          43       0      0.      1.      0.       0        
GRID          45       0 0.20979 0.75279      0.       0        
GRID          48       0  0.3382 0.23915      0.       0        
GRID          51       0 0.51756 0.35279      0.       0        
GRID          54       0 0.14019 0.38942      0.       0        
GRID          55       0 0.18057 0.23072      0.       0        
GRID          56       0  0.3037 0.13408      0.       0        
GRID          57       0 0.65981 0.38942      0.       0        
GRID          58       0 0.61943 0.23072      0.       0        
GRID          59       0  0.4963 0.13408      0.       0        
GRID          61       0     0.4      1.      0.       0        
GRID          62       0     0.6      1.      0.       0        
GRID          63       0 0.59021 0.75279      0.       0        
GRID          64       0  0.5618 0.52977      0.       0        
GRID          66       0  0.4618 0.23915      0.       0        
GRID          67       0     0.4     0.2      0.       0        
GRID          69       0 0.28244 0.35279      0.       0        
GRID          70       0  0.2382 0.52977      0.       0        
GRID          72       0 0.39999 0.64597      0.       0        
GRID          73       0     0.4 0.45267      0.       0        
GRID          74       0     0.4 0.30605      0.       0        
GRID          75       0     0.4 0.85848      0.       0        
GRID          76       0     0.8      1.      0.       0        
GRID          77       0     0.8 1.16667      0.       0        
GRID          78       0     0.8 1.33333      0.       0        
GRID          79       0     0.8     1.5      0.       0        
GRID          80       0     0.8 1.66667      0.       0        
GRID          81       0     0.8 1.83333      0.       0        
GRID          82       0     0.8      2.      0.       0        
GRID          83       0     0.6      2.      0.       0        
GRID          84       0     0.4      2.      0.       0        
GRID          85       0     0.2      2.      0.       0        
GRID          86       0      0.      2.      0.       0        
GRID          87       0      0. 1.83333      0.       0        
GRID          88       0      0. 1.66667      0.       0        
GRID          89       0      0.     1.5      0.       0        
GRID          90       0      0. 1.33333      0.       0        
GRID          91       0      0. 1.16667      0.       0        
GRID          93       0     0.2      1.      0.       0        
GRID          94       0 0.20979 1.24721      0.       0        
GRID          95       0  0.3382 1.76085      0.       0        
GRID          96       0 0.51756 1.64721      0.       0        
GRID          97       0 0.14019 1.61058      0.       0        
GRID          98       0 0.18057 1.76928      0.       0        
GRID          99       0  0.3037 1.86592      0.       0        
GRID         100       0 0.65981 1.61058      0.       0        
GRID         101       0 0.61943 1.76928      0.       0        
GRID         102       0  0.4963 1.86592      0.       0        
GRID         105       0 0.59021 1.24721      0.       0        
GRID         106       0  0.5618 1.47023      0.       0        
GRID         107       0  0.4618 1.76085      0.       0        
GRID         108       0     0.4     1.8      0.       0        
GRID         109       0 0.28244 1.64721      0.       0        
GRID         110       0  0.2382 1.47023      0.       0        
GRID         111       0 0.39999 1.35403      0.       0        
GRID         112       0     0.4 1.54733      0.       0        
GRID         113       0     0.4 1.69395      0.       0        
GRID         114       0     0.4 1.14152      0.       0        
CTRIA3        31       1      93      43      42                        
CTRIA3        32       1      45      93      42                        
CTRIA3        33       1      45      42      41                        
CTRIA3        34       1      70      45      41                        
CTRIA3        35       1      70      41      40                        
CTRIA3        36       1      54      70      40                        
CTRIA3        37       1      54      40      39                        
CTRIA3        38       1      55      54      39                        
CTRIA3        39       1      55      39      38                        
CTRIA3        40       1      38      37      36                        
CTRIA3        41       1      55      38      36                        
CTRIA3        42       1      69      70      54                        
CTRIA3        43       1      69      54      55                        
CTRIA3        44       1      48      69      55                        
CTRIA3        45       1      56      48      55                        
CTRIA3        46       1      56      55      36                        
CTRIA3        47       1      67      48      56                        
CTRIA3        48       1      56      36      35                        
CTRIA3        49       1      67      56      35                        
CTRIA3        50       1      57      64      51                        
CTRIA3        51       1      58      57      51                        
CTRIA3        52       1      59      66      67                        
CTRIA3        53       1      59      67      35                        
CTRIA3        54       1      59      35      34                        
CTRIA3        55       1      58      51      66                        
CTRIA3        56       1      58      66      59                        
CTRIA3        57       1      58      59      34                        
CTRIA3        58       1      34      33      32                        
CTRIA3        59       1      58      34      32                        
CTRIA3        60       1      58      32      31                        
CTRIA3        61       1      57      58      31                        
CTRIA3        62       1      57      31      30                        
CTRIA3        63       1      64      57      30                        
CTRIA3        64       1      64      30      29                        
CTRIA3        65       1      63      64      29                        
CTRIA3        66       1      63      29      28                        
CTRIA3        67       1      28      76      62                        
CTRIA3        68       1      63      28      62                        
CTRIA3        69       2      72      45      70                        
CTRIA3        70       2      73      72      70                        
CTRIA3        71       2      73      70      69                        
CTRIA3        72       2      74      73      69                        
CTRIA3        73       2      74      69      48                        
CTRIA3        74       2      74      48      67                        
CTRIA3        75       2      74      67      66                        
CTRIA3        76       2      74      66      51                        
CTRIA3        77       2      73      74      51                        
CTRIA3        78       2      73      51      64                        
CTRIA3        79       2      72      73      64                        
CTRIA3        80       2      72      64      63                        
CTRIA3        81       2      75      72      63                        
CTRIA3        82       2      75      63      62                        
CTRIA3        83       2      75      62      61                        
CTRIA3        84       2      45      72      75                        
CTRIA3        85       2      75      61      93                        
CTRIA3        86       2      45      75      93                        
CTRIA3        87       1      93      91      43                        
CTRIA3        88       1      94      91      93                        
CTRIA3        89       1      94      90      91                        
CTRIA3        90       1     110      90      94                        
CTRIA3        91       1     110      89      90                        
CTRIA3        92       1      97      89     110                        
CTRIA3        93       1      97      88      89                        
CTRIA3        94       1      98      88      97                        
CTRIA3        95       1      98      87      88                        
CTRIA3        96       1      87      85      86                        
CTRIA3        97       1      98      85      87                        
CTRIA3        98       1     109      97     110                        
CTRIA3        99       1     109      98      97                        
CTRIA3       100       1      95      98     109                        
CTRIA3       101       1      99      98      95                        
CTRIA3       102       1      99      85      98                        
CTRIA3       103       1     108      99      95                        
CTRIA3       104       1      99      84      85                        
CTRIA3       105       1     108      84      99                        
CTRIA3       106       1     100      96     106                        
CTRIA3       107       1     101      96     100                        
CTRIA3       108       1     102     108     107                        
CTRIA3       109       1     102      84     108                        
CTRIA3       110       1     102      83      84                        
CTRIA3       111       1     101     107      96                        
CTRIA3       112       1     101     102     107                        
CTRIA3       113       1     101      83     102                        
CTRIA3       114       1      83      81      82                        
CTRIA3       115       1     101      81      83                        
CTRIA3       116       1     101      80      81                        
CTRIA3       117       1     100      80     101                        
CTRIA3       118       1     100      79      80                        
CTRIA3       119       1     106      79     100                        
CTRIA3       120       1     106      78      79                        
CTRIA3       121       1     105      78     106                        
CTRIA3       122       1     105      77      78                        
CTRIA3       123       1      77      62      76                        
CTRIA3       124       1     105      62      77                        
CTRIA3       125       2     111     110      94                        
CTRIA3       126       2     112     110     111                        
CTRIA3       127       2     112     109     110                        
CTRIA3       128       2     113     109     112                        
CTRIA3       129       2     113      95     109                        
CTRIA3       130       2     113     108      95                        
CTRIA3       131       2     113     107     108                        
CTRIA3       132       2     113      96     107                        
CTRIA3       133       2     112      96     113                        
CTRIA3       134       2     112     106      96                        
CTRIA3       135       2     111     106     112                        
CTRIA3       136       2     111     105     106                        
CTRIA3       137       2     114     105     111                        
CTRIA3       138       2     114      62     105                        
CTRIA3       139       2     114      61      62                        
CTRIA3       140       2      94     114     111                        
CTRIA3       141       2     114      93      61                        
CTRIA3       142       2      94      93     114                        
ENDDATA
