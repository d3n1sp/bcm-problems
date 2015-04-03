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
$   Date       : Fri Apr 09 10:29:01 2010
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
$ FEMAP Property 1 : Untitled
PSOLID         1       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0      0.      0.      1.       0        
GRID           2       0    0.25    0.25      1.       0        
GRID           3       0    0.25      0.      1.       0        
GRID           4       0      0.      0.    0.75       0        
GRID           5       0    0.25      0.    0.75       0        
GRID           6       0      0.      0.     0.5       0        
GRID           7       0      0.    0.25      1.       0        
GRID           8       0      0.    0.25    0.75       0        
GRID           9       0    0.25    0.25     0.5       0        
GRID          10       0     0.5     0.5      1.       0        
GRID          11       0     0.5    0.25      1.       0        
GRID          12       0     0.5      0.    0.75       0        
GRID          13       0    0.25     0.5    0.75       0        
GRID          14       0     0.5      0.      0.       0        
GRID          15       0    0.25      0.      0.       0        
GRID          16       0      0.      0.      0.       0        
GRID          17       0    0.25    0.25      0.       0        
GRID          18       0      0.    0.25      0.       0        
GRID          19       0      0.      0.    0.25       0        
GRID          20       0      0.    0.25     0.5       0        
GRID          21       0      0.    0.25    0.25       0        
GRID          22       0    0.25      0.     0.5       0        
GRID          23       0    0.25      0.    0.25       0        
GRID          24       0     0.5    0.25      0.       0        
GRID          25       0      0.     0.5    0.25       0        
GRID          26       0    0.25     0.5     0.5       0        
GRID          27       0     0.5     0.5    0.25       0        
GRID          28       0    0.25     0.5    0.25       0        
GRID          29       0      1.     0.5      1.       0        
GRID          30       0      1.    0.25      1.       0        
GRID          31       0      1.      0.      1.       0        
GRID          32       0      1.      0.    0.75       0        
GRID          33       0      1.      0.     0.5       0        
GRID          34       0      1.    0.25    0.75       0        
GRID          35       0     0.5      0.      1.       0        
GRID          36       0    0.75      0.      1.       0        
GRID          37       0    0.75    0.25      1.       0        
GRID          38       0     0.5      0.     0.5       0        
GRID          39       0    0.75      0.    0.75       0        
GRID          40       0      1.    0.25     0.5       0        
GRID          41       0    0.75      0.     0.5       0        
GRID          42       0     0.5    0.25    0.75       0        
GRID          43       0      1.     0.5     0.5       0        
GRID          44       0     0.5     0.5     0.5       0        
GRID          45       0     0.5     0.5    0.75       0        
GRID          46       0    0.75     0.5    0.75       0        
GRID          47       0      1.    0.25    0.25       0        
GRID          48       0     0.5    0.25    0.25       0        
GRID          49       0     0.5      0.    0.25       0        
GRID          50       0    0.75      0.      0.       0        
GRID          51       0      1.      0.    0.25       0        
GRID          52       0    0.75      0.    0.25       0        
GRID          53       0      1.     0.5      0.       0        
GRID          54       0      1.    0.25      0.       0        
GRID          55       0      1.      0.      0.       0        
GRID          56       0    0.75    0.25      0.       0        
GRID          57       0    0.75     0.5     0.5       0        
GRID          58       0     0.5    0.25     0.5       0        
GRID          59       0    0.75    0.25     0.5       0        
GRID          60       0    0.75     0.5      0.       0        
GRID          61       0    0.25      1.     0.5       0        
GRID          62       0      0.      1.    0.75       0        
GRID          63       0      0.      1.      1.       0        
GRID          64       0    0.25      1.    0.75       0        
GRID          65       0      0.     0.5      1.       0        
GRID          66       0      0.    0.75    0.75       0        
GRID          67       0    0.25    0.75     0.5       0        
GRID          68       0      0.    0.75      1.       0        
GRID          69       0    0.25      1.      1.       0        
GRID          70       0    0.25    0.75      1.       0        
GRID          71       0     0.5      1.      1.       0        
GRID          72       0      0.     0.5    0.75       0        
GRID          73       0      0.     0.5     0.5       0        
GRID          74       0    0.25     0.5      1.       0        
GRID          75       0     0.5      1.     0.5       0        
GRID          76       0      0.      1.      0.       0        
GRID          77       0      0.      1.    0.25       0        
GRID          78       0    0.25      1.    0.25       0        
GRID          79       0      0.      1.     0.5       0        
GRID          80       0      0.    0.75     0.5       0        
GRID          81       0     0.5    0.75     0.5       0        
GRID          82       0      0.     0.5      0.       0        
GRID          83       0      0.    0.75      0.       0        
GRID          84       0    0.25      1.      0.       0        
GRID          85       0    0.25    0.75      0.       0        
GRID          86       0      0.    0.75    0.25       0        
GRID          87       0     0.5    0.75    0.25       0        
GRID          88       0    0.25     0.5      0.       0        
GRID          89       0     0.5     0.5      0.       0        
GRID          90       0     0.5      1.    0.75       0        
GRID          91       0    0.75      1.      1.       0        
GRID          92       0    0.75      1.    0.75       0        
GRID          93       0      1.    0.75      1.       0        
GRID          94       0      1.     0.5    0.75       0        
GRID          95       0      1.      1.    0.75       0        
GRID          96       0      1.    0.75    0.75       0        
GRID          97       0     0.5    0.75      1.       0        
GRID          98       0     0.5    0.75    0.75       0        
GRID          99       0    0.75     0.5      1.       0        
GRID         100       0      1.      1.      1.       0        
GRID         101       0    0.75    0.75      1.       0        
GRID         102       0     0.5      1.      0.       0        
GRID         103       0    0.75      1.     0.5       0        
GRID         104       0      1.      1.    0.25       0        
GRID         105       0    0.75      1.    0.25       0        
GRID         106       0    0.75     0.5    0.25       0        
GRID         107       0      1.    0.75      0.       0        
GRID         108       0    0.75      1.      0.       0        
GRID         109       0    0.75    0.75      0.       0        
GRID         110       0      1.      1.     0.5       0        
GRID         111       0    0.75    0.75     0.5       0        
GRID         112       0      1.    0.75     0.5       0        
GRID         113       0      1.     0.5    0.25       0        
GRID         114       0      1.      1.      0.       0        
GRID         115       0      1.    0.75    0.25       0        
GRID         116       0     0.5    0.75      0.       0        
GRID         117       0     0.5      1.    0.25       0        
GRID         118       0    0.25    0.25    0.75       0        
GRID         119       0    0.25    0.25    0.25       0        
GRID         120       0    0.75    0.25    0.75       0        
GRID         121       0    0.75    0.25    0.25       0        
GRID         122       0    0.25    0.75    0.75       0        
GRID         123       0    0.25    0.75    0.25       0        
GRID         124       0    0.75    0.75    0.75       0        
GRID         125       0    0.75    0.75    0.25       0        
CHEXA          1       1      38      22       5      12      58       9+EL    1
+EL    1     118      42                                                        
CHEXA          2       1      12       5       3      35      42     118+EL    2
+EL    2       2      11                                                        
CHEXA          3       1      22       6       4       5       9      20+EL    3
+EL    3       8     118                                                        
CHEXA          4       1       5       4       1       3     118       8+EL    4
+EL    4       7       2                                                        
CHEXA          5       1      58       9     118      42      44      26+EL    5
+EL    5      13      45                                                        
CHEXA          6       1      42     118       2      11      45      13+EL    6
+EL    6      74      10                                                        
CHEXA          7       1       9      20       8     118      26      73+EL    7
+EL    7      72      13                                                        
CHEXA          8       1     118       8       7       2      13      72+EL    8
+EL    8      65      74                                                        
CHEXA          9       1      14      15      23      49      24      17+EL    9
+EL    9     119      48                                                        
CHEXA         10       1      49      23      22      38      48     119+EL    A
+EL    A       9      58                                                        
CHEXA         11       1      15      16      19      23      17      18+EL    B
+EL    B      21     119                                                        
CHEXA         12       1      23      19       6      22     119      21+EL    C
+EL    C      20       9                                                        
CHEXA         13       1      24      17     119      48      89      88+EL    D
+EL    D      28      27                                                        
CHEXA         14       1      48     119       9      58      27      28+EL    E
+EL    E      26      44                                                        
CHEXA         15       1      17      18      21     119      88      82+EL    F
+EL    F      25      28                                                        
CHEXA         16       1     119      21      20       9      28      25+EL    G
+EL    G      73      26                                                        
CHEXA         17       1      29      30      37      99      94      34+EL    H
+EL    H     120      46                                                        
CHEXA         18       1      99      37      11      10      46     120+EL    I
+EL    I      42      45                                                        
CHEXA         19       1      30      31      36      37      34      32+EL    J
+EL    J      39     120                                                        
CHEXA         20       1      37      36      35      11     120      39+EL    K
+EL    K      12      42                                                        
CHEXA         21       1      94      34     120      46      43      40+EL    L
+EL    L      59      57                                                        
CHEXA         22       1      46     120      42      45      57      59+EL    M
+EL    M      58      44                                                        
CHEXA         23       1      34      32      39     120      40      33+EL    N
+EL    N      41      59                                                        
CHEXA         24       1     120      39      12      42      59      41+EL    O
+EL    O      38      58                                                        
CHEXA         25       1      44      58      48      27      57      59+EL    P
+EL    P     121     106                                                        
CHEXA         26       1      27      48      24      89     106     121+EL    Q
+EL    Q      56      60                                                        
CHEXA         27       1      58      38      49      48      59      41+EL    R
+EL    R      52     121                                                        
CHEXA         28       1      48      49      14      24     121      52+EL    S
+EL    S      50      56                                                        
CHEXA         29       1      57      59     121     106      43      40+EL    T
+EL    T      47     113                                                        
CHEXA         30       1     106     121      56      60     113      47+EL    U
+EL    U      54      53                                                        
CHEXA         31       1      59      41      52     121      40      33+EL    V
+EL    V      51      47                                                        
CHEXA         32       1     121      52      50      56      47      51+EL    W
+EL    W      55      54                                                        
CHEXA         33       1      79      62      66      80      61      64+EL    X
+EL    X     122      67                                                        
CHEXA         34       1      80      66      72      73      67     122+EL    Y
+EL    Y      13      26                                                        
CHEXA         35       1      62      63      68      66      64      69+EL    Z
+EL    Z      70     122                                                        
CHEXA         36       1      66      68      65      72     122      70+EL   10
+EL   10      74      13                                                        
CHEXA         37       1      61      64     122      67      75      90+EL   11
+EL   11      98      81                                                        
CHEXA         38       1      67     122      13      26      81      98+EL   12
+EL   12      45      44                                                        
CHEXA         39       1      64      69      70     122      90      71+EL   13
+EL   13      97      98                                                        
CHEXA         40       1     122      70      74      13      98      97+EL   14
+EL   14      10      45                                                        
CHEXA         41       1      79      61      67      80      77      78+EL   15
+EL   15     123      86                                                        
CHEXA         42       1      80      67      26      73      86     123+EL   16
+EL   16      28      25                                                        
CHEXA         43       1      61      75      81      67      78     117+EL   17
+EL   17      87     123                                                        
CHEXA         44       1      67      81      44      26     123      87+EL   18
+EL   18      27      28                                                        
CHEXA         45       1      77      78     123      86      76      84+EL   19
+EL   19      85      83                                                        
CHEXA         46       1      86     123      28      25      83      85+EL   1A
+EL   1A      88      82                                                        
CHEXA         47       1      78     117      87     123      84     102+EL   1B
+EL   1B     116      85                                                        
CHEXA         48       1     123      87      27      28      85     116+EL   1C
+EL   1C      89      88                                                        
CHEXA         49       1      43      57      46      94     112     111+EL   1D
+EL   1D     124      96                                                        
CHEXA         50       1      94      46      99      29      96     124+EL   1E
+EL   1E     101      93                                                        
CHEXA         51       1      57      44      45      46     111      81+EL   1F
+EL   1F      98     124                                                        
CHEXA         52       1      46      45      10      99     124      98+EL   1G
+EL   1G      97     101                                                        
CHEXA         53       1     112     111     124      96     110     103+EL   1H
+EL   1H      92      95                                                        
CHEXA         54       1      96     124     101      93      95      92+EL   1I
+EL   1I      91     100                                                        
CHEXA         55       1     111      81      98     124     103      75+EL   1J
+EL   1J      90      92                                                        
CHEXA         56       1     124      98      97     101      92      90+EL   1K
+EL   1K      71      91                                                        
CHEXA         57       1      89      27     106      60     116      87+EL   1L
+EL   1L     125     109                                                        
CHEXA         58       1      60     106     113      53     109     125+EL   1M
+EL   1M     115     107                                                        
CHEXA         59       1      27      44      57     106      87      81+EL   1N
+EL   1N     111     125                                                        
CHEXA         60       1     106      57      43     113     125     111+EL   1O
+EL   1O     112     115                                                        
CHEXA         61       1     116      87     125     109     102     117+EL   1P
+EL   1P     105     108                                                        
CHEXA         62       1     109     125     115     107     108     105+EL   1Q
+EL   1Q     104     114                                                        
CHEXA         63       1      87      81     111     125     117      75+EL   1R
+EL   1R     103     105                                                        
CHEXA         64       1     125     111     112     115     105     103+EL   1S
+EL   1S     110     104                                                        
ENDDATA
