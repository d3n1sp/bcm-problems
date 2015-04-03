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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Parametric\Solid\UnitBox.MOD
$   Date       : Mon Oct 10 15:30:01 2005
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
$ FEMAP Property 1 : UnitPlane
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : MaterialNull
MAT1           1                              0.      0.      0.        
GRID           5       0      0.      1.      0.       0        
GRID           6       0      0.     0.9      0.       0        
GRID           7       0      0.     0.8      0.       0        
GRID           8       0      0.     0.7      0.       0        
GRID           9       0      0.     0.6      0.       0        
GRID          10       0      0.     0.5      0.       0        
GRID          11       0      0.     0.4      0.       0        
GRID          12       0      0.     0.3      0.       0        
GRID          13       0      0.     0.2      0.       0        
GRID          14       0      0.     0.1      0.       0        
GRID          15       0      0.      0.      0.       0        
GRID          16       0     0.1      0.      0.       0        
GRID          17       0     0.2      0.      0.       0        
GRID          18       0     0.3      0.      0.       0        
GRID          19       0     0.4      0.      0.       0        
GRID          20       0     0.5      0.      0.       0        
GRID          21       0     0.6      0.      0.       0        
GRID          22       0     0.7      0.      0.       0        
GRID          23       0     0.8      0.      0.       0        
GRID          24       0     0.9      0.      0.       0        
GRID          25       0      1.      0.      0.       0        
GRID          26       0      1.     0.1      0.       0        
GRID          27       0      1.     0.2      0.       0        
GRID          28       0      1.     0.3      0.       0        
GRID          29       0      1.     0.4      0.       0        
GRID          30       0      1.     0.5      0.       0        
GRID          31       0      1.     0.6      0.       0        
GRID          32       0      1.     0.7      0.       0        
GRID          33       0      1.     0.8      0.       0        
GRID          34       0      1.     0.9      0.       0        
GRID          35       0      1.      1.      0.       0        
GRID          36       0     0.9      1.      0.       0        
GRID          37       0     0.8      1.      0.       0        
GRID          38       0     0.7      1.      0.       0        
GRID          39       0     0.6      1.      0.       0        
GRID          40       0     0.5      1.      0.       0        
GRID          41       0     0.4      1.      0.       0        
GRID          42       0     0.3      1.      0.       0        
GRID          43       0     0.2      1.      0.       0        
GRID          44       0     0.1      1.      0.       0        
GRID          45       0     0.1     0.9      0.       0        
GRID          46       0     0.1     0.8      0.       0        
GRID          47       0     0.1     0.7      0.       0        
GRID          48       0     0.1     0.6      0.       0        
GRID          49       0     0.1     0.5      0.       0        
GRID          50       0     0.1     0.4      0.       0        
GRID          51       0     0.1     0.3      0.       0        
GRID          52       0     0.1     0.2      0.       0        
GRID          53       0     0.1     0.1      0.       0        
GRID          54       0     0.2     0.9      0.       0        
GRID          55       0     0.2     0.8      0.       0        
GRID          56       0     0.2     0.7      0.       0        
GRID          57       0     0.2     0.6      0.       0        
GRID          58       0     0.2     0.5      0.       0        
GRID          59       0     0.2     0.4      0.       0        
GRID          60       0     0.2     0.3      0.       0        
GRID          61       0     0.2     0.2      0.       0        
GRID          62       0     0.2     0.1      0.       0        
GRID          63       0     0.3     0.9      0.       0        
GRID          64       0     0.3     0.8      0.       0        
GRID          65       0     0.3     0.7      0.       0        
GRID          66       0     0.3     0.6      0.       0        
GRID          67       0     0.3     0.5      0.       0        
GRID          68       0     0.3     0.4      0.       0        
GRID          69       0     0.3     0.3      0.       0        
GRID          70       0     0.3     0.2      0.       0        
GRID          71       0     0.3     0.1      0.       0        
GRID          72       0     0.4     0.9      0.       0        
GRID          73       0     0.4     0.8      0.       0        
GRID          74       0     0.4     0.7      0.       0        
GRID          75       0     0.4     0.6      0.       0        
GRID          76       0     0.4     0.5      0.       0        
GRID          77       0     0.4     0.4      0.       0        
GRID          78       0     0.4     0.3      0.       0        
GRID          79       0     0.4     0.2      0.       0        
GRID          80       0     0.4     0.1      0.       0        
GRID          81       0     0.5     0.9      0.       0        
GRID          82       0     0.5     0.8      0.       0        
GRID          83       0     0.5     0.7      0.       0        
GRID          84       0     0.5     0.6      0.       0        
GRID          85       0     0.5     0.5      0.       0        
GRID          86       0     0.5     0.4      0.       0        
GRID          87       0     0.5     0.3      0.       0        
GRID          88       0     0.5     0.2      0.       0        
GRID          89       0     0.5     0.1      0.       0        
GRID          90       0     0.6     0.9      0.       0        
GRID          91       0     0.6     0.8      0.       0        
GRID          92       0     0.6     0.7      0.       0        
GRID          93       0     0.6     0.6      0.       0        
GRID          94       0     0.6     0.5      0.       0        
GRID          95       0     0.6     0.4      0.       0        
GRID          96       0     0.6     0.3      0.       0        
GRID          97       0     0.6     0.2      0.       0        
GRID          98       0     0.6     0.1      0.       0        
GRID          99       0     0.7     0.9      0.       0        
GRID         100       0     0.7     0.8      0.       0        
GRID         101       0     0.7     0.7      0.       0        
GRID         102       0     0.7     0.6      0.       0        
GRID         103       0     0.7     0.5      0.       0        
GRID         104       0     0.7     0.4      0.       0        
GRID         105       0     0.7     0.3      0.       0        
GRID         106       0     0.7     0.2      0.       0        
GRID         107       0     0.7     0.1      0.       0        
GRID         108       0     0.8     0.9      0.       0        
GRID         109       0     0.8     0.8      0.       0        
GRID         110       0     0.8     0.7      0.       0        
GRID         111       0     0.8     0.6      0.       0        
GRID         112       0     0.8     0.5      0.       0        
GRID         113       0     0.8     0.4      0.       0        
GRID         114       0     0.8     0.3      0.       0        
GRID         115       0     0.8     0.2      0.       0        
GRID         116       0     0.8     0.1      0.       0        
GRID         117       0     0.9     0.9      0.       0        
GRID         118       0     0.9     0.8      0.       0        
GRID         119       0     0.9     0.7      0.       0        
GRID         120       0     0.9     0.6      0.       0        
GRID         121       0     0.9     0.5      0.       0        
GRID         122       0     0.9     0.4      0.       0        
GRID         123       0     0.9     0.3      0.       0        
GRID         124       0     0.9     0.2      0.       0        
GRID         125       0     0.9     0.1      0.       0        
CQUAD4         2       1       5       6      45      44                
CQUAD4         3       1       6       7      46      45                
CQUAD4         4       1       7       8      47      46                
CQUAD4         5       1       8       9      48      47                
CQUAD4         6       1       9      10      49      48                
CQUAD4         7       1      10      11      50      49                
CQUAD4         8       1      11      12      51      50                
CQUAD4         9       1      12      13      52      51                
CQUAD4        10       1      13      14      53      52                
CQUAD4        11       1      14      15      16      53                
CQUAD4        12       1      44      45      54      43                
CQUAD4        13       1      45      46      55      54                
CQUAD4        14       1      46      47      56      55                
CQUAD4        15       1      47      48      57      56                
CQUAD4        16       1      48      49      58      57                
CQUAD4        17       1      49      50      59      58                
CQUAD4        18       1      50      51      60      59                
CQUAD4        19       1      51      52      61      60                
CQUAD4        20       1      52      53      62      61                
CQUAD4        21       1      53      16      17      62                
CQUAD4        22       1      43      54      63      42                
CQUAD4        23       1      54      55      64      63                
CQUAD4        24       1      55      56      65      64                
CQUAD4        25       1      56      57      66      65                
CQUAD4        26       1      57      58      67      66                
CQUAD4        27       1      58      59      68      67                
CQUAD4        28       1      59      60      69      68                
CQUAD4        29       1      60      61      70      69                
CQUAD4        30       1      61      62      71      70                
CQUAD4        31       1      62      17      18      71                
CQUAD4        32       1      42      63      72      41                
CQUAD4        33       1      63      64      73      72                
CQUAD4        34       1      64      65      74      73                
CQUAD4        35       1      65      66      75      74                
CQUAD4        36       1      66      67      76      75                
CQUAD4        37       1      67      68      77      76                
CQUAD4        38       1      68      69      78      77                
CQUAD4        39       1      69      70      79      78                
CQUAD4        40       1      70      71      80      79                
CQUAD4        41       1      71      18      19      80                
CQUAD4        42       1      41      72      81      40                
CQUAD4        43       1      72      73      82      81                
CQUAD4        44       1      73      74      83      82                
CQUAD4        45       1      74      75      84      83                
CQUAD4        46       1      75      76      85      84                
CQUAD4        47       1      76      77      86      85                
CQUAD4        48       1      77      78      87      86                
CQUAD4        49       1      78      79      88      87                
CQUAD4        50       1      79      80      89      88                
CQUAD4        51       1      80      19      20      89                
CQUAD4        52       1      40      81      90      39                
CQUAD4        53       1      81      82      91      90                
CQUAD4        54       1      82      83      92      91                
CQUAD4        55       1      83      84      93      92                
CQUAD4        56       1      84      85      94      93                
CQUAD4        57       1      85      86      95      94                
CQUAD4        58       1      86      87      96      95                
CQUAD4        59       1      87      88      97      96                
CQUAD4        60       1      88      89      98      97                
CQUAD4        61       1      89      20      21      98                
CQUAD4        62       1      39      90      99      38                
CQUAD4        63       1      90      91     100      99                
CQUAD4        64       1      91      92     101     100                
CQUAD4        65       1      92      93     102     101                
CQUAD4        66       1      93      94     103     102                
CQUAD4        67       1      94      95     104     103                
CQUAD4        68       1      95      96     105     104                
CQUAD4        69       1      96      97     106     105                
CQUAD4        70       1      97      98     107     106                
CQUAD4        71       1      98      21      22     107                
CQUAD4        72       1      38      99     108      37                
CQUAD4        73       1      99     100     109     108                
CQUAD4        74       1     100     101     110     109                
CQUAD4        75       1     101     102     111     110                
CQUAD4        76       1     102     103     112     111                
CQUAD4        77       1     103     104     113     112                
CQUAD4        78       1     104     105     114     113                
CQUAD4        79       1     105     106     115     114                
CQUAD4        80       1     106     107     116     115                
CQUAD4        81       1     107      22      23     116                
CQUAD4        82       1      37     108     117      36                
CQUAD4        83       1     108     109     118     117                
CQUAD4        84       1     109     110     119     118                
CQUAD4        85       1     110     111     120     119                
CQUAD4        86       1     111     112     121     120                
CQUAD4        87       1     112     113     122     121                
CQUAD4        88       1     113     114     123     122                
CQUAD4        89       1     114     115     124     123                
CQUAD4        90       1     115     116     125     124                
CQUAD4        91       1     116      23      24     125                
CQUAD4        92       1      36     117      34      35                
CQUAD4        93       1     117     118      33      34                
CQUAD4        94       1     118     119      32      33                
CQUAD4        95       1     119     120      31      32                
CQUAD4        96       1     120     121      30      31                
CQUAD4        97       1     121     122      29      30                
CQUAD4        98       1     122     123      28      29                
CQUAD4        99       1     123     124      27      28                
CQUAD4       100       1     124     125      26      27                
CQUAD4       101       1     125      24      25      26                
ENDDATA
