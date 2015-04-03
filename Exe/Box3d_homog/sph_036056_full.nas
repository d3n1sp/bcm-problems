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
$   Date       : Mon Apr 03 09:58:05 2006
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
$ FEMAP Property 2 : Untitled
PSOLID         2       1       0        
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           5       0 0.29183 0.29183 0.29183       0        
GRID           7       0 0.75496     0.5 0.24504       0        
GRID           9       0      0.      0.      0.       0        
GRID          10       0      0.     0.5     0.5       0        
GRID          13       0     0.5      0.     0.5       0        
GRID          17       0     0.5 0.24504 0.24504       0        
GRID          19       0     0.5 0.24504 0.75496       0        
GRID          20       0 0.29183 0.29183 0.70817       0        
GRID          21       0 0.70817 0.29183 0.29183       0        
GRID          22       0     0.5 0.13944     0.5       0        
GRID          23       0 0.24504 0.24504     0.5       0        
GRID          24       0      0.      0.      1.       0        
GRID          26       0     0.5     0.5      1.       0        
GRID          28       0      1.      0.      0.       0        
GRID          29       0      0.      0.     0.5       0        
GRID          30       0     0.5      0.      1.       0        
GRID          32       0     0.5     0.5      0.       0        
GRID          33       0      1.     0.5      0.       0        
GRID          35       0     0.5      0.      0.       0        
GRID          37       0     0.5     0.5 0.86056       0        
GRID          39       0 0.70817 0.29183 0.70817       0        
GRID          41       0 0.75496 0.24504     0.5       0        
GRID          42       0      1.      0.      1.       0        
GRID          45       0      1.      0.     0.5       0        
GRID          47       0 0.24504     0.5 0.24504       0        
GRID          48       0 0.13944     0.5     0.5       0        
GRID          49       0 0.29183 0.70817 0.29183       0        
GRID          51       0 0.86056     0.5     0.5       0        
GRID          52       0      0.      1.      0.       0        
GRID          54       0      0.     0.5      0.       0        
GRID          55       0     0.5      1.     0.5       0        
GRID          56       0     0.5     0.5 0.13944       0        
GRID          57       0 0.24504     0.5 0.75496       0        
GRID          58       0     0.5 0.75496 0.24504       0        
GRID          59       0     0.5 0.75496 0.75496       0        
GRID          60       0 0.29183 0.70817 0.70817       0        
GRID          61       0 0.70817 0.70817 0.29183       0        
GRID          62       0     0.5 0.86056     0.5       0        
GRID          63       0 0.24504 0.75496     0.5       0        
GRID          64       0      0.      1.      1.       0        
GRID          65       0     0.5     0.5     0.5       0        
GRID          67       0      0.     0.5      1.       0        
GRID          68       0      1.      1.      0.       0        
GRID          69       0      0.      1.     0.5       0        
GRID          70       0     0.5      1.      1.       0        
GRID          71       0      1.     0.5     0.5       0        
GRID          74       0     0.5      1.      0.       0        
GRID          75       0 0.75496     0.5 0.75496       0        
GRID          77       0 0.70817 0.70817 0.70817       0        
GRID          78       0 0.75496 0.75496     0.5       0        
GRID          79       0      1.      1.      1.       0        
GRID          80       0      1.     0.5      1.       0        
GRID          81       0      1.      1.     0.5       0        
CTETRA         1       2       5      48      65      47                        
CTETRA         2       2       5      65      48      23                        
CTETRA         3       2       5      56      65      17                        
CTETRA         4       2       5      65      56      47                        
CTETRA         5       2       5      22      65      23                        
CTETRA         6       2       5      65      22      17                        
CPENTA         7       1       5      48      47       9      10      54        
CPENTA         8       1       9      10      29       5      48      23        
CPENTA         9       1       5      56      17       9      32      35        
CPENTA        10       1       9      32      54       5      56      47        
CPENTA        11       1       5      22      23       9      13      29        
CPENTA        12       1       9      13      35       5      22      17        
CTETRA        13       2      65      48      20      57                        
CTETRA        14       2      48      65      20      23                        
CTETRA        15       2      65      37      20      19                        
CTETRA        16       2      37      65      20      57                        
CTETRA        17       2      65      22      20      23                        
CTETRA        18       2      22      65      20      19                        
CPENTA        19       1      57      48      20      67      10      24        
CPENTA        20       1      29      10      24      23      48      20        
CPENTA        21       1      19      37      20      30      26      24        
CPENTA        22       1      67      26      24      57      37      20        
CPENTA        23       1      23      22      20      29      13      24        
CPENTA        24       1      30      13      24      19      22      20        
CTETRA        25       2      65      51      21       7                        
CTETRA        26       2      51      65      21      41                        
CTETRA        27       2      65      56      21      17                        
CTETRA        28       2      56      65      21       7                        
CTETRA        29       2      65      22      21      41                        
CTETRA        30       2      22      65      21      17                        
CPENTA        31       1       7      51      21      33      71      28        
CPENTA        32       1      45      71      28      41      51      21        
CPENTA        33       1      17      56      21      35      32      28        
CPENTA        34       1      33      32      28       7      56      21        
CPENTA        35       1      41      22      21      45      13      28        
CPENTA        36       1      35      13      28      17      22      21        
CTETRA        37       2      39      51      65      75                        
CTETRA        38       2      39      65      51      41                        
CTETRA        39       2      39      37      65      19                        
CTETRA        40       2      39      65      37      75                        
CTETRA        41       2      39      22      65      41                        
CTETRA        42       2      39      65      22      19                        
CPENTA        43       1      39      51      75      42      71      80        
CPENTA        44       1      42      71      45      39      51      41        
CPENTA        45       1      39      37      19      42      26      30        
CPENTA        46       1      42      26      80      39      37      75        
CPENTA        47       1      39      22      41      42      13      45        
CPENTA        48       1      42      13      30      39      22      19        
CTETRA        49       2      65      48      49      47                        
CTETRA        50       2      48      65      49      63                        
CTETRA        51       2      65      56      49      58                        
CTETRA        52       2      56      65      49      47                        
CTETRA        53       2      65      62      49      63                        
CTETRA        54       2      62      65      49      58                        
CPENTA        55       1      47      48      49      54      10      52        
CPENTA        56       1      69      10      52      63      48      49        
CPENTA        57       1      58      56      49      74      32      52        
CPENTA        58       1      54      32      52      47      56      49        
CPENTA        59       1      63      62      49      69      55      52        
CPENTA        60       1      74      55      52      58      62      49        
CTETRA        61       2      60      48      65      57                        
CTETRA        62       2      60      65      48      63                        
CTETRA        63       2      60      37      65      59                        
CTETRA        64       2      60      65      37      57                        
CTETRA        65       2      60      62      65      63                        
CTETRA        66       2      60      65      62      59                        
CPENTA        67       1      60      48      57      64      10      67        
CPENTA        68       1      64      10      69      60      48      63        
CPENTA        69       1      60      37      59      64      26      70        
CPENTA        70       1      64      26      67      60      37      57        
CPENTA        71       1      60      62      63      64      55      69        
CPENTA        72       1      64      55      70      60      62      59        
CTETRA        73       2      61      51      65       7                        
CTETRA        74       2      61      65      51      78                        
CTETRA        75       2      61      56      65      58                        
CTETRA        76       2      61      65      56       7                        
CTETRA        77       2      61      62      65      78                        
CTETRA        78       2      61      65      62      58                        
CPENTA        79       1      61      51       7      68      71      33        
CPENTA        80       1      68      71      81      61      51      78        
CPENTA        81       1      61      56      58      68      32      74        
CPENTA        82       1      68      32      33      61      56       7        
CPENTA        83       1      61      62      78      68      55      81        
CPENTA        84       1      68      55      74      61      62      58        
CTETRA        85       2      65      51      77      75                        
CTETRA        86       2      51      65      77      78                        
CTETRA        87       2      65      37      77      59                        
CTETRA        88       2      37      65      77      75                        
CTETRA        89       2      65      62      77      78                        
CTETRA        90       2      62      65      77      59                        
CPENTA        91       1      75      51      77      80      71      79        
CPENTA        92       1      81      71      79      78      51      77        
CPENTA        93       1      59      37      77      70      26      79        
CPENTA        94       1      80      26      79      75      37      77        
CPENTA        95       1      78      62      77      81      55      79        
CPENTA        96       1      70      55      79      59      62      77        
ENDDATA
