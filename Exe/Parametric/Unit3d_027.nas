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
$   Date       : Fri Apr 09 10:38:25 2010
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
GRID           1       0 0.33333      0.      1.       0        
GRID           2       0 0.66667      0.      1.       0        
GRID           3       0      1. 0.33333      1.       0        
GRID           4       0      1.      1.      1.       0        
GRID           5       0 0.33333      1.      1.       0        
GRID           6       0      0. 0.66667      1.       0        
GRID           7       0 0.33333 0.33333      1.       0        
GRID           8       0 0.66667 0.33333      1.       0        
GRID           9       0 0.33333 0.66667      1.       0        
GRID          10       0 0.66667 0.66667      1.       0        
GRID          11       0      1.      0.      1.       0        
GRID          12       0 0.33333      0.      0.       0        
GRID          13       0 0.66667      0.      0.       0        
GRID          14       0 0.66667      0. 0.33333       0        
GRID          15       0 0.66667      0. 0.66667       0        
GRID          16       0 0.33333      0. 0.33333       0        
GRID          17       0 0.33333      0. 0.66667       0        
GRID          18       0      0.      0. 0.33333       0        
GRID          19       0      0.      0. 0.66667       0        
GRID          20       0      0.      0.      1.       0        
GRID          21       0      0. 0.33333      1.       0        
GRID          22       0      0.      1.      1.       0        
GRID          23       0      0.      1. 0.66667       0        
GRID          24       0      0.      1.      0.       0        
GRID          25       0      0. 0.66667      0.       0        
GRID          26       0      0. 0.33333 0.33333       0        
GRID          27       0      0. 0.33333 0.66667       0        
GRID          28       0      0. 0.66667 0.33333       0        
GRID          29       0      0. 0.66667 0.66667       0        
GRID          30       0      0.      1. 0.33333       0        
GRID          31       0 0.66667      1.      1.       0        
GRID          32       0      1.      1. 0.66667       0        
GRID          33       0 0.33333      1. 0.33333       0        
GRID          34       0 0.33333      1. 0.66667       0        
GRID          35       0 0.66667      1. 0.33333       0        
GRID          36       0 0.66667      1. 0.66667       0        
GRID          37       0      1. 0.66667      0.       0        
GRID          38       0      1. 0.33333      0.       0        
GRID          39       0      0.      0.      0.       0        
GRID          40       0      0. 0.33333      0.       0        
GRID          41       0 0.33333      1.      0.       0        
GRID          42       0 0.66667      1.      0.       0        
GRID          43       0 0.66667 0.66667      0.       0        
GRID          44       0 0.66667 0.33333      0.       0        
GRID          45       0 0.33333 0.66667      0.       0        
GRID          46       0 0.33333 0.33333      0.       0        
GRID          47       0      1.      0.      0.       0        
GRID          48       0      1.      1.      0.       0        
GRID          49       0      1.      1. 0.33333       0        
GRID          50       0      1. 0.66667      1.       0        
GRID          51       0      1.      0. 0.66667       0        
GRID          52       0      1.      0. 0.33333       0        
GRID          53       0      1. 0.33333 0.33333       0        
GRID          54       0      1. 0.66667 0.33333       0        
GRID          55       0      1. 0.33333 0.66667       0        
GRID          56       0      1. 0.66667 0.66667       0        
GRID          57       0 0.33333 0.33333 0.66667       0        
GRID          58       0 0.66667 0.33333 0.66667       0        
GRID          59       0 0.33333 0.66667 0.66667       0        
GRID          60       0 0.66667 0.66667 0.66667       0        
GRID          61       0 0.33333 0.33333 0.33333       0        
GRID          62       0 0.66667 0.33333 0.33333       0        
GRID          63       0 0.33333 0.66667 0.33333       0        
GRID          64       0 0.66667 0.66667 0.33333       0        
CHEXA          1       1      20      21       7       1      19      27+EL    1
+EL    1      57      17                                                        
CHEXA          2       1       1       7       8       2      17      57+EL    2
+EL    2      58      15                                                        
CHEXA          3       1       2       8       3      11      15      58+EL    3
+EL    3      55      51                                                        
CHEXA          4       1      21       6       9       7      27      29+EL    4
+EL    4      59      57                                                        
CHEXA          5       1       7       9      10       8      57      59+EL    5
+EL    5      60      58                                                        
CHEXA          6       1       8      10      50       3      58      60+EL    6
+EL    6      56      55                                                        
CHEXA          7       1       6      22       5       9      29      23+EL    7
+EL    7      34      59                                                        
CHEXA          8       1       9       5      31      10      59      34+EL    8
+EL    8      36      60                                                        
CHEXA          9       1      10      31       4      50      60      36+EL    9
+EL    9      32      56                                                        
CHEXA         10       1      19      27      57      17      18      26+EL    A
+EL    A      61      16                                                        
CHEXA         11       1      17      57      58      15      16      61+EL    B
+EL    B      62      14                                                        
CHEXA         12       1      15      58      55      51      14      62+EL    C
+EL    C      53      52                                                        
CHEXA         13       1      27      29      59      57      26      28+EL    D
+EL    D      63      61                                                        
CHEXA         14       1      57      59      60      58      61      63+EL    E
+EL    E      64      62                                                        
CHEXA         15       1      58      60      56      55      62      64+EL    F
+EL    F      54      53                                                        
CHEXA         16       1      29      23      34      59      28      30+EL    G
+EL    G      33      63                                                        
CHEXA         17       1      59      34      36      60      63      33+EL    H
+EL    H      35      64                                                        
CHEXA         18       1      60      36      32      56      64      35+EL    I
+EL    I      49      54                                                        
CHEXA         19       1      18      26      61      16      39      40+EL    J
+EL    J      46      12                                                        
CHEXA         20       1      16      61      62      14      12      46+EL    K
+EL    K      44      13                                                        
CHEXA         21       1      14      62      53      52      13      44+EL    L
+EL    L      38      47                                                        
CHEXA         22       1      26      28      63      61      40      25+EL    M
+EL    M      45      46                                                        
CHEXA         23       1      61      63      64      62      46      45+EL    N
+EL    N      43      44                                                        
CHEXA         24       1      62      64      54      53      44      43+EL    O
+EL    O      37      38                                                        
CHEXA         25       1      28      30      33      63      25      24+EL    P
+EL    P      41      45                                                        
CHEXA         26       1      63      33      35      64      45      41+EL    Q
+EL    Q      42      43                                                        
CHEXA         27       1      64      35      49      54      43      42+EL    R
+EL    R      48      37                                                        
ENDDATA
