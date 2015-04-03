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
$   Date       : Tue Nov 22 17:55:45 2005
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
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1                              0.      0.      0.        
GRID           1       0     -1.      0.      0.       0        
GRID           2       0      0.      1. 2.59809       0        
GRID           3       0      0. 1.73206     1.5       0        
GRID           4       0      0.      2.      0.       0        
GRID           5       0    -0.5 1.73206      0.       0        
GRID           6       0-0.86603      1.      0.       0        
GRID           7       0    -0.5 1.22474 1.83711       0        
GRID           8       0      0.      0.     -3.       0        
GRID           9       0      0.     -1.-2.59809       0        
GRID          10       0      0.-1.73206    -1.5       0        
GRID          11       0      0.     -2.      0.       0        
GRID          12       0-0.86603     -1.      0.       0        
GRID          13       0    -0.5-1.22474-1.83711       0        
GRID          14       0    -0.5-1.73206      0.       0        
GRID          15       0      0.-1.73206     1.5       0        
GRID          16       0      0.     -1. 2.59809       0        
GRID          17       0      0.      0.      3.       0        
GRID          18       0    -0.5      0. 2.59809       0        
GRID          19       0-0.86603      0.     1.5       0        
GRID          20       0    -0.5-1.22474 1.83711       0        
GRID          21       0      0. 1.73206    -1.5       0        
GRID          22       0      0.      1.-2.59809       0        
GRID          23       0    -0.5      0.-2.59809       0        
GRID          24       0-0.86603      0.    -1.5       0        
GRID          25       0    -0.5 1.22474-1.83711       0        
GRID          26       0     0.5 1.73206      0.       0        
GRID          27       0     0.5      0.-2.59809       0        
GRID          28       0     0.5 1.22474-1.83711       0        
GRID          29       0 0.86603     -1.      0.       0        
GRID          30       0 0.86603      0.     1.5       0        
GRID          31       0     0.5      0. 2.59809       0        
GRID          32       0     0.5-1.22474 1.83711       0        
GRID          33       0 0.86603      0.    -1.5       0        
GRID          34       0     0.5-1.73206      0.       0        
GRID          35       0     0.5-1.22474-1.83711       0        
GRID          36       0      1.      0.      0.       0        
GRID          37       0 0.86603      1.      0.       0        
GRID          38       0     0.5 1.22474 1.83711       0        
CTRIA3         1       1       3       4       5                        
CTRIA3         2       1       3       5       7                        
CTRIA3         3       1      18      17       2                        
CTRIA3         4       1       2       3       7                        
CTRIA3         5       1      18       2       7                        
CTRIA3         6       1      19      18       7                        
CTRIA3         7       1       7       5       6                        
CTRIA3         8       1      19       7       6                        
CTRIA3         9       1       1      19       6                        
CTRIA3        10       1      10      11      14                        
CTRIA3        11       1      10      14      13                        
CTRIA3        12       1       9      10      13                        
CTRIA3        13       1      23       8       9                        
CTRIA3        14       1      23       9      13                        
CTRIA3        15       1      13      14      12                        
CTRIA3        16       1      24      23      13                        
CTRIA3        17       1      24      13      12                        
CTRIA3        18       1       1      24      12                        
CTRIA3        19       1      16      17      18                        
CTRIA3        20       1      16      18      20                        
CTRIA3        21       1      15      16      20                        
CTRIA3        22       1      14      11      15                        
CTRIA3        23       1      14      15      20                        
CTRIA3        24       1      20      18      19                        
CTRIA3        25       1      12      14      20                        
CTRIA3        26       1      12      20      19                        
CTRIA3        27       1       1      12      19                        
CTRIA3        28       1      22       8      23                        
CTRIA3        29       1      22      23      25                        
CTRIA3        30       1      21      22      25                        
CTRIA3        31       1       5       4      21                        
CTRIA3        32       1       5      21      25                        
CTRIA3        33       1      25      23      24                        
CTRIA3        34       1       6       5      25                        
CTRIA3        35       1       6      25      24                        
CTRIA3        36       1       1       6      24                        
CTRIA3        37       1      37      36      33                        
CTRIA3        38       1      33      27      28                        
CTRIA3        39       1      37      33      28                        
CTRIA3        40       1      26      37      28                        
CTRIA3        41       1      21       4      26                        
CTRIA3        42       1      21      26      28                        
CTRIA3        43       1      22      21      28                        
CTRIA3        44       1      22      28      27                        
CTRIA3        45       1       8      22      27                        
CTRIA3        46       1      29      36      30                        
CTRIA3        47       1      30      31      32                        
CTRIA3        48       1      29      30      32                        
CTRIA3        49       1      34      29      32                        
CTRIA3        50       1      15      11      34                        
CTRIA3        51       1      15      34      32                        
CTRIA3        52       1      16      15      32                        
CTRIA3        53       1      16      32      31                        
CTRIA3        54       1      17      16      31                        
CTRIA3        55       1      33      36      29                        
CTRIA3        56       1      29      34      35                        
CTRIA3        57       1      33      29      35                        
CTRIA3        58       1      27      33      35                        
CTRIA3        59       1       9       8      27                        
CTRIA3        60       1       9      27      35                        
CTRIA3        61       1      10       9      35                        
CTRIA3        62       1      10      35      34                        
CTRIA3        63       1      11      10      34                        
CTRIA3        64       1      30      36      37                        
CTRIA3        65       1      37      26      38                        
CTRIA3        66       1      30      37      38                        
CTRIA3        67       1      31      30      38                        
CTRIA3        68       1       2      17      31                        
CTRIA3        69       1       2      31      38                        
CTRIA3        70       1       3       2      38                        
CTRIA3        71       1       3      38      26                        
CTRIA3        72       1       4       3      26                        
ENDDATA
