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
$   Date       : Sat Apr 01 10:13:29 2006
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
GRID           1       0     0.5     0.5    0.25       0        
GRID           2       0 0.35566 0.35566 0.35566       0        
GRID           3       0      0.      0.      0.       0        
GRID           4       0      0.     0.5     0.5       0        
GRID           5       0    0.25     0.5     0.5       0        
GRID           6       0 0.35566 0.35566 0.64434       0        
GRID           7       0     0.5    0.25     0.5       0        
GRID           8       0 0.32322 0.32322     0.5       0        
GRID           9       0      0.      0.      1.       0        
GRID          10       0     0.5     0.5      1.       0        
GRID          11       0      0.     0.5      1.       0        
GRID          12       0      0.      0.     0.5       0        
GRID          13       0 0.67678     0.5 0.32322       0        
GRID          14       0     0.5 0.32322 0.32322       0        
GRID          15       0 0.64434 0.35566 0.35566       0        
GRID          16       0      1.      0.      0.       0        
GRID          17       0     0.5      0.     0.5       0        
GRID          18       0     0.5      0.      0.       0        
GRID          19       0    0.75     0.5     0.5       0        
GRID          20       0     0.5 0.32322 0.67678       0        
GRID          21       0 0.64434 0.35566 0.64434       0        
GRID          22       0     0.5     0.5     0.5       0        
GRID          23       0 0.67678 0.32322     0.5       0        
GRID          24       0      1.      0.      1.       0        
GRID          25       0      1.      0.     0.5       0        
GRID          26       0     0.5      0.      1.       0        
GRID          27       0 0.32322     0.5 0.32322       0        
GRID          28       0 0.35566 0.64434 0.35566       0        
GRID          29       0      0.      1.      0.       0        
GRID          30       0      0.     0.5      0.       0        
GRID          31       0 0.32322     0.5 0.67678       0        
GRID          32       0 0.35566 0.64434 0.64434       0        
GRID          33       0     0.5    0.75     0.5       0        
GRID          34       0 0.32322 0.67678     0.5       0        
GRID          35       0      0.      1.      1.       0        
GRID          36       0      0.      1.     0.5       0        
GRID          37       0     0.5 0.67678 0.32322       0        
GRID          38       0 0.64434 0.64434 0.35566       0        
GRID          39       0      1.      1.      0.       0        
GRID          40       0      1.     0.5     0.5       0        
GRID          41       0     0.5     0.5      0.       0        
GRID          42       0      1.     0.5      0.       0        
GRID          43       0     0.5      1.     0.5       0        
GRID          44       0     0.5      1.      0.       0        
GRID          45       0 0.67678     0.5 0.67678       0        
GRID          46       0     0.5     0.5    0.75       0        
GRID          47       0     0.5 0.67678 0.67678       0        
GRID          48       0 0.64434 0.64434 0.64434       0        
GRID          49       0 0.67678 0.67678     0.5       0        
GRID          50       0      1.      1.      1.       0        
GRID          51       0      1.     0.5      1.       0        
GRID          52       0      1.      1.     0.5       0        
GRID          53       0     0.5      1.      1.       0        
CTETRA         1       2       2       5      22      27                        
CTETRA         2       2       2      22       5       8                        
CTETRA         3       2       2       1      22      14                        
CTETRA         4       2       2      22       1      27                        
CTETRA         5       2       2       7      22       8                        
CTETRA         6       2       2      22       7      14                        
CPENTA         7       1       2       5      27       3       4      30        
CPENTA         8       1       3       4      12       2       5       8        
CPENTA         9       1       2       1      14       3      41      18        
CPENTA        10       1       3      41      30       2       1      27        
CPENTA        11       1       2       7       8       3      17      12        
CPENTA        12       1       3      17      18       2       7      14        
CTETRA        13       2      22       5       6      31                        
CTETRA        14       2       5      22       6       8                        
CTETRA        15       2      22      46       6      20                        
CTETRA        16       2      46      22       6      31                        
CTETRA        17       2      22       7       6       8                        
CTETRA        18       2       7      22       6      20                        
CPENTA        19       1      31       5       6      11       4       9        
CPENTA        20       1      12       4       9       8       5       6        
CPENTA        21       1      20      46       6      26      10       9        
CPENTA        22       1      11      10       9      31      46       6        
CPENTA        23       1       8       7       6      12      17       9        
CPENTA        24       1      26      17       9      20       7       6        
CTETRA        25       2      22      19      15      13                        
CTETRA        26       2      19      22      15      23                        
CTETRA        27       2      22       1      15      14                        
CTETRA        28       2       1      22      15      13                        
CTETRA        29       2      22       7      15      23                        
CTETRA        30       2       7      22      15      14                        
CPENTA        31       1      13      19      15      42      40      16        
CPENTA        32       1      25      40      16      23      19      15        
CPENTA        33       1      14       1      15      18      41      16        
CPENTA        34       1      42      41      16      13       1      15        
CPENTA        35       1      23       7      15      25      17      16        
CPENTA        36       1      18      17      16      14       7      15        
CTETRA        37       2      21      19      22      45                        
CTETRA        38       2      21      22      19      23                        
CTETRA        39       2      21      46      22      20                        
CTETRA        40       2      21      22      46      45                        
CTETRA        41       2      21       7      22      23                        
CTETRA        42       2      21      22       7      20                        
CPENTA        43       1      21      19      45      24      40      51        
CPENTA        44       1      24      40      25      21      19      23        
CPENTA        45       1      21      46      20      24      10      26        
CPENTA        46       1      24      10      51      21      46      45        
CPENTA        47       1      21       7      23      24      17      25        
CPENTA        48       1      24      17      26      21       7      20        
CTETRA        49       2      22       5      28      27                        
CTETRA        50       2       5      22      28      34                        
CTETRA        51       2      22       1      28      37                        
CTETRA        52       2       1      22      28      27                        
CTETRA        53       2      22      33      28      34                        
CTETRA        54       2      33      22      28      37                        
CPENTA        55       1      27       5      28      30       4      29        
CPENTA        56       1      36       4      29      34       5      28        
CPENTA        57       1      37       1      28      44      41      29        
CPENTA        58       1      30      41      29      27       1      28        
CPENTA        59       1      34      33      28      36      43      29        
CPENTA        60       1      44      43      29      37      33      28        
CTETRA        61       2      32       5      22      31                        
CTETRA        62       2      32      22       5      34                        
CTETRA        63       2      32      46      22      47                        
CTETRA        64       2      32      22      46      31                        
CTETRA        65       2      32      33      22      34                        
CTETRA        66       2      32      22      33      47                        
CPENTA        67       1      32       5      31      35       4      11        
CPENTA        68       1      35       4      36      32       5      34        
CPENTA        69       1      32      46      47      35      10      53        
CPENTA        70       1      35      10      11      32      46      31        
CPENTA        71       1      32      33      34      35      43      36        
CPENTA        72       1      35      43      53      32      33      47        
CTETRA        73       2      38      19      22      13                        
CTETRA        74       2      38      22      19      49                        
CTETRA        75       2      38       1      22      37                        
CTETRA        76       2      38      22       1      13                        
CTETRA        77       2      38      33      22      49                        
CTETRA        78       2      38      22      33      37                        
CPENTA        79       1      38      19      13      39      40      42        
CPENTA        80       1      39      40      52      38      19      49        
CPENTA        81       1      38       1      37      39      41      44        
CPENTA        82       1      39      41      42      38       1      13        
CPENTA        83       1      38      33      49      39      43      52        
CPENTA        84       1      39      43      44      38      33      37        
CTETRA        85       2      22      19      48      45                        
CTETRA        86       2      19      22      48      49                        
CTETRA        87       2      22      46      48      47                        
CTETRA        88       2      46      22      48      45                        
CTETRA        89       2      22      33      48      49                        
CTETRA        90       2      33      22      48      47                        
CPENTA        91       1      45      19      48      51      40      50        
CPENTA        92       1      52      40      50      49      19      48        
CPENTA        93       1      47      46      48      53      10      50        
CPENTA        94       1      51      10      50      45      46      48        
CPENTA        95       1      49      33      48      52      43      50        
CPENTA        96       1      53      43      50      47      33      48        
ENDDATA
