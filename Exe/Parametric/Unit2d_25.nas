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
$   Date       : Sun Apr 04 16:43:24 2010
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
GRID           1       0      0.      1.      0.       0        
GRID           2       0      0.     0.8      0.       0        
GRID           3       0      0.     0.6      0.       0        
GRID           4       0      0.     0.4      0.       0        
GRID           5       0      0.     0.2      0.       0        
GRID           6       0      0.      0.      0.       0        
GRID           7       0     0.2      0.      0.       0        
GRID           8       0     0.4      0.      0.       0        
GRID           9       0     0.6      0.      0.       0        
GRID          10       0     0.8      0.      0.       0        
GRID          11       0      1.      0.      0.       0        
GRID          12       0      1.     0.2      0.       0        
GRID          13       0      1.     0.4      0.       0        
GRID          14       0      1.     0.6      0.       0        
GRID          15       0      1.     0.8      0.       0        
GRID          16       0      1.      1.      0.       0        
GRID          17       0     0.8      1.      0.       0        
GRID          18       0     0.6      1.      0.       0        
GRID          19       0     0.4      1.      0.       0        
GRID          20       0     0.2      1.      0.       0        
GRID          21       0     0.2     0.8      0.       0        
GRID          22       0     0.2     0.6      0.       0        
GRID          23       0     0.2     0.4      0.       0        
GRID          24       0     0.2     0.2      0.       0        
GRID          25       0     0.4     0.8      0.       0        
GRID          26       0     0.4     0.6      0.       0        
GRID          27       0     0.4     0.4      0.       0        
GRID          28       0     0.4     0.2      0.       0        
GRID          29       0     0.6     0.8      0.       0        
GRID          30       0     0.6     0.6      0.       0        
GRID          31       0     0.6     0.4      0.       0        
GRID          32       0     0.6     0.2      0.       0        
GRID          33       0     0.8     0.8      0.       0        
GRID          34       0     0.8     0.6      0.       0        
GRID          35       0     0.8     0.4      0.       0        
GRID          36       0     0.8     0.2      0.       0        
CQUAD4         1       1       1       2      21      20                
CQUAD4         2       1       2       3      22      21                
CQUAD4         3       1       3       4      23      22                
CQUAD4         4       1       4       5      24      23                
CQUAD4         5       1       5       6       7      24                
CQUAD4         6       1      20      21      25      19                
CQUAD4         7       1      21      22      26      25                
CQUAD4         8       1      22      23      27      26                
CQUAD4         9       1      23      24      28      27                
CQUAD4        10       1      24       7       8      28                
CQUAD4        11       1      19      25      29      18                
CQUAD4        12       1      25      26      30      29                
CQUAD4        13       1      26      27      31      30                
CQUAD4        14       1      27      28      32      31                
CQUAD4        15       1      28       8       9      32                
CQUAD4        16       1      18      29      33      17                
CQUAD4        17       1      29      30      34      33                
CQUAD4        18       1      30      31      35      34                
CQUAD4        19       1      31      32      36      35                
CQUAD4        20       1      32       9      10      36                
CQUAD4        21       1      17      33      15      16                
CQUAD4        22       1      33      34      14      15                
CQUAD4        23       1      34      35      13      14                
CQUAD4        24       1      35      36      12      13                
CQUAD4        25       1      36      10      11      12                
ENDDATA
