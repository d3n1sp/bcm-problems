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
$   Date       : Tue Feb 09 23:24:16 2010
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.5 0.5 0.25 0.25 0.
$ FEMAP Property 1 : Untitled
PSHELL         1       1      0.       1      1.       1 0.83333      0.
$ FEMAP Property 2 : Untitled
PSHELL         2       1      0.       1      1.       1 0.83333      0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID           1       0 0.32322 0.67678      0.       0        
GRID           2       0      0.      1.      0.       0        
GRID           3       0     0.5    0.75      0.       0        
GRID           4       0 0.67678 0.67678      0.       0        
GRID           5       0      1.      1.      0.       0        
GRID           6       0     0.5      1.      0.       0        
GRID           7       0    0.75     0.5      0.       0        
GRID           8       0 0.32322 0.32322      0.       0        
GRID           9       0      0.      0.      0.       0        
GRID          10       0     0.5      0.      0.       0        
GRID          11       0    0.25     0.5      0.       0        
GRID          12       0      0.     0.5      0.       0        
GRID          13       0     0.5    0.25      0.       0        
GRID          14       0 0.67678 0.32322      0.       0        
GRID          15       0      1.      0.      0.       0        
GRID          16       0     0.5     0.5      0.       0        
GRID          17       0      1.     0.5      0.       0        
CTRIA3         1       2      11      16       1                        
CTRIA3         2       2      16       3       1                        
CQUAD4         3       1      12      11       1       2                
CQUAD4         4       1       1       3       6       2                
CTRIA3         5       2      16       7       4                        
CTRIA3         6       2       3      16       4                        
CQUAD4         7       1       7      17       5       4                
CQUAD4         8       1       3       4       5       6                
CTRIA3         9       2      16      11       8                        
CTRIA3        10       2      13      16       8                        
CQUAD4        11       1      11      12       9       8                
CQUAD4        12       1      13       8       9      10                
CTRIA3        13       2       7      16      14                        
CTRIA3        14       2      16      13      14                        
CQUAD4        15       1      17       7      14      15                
CQUAD4        16       1      14      13      10      15                
ENDDATA
