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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Parametric\Solid\sphere1.MOD
$   Date       : Sat Nov 19 18:52:25 2005
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
$ FEMAP Property 1 : Sphere surface
PSHELL         1       1      0.       1               1              0.
$ FEMAP Material 1 : Mat1
MAT1           1                              0.      0.      0.        
GRID          57       0     -1.      0.      0.       0        
GRID          58       0-0.70711      0. 0.70711       0        
GRID          60       0      0. 0.70711 0.70711       0        
GRID          62       0-0.70711 0.70711      0.       0        
GRID          66       0      0.      0.     -1.       0        
GRID          67       0      0.-0.70711-0.70711       0        
GRID          72       0-0.70711-0.70711      0.       0        
GRID          73       0      0.     -1.      0.       0        
GRID          74       0      0.-0.70711 0.70711       0        
GRID          75       0      0.      0.      1.       0        
GRID          80       0      0.      1.      0.       0        
GRID          81       0      0. 0.70711-0.70711       0        
GRID          83       0-0.70711      0.-0.70711       0        
GRID          88       0 0.70711 0.70711      0.       0        
GRID          95       0 0.70711-0.70711      0.       0        
GRID          96       0      1.      0.      0.       0        
GRID          97       0 0.70711      0. 0.70711       0        
GRID         102       0 0.70711      0.-0.70711       0        
CTRIA3        33       1      58      75      60                        
CTRIA3        34       1      60      80      62                        
CTRIA3        35       1      58      60      62                        
CTRIA3        36       1      57      58      62                        
CTRIA3        45       1      83      66      67                        
CTRIA3        46       1      67      73      72                        
CTRIA3        47       1      83      67      72                        
CTRIA3        48       1      57      83      72                        
CTRIA3        57       1      74      75      58                        
CTRIA3        58       1      72      73      74                        
CTRIA3        59       1      72      74      58                        
CTRIA3        60       1      57      72      58                        
CTRIA3        69       1      81      66      83                        
CTRIA3        70       1      62      80      81                        
CTRIA3        71       1      62      81      83                        
CTRIA3        72       1      57      62      83                        
CTRIA3        81       1      88      96     102                        
CTRIA3        82       1      81      80      88                        
CTRIA3        83       1      81      88     102                        
CTRIA3        84       1      66      81     102                        
CTRIA3        93       1      95      96      97                        
CTRIA3        94       1      74      73      95                        
CTRIA3        95       1      74      95      97                        
CTRIA3        96       1      75      74      97                        
CTRIA3       105       1     102      96      95                        
CTRIA3       106       1      67      66     102                        
CTRIA3       107       1      67     102      95                        
CTRIA3       108       1      73      67      95                        
CTRIA3       117       1      97      96      88                        
CTRIA3       118       1      60      75      97                        
CTRIA3       119       1      60      97      88                        
CTRIA3       120       1      80      60      88                        
ENDDATA
