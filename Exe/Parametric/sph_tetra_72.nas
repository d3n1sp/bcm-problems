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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box3d_fulleren\Solid_sph\sphere1.MOD
$   Date       : Sat Nov 19 14:30:44 2005
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
GRID         145       0     -1.      0.      0.       0        
GRID         149       0      0.     0.5 0.86603       0        
GRID         150       0      0. 0.86603     0.5       0        
GRID         151       0      0.      1.      0.       0        
GRID         152       0    -0.5 0.86603      0.       0        
GRID         153       0-0.86603     0.5      0.       0        
GRID         158       0    -0.5 0.61237 0.61237       0        
GRID         162       0      0.      0.     -1.       0        
GRID         163       0      0.    -0.5-0.86603       0        
GRID         164       0      0.-0.86603    -0.5       0        
GRID         165       0      0.     -1.      0.       0        
GRID         167       0-0.86603    -0.5      0.       0        
GRID         172       0    -0.5-0.61237-0.61237       0        
GRID         175       0    -0.5-0.86603      0.       0        
GRID         177       0      0.-0.86603     0.5       0        
GRID         178       0      0.    -0.5 0.86603       0        
GRID         179       0      0.      0.      1.       0        
GRID         180       0    -0.5      0. 0.86603       0        
GRID         181       0-0.86603      0.     0.5       0        
GRID         186       0    -0.5-0.61237 0.61237       0        
GRID         191       0      0. 0.86603    -0.5       0        
GRID         192       0      0.     0.5-0.86603       0        
GRID         194       0    -0.5      0.-0.86603       0        
GRID         195       0-0.86603      0.    -0.5       0        
GRID         200       0    -0.5 0.61237-0.61237       0        
GRID         205       0     0.5 0.86603      0.       0        
GRID         209       0     0.5      0.-0.86603       0        
GRID         214       0     0.5 0.61237-0.61237       0        
GRID         220       0 0.86603    -0.5      0.       0        
GRID         222       0 0.86603      0.     0.5       0        
GRID         223       0     0.5      0. 0.86603       0        
GRID         228       0     0.5-0.61237 0.61237       0        
GRID         234       0 0.86603      0.    -0.5       0        
GRID         237       0     0.5-0.86603      0.       0        
GRID         242       0     0.5-0.61237-0.61237       0        
GRID         249       0      1.      0.      0.       0        
GRID         250       0 0.86603     0.5      0.       0        
GRID         256       0     0.5 0.61237 0.61237       0        
CTRIA3        91       1     150     151     152                        
CTRIA3        92       1     150     152     158                        
CTRIA3        93       1     180     179     149                        
CTRIA3        94       1     149     150     158                        
CTRIA3        95       1     180     149     158                        
CTRIA3        96       1     181     180     158                        
CTRIA3        97       1     158     152     153                        
CTRIA3        98       1     181     158     153                        
CTRIA3        99       1     145     181     153                        
CTRIA3       118       1     164     165     175                        
CTRIA3       119       1     164     175     172                        
CTRIA3       120       1     163     164     172                        
CTRIA3       121       1     194     162     163                        
CTRIA3       122       1     194     163     172                        
CTRIA3       123       1     172     175     167                        
CTRIA3       124       1     195     194     172                        
CTRIA3       125       1     195     172     167                        
CTRIA3       126       1     145     195     167                        
CTRIA3       145       1     178     179     180                        
CTRIA3       146       1     178     180     186                        
CTRIA3       147       1     177     178     186                        
CTRIA3       148       1     175     165     177                        
CTRIA3       149       1     175     177     186                        
CTRIA3       150       1     186     180     181                        
CTRIA3       151       1     167     175     186                        
CTRIA3       152       1     167     186     181                        
CTRIA3       153       1     145     167     181                        
CTRIA3       172       1     192     162     194                        
CTRIA3       173       1     192     194     200                        
CTRIA3       174       1     191     192     200                        
CTRIA3       175       1     152     151     191                        
CTRIA3       176       1     152     191     200                        
CTRIA3       177       1     200     194     195                        
CTRIA3       178       1     153     152     200                        
CTRIA3       179       1     153     200     195                        
CTRIA3       180       1     145     153     195                        
CTRIA3       199       1     250     249     234                        
CTRIA3       200       1     234     209     214                        
CTRIA3       201       1     250     234     214                        
CTRIA3       202       1     205     250     214                        
CTRIA3       203       1     191     151     205                        
CTRIA3       204       1     191     205     214                        
CTRIA3       205       1     192     191     214                        
CTRIA3       206       1     192     214     209                        
CTRIA3       207       1     162     192     209                        
CTRIA3       226       1     220     249     222                        
CTRIA3       227       1     222     223     228                        
CTRIA3       228       1     220     222     228                        
CTRIA3       229       1     237     220     228                        
CTRIA3       230       1     177     165     237                        
CTRIA3       231       1     177     237     228                        
CTRIA3       232       1     178     177     228                        
CTRIA3       233       1     178     228     223                        
CTRIA3       234       1     179     178     223                        
CTRIA3       253       1     234     249     220                        
CTRIA3       254       1     220     237     242                        
CTRIA3       255       1     234     220     242                        
CTRIA3       256       1     209     234     242                        
CTRIA3       257       1     163     162     209                        
CTRIA3       258       1     163     209     242                        
CTRIA3       259       1     164     163     242                        
CTRIA3       260       1     164     242     237                        
CTRIA3       261       1     165     164     237                        
CTRIA3       280       1     222     249     250                        
CTRIA3       281       1     250     205     256                        
CTRIA3       282       1     222     250     256                        
CTRIA3       283       1     223     222     256                        
CTRIA3       284       1     149     179     223                        
CTRIA3       285       1     149     223     256                        
CTRIA3       286       1     150     149     256                        
CTRIA3       287       1     150     256     205                        
CTRIA3       288       1     151     150     205                        
ENDDATA
