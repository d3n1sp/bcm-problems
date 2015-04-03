ID c:\Users\Dima\Lame3d2\IDE,FEMAP
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
$   From Model : c:\Users\Dima\Lame3d2\IDENT_17may2006\Ident\Exe\box2d_fulleren\Solid\model.MOD
$   Date       : Sat May 20 06:39:28 2006
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.73705 0.73705 0.623 0.623 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID          24       0      0.  1.4741      0.       0        
GRID          25       0      0. 1.22842      0.       0        
GRID          26       0      0. 0.98273      0.       0        
GRID          27       0      0. 0.73705      0.       0        
GRID          28       0      0. 0.49137      0.       0        
GRID          29       0      0. 0.24568      0.       0        
GRID          30       0      0.      0.      0.       0        
GRID          31       0 0.24568      0.      0.       0        
GRID          32       0 0.49137      0.      0.       0        
GRID          33       0 0.73705      0.      0.       0        
GRID          34       0 0.98273      0.      0.       0        
GRID          35       0 1.22842      0.      0.       0        
GRID          36       0  1.4741      0.      0.       0        
GRID          37       0  1.4741 0.24568      0.       0        
GRID          38       0  1.4741 0.49137      0.       0        
GRID          39       0  1.4741 0.73705      0.       0        
GRID          40       0  1.4741 0.98273      0.       0        
GRID          41       0  1.4741 1.22842      0.       0        
GRID          42       0  1.4741  1.4741      0.       0        
GRID          43       0 1.22842  1.4741      0.       0        
GRID          44       0 0.98273  1.4741      0.       0        
GRID          45       0 0.73705  1.4741      0.       0        
GRID          46       0 0.49137  1.4741      0.       0        
GRID          47       0 0.24568  1.4741      0.       0        
GRID          49       0 0.16147 0.97546      0.       0        
GRID          50       0 0.29652 1.17758      0.       0        
GRID          52       0 0.73705 1.36005      0.       0        
GRID          53       0 0.97546 1.31263      0.       0        
GRID          55       0 1.31263 0.97546      0.       0        
GRID          56       0 1.36005 0.73705      0.       0        
GRID          58       0 1.17758 0.29652      0.       0        
GRID          61       0 0.49864 0.16147      0.       0        
GRID          64       0 1.33855 1.33855      0.       0        
GRID          65       0 1.33855 0.13555      0.       0        
GRID          66       0 0.13555 0.13555      0.       0        
GRID          67       0 0.13555 1.33855      0.       0        
GRID          70       0 1.17758 1.17758      0.       0        
GRID          73       0 0.49864 1.31263      0.       0        
GRID          76       0 0.11405 0.73705      0.       0        
GRID          77       0 0.16147 0.49864      0.       0        
GRID          78       0 0.29652 0.29652      0.       0        
GRID          80       0 0.73705 0.11405      0.       0        
GRID          81       0 0.97546 0.16147      0.       0        
GRID          83       0 1.31263 0.49864      0.       0        
GRID          84       0 0.35494 0.57878      0.       0        
GRID          85       0  0.6076 0.68343      0.       0        
GRID          86       0  0.8665 0.79067      0.       0        
GRID          87       0 1.11916 0.89532      0.       0        
GRID          88       0 0.91676 1.04636      0.       0        
GRID          89       0 0.70601 1.17141      0.       0        
GRID          90       0 0.65669 0.93107      0.       0        
GRID          91       0 0.39126 0.82869      0.       0        
GRID          92       0 0.45186 1.06614      0.       0        
GRID          93       0 1.02224 0.40796      0.       0        
GRID          94       0 1.08284 0.64541      0.       0        
GRID          95       0 0.81741 0.54303      0.       0        
GRID          96       0 0.76809 0.30269      0.       0        
GRID          97       0 0.55734 0.42774      0.       0        
CTRIA3       124       1      44      45      52                        
CTRIA3       125       1      44      52      53                        
CTRIA3       126       1      43      44      53                        
CTRIA3       127       1      43      53      70                        
CTRIA3       128       1      43      70      64                        
CTRIA3       129       1      42      43      64                        
CTRIA3       130       1      41      42      64                        
CTRIA3       131       1      41      64      70                        
CTRIA3       132       1      41      70      55                        
CTRIA3       133       1      40      41      55                        
CTRIA3       134       1      40      55      56                        
CTRIA3       135       1      39      40      56                        
CTRIA3       136       1      38      39      56                        
CTRIA3       137       1      38      56      83                        
CTRIA3       138       1      37      38      83                        
CTRIA3       139       1      37      83      58                        
CTRIA3       140       1      37      58      65                        
CTRIA3       141       1      36      37      65                        
CTRIA3       142       1      35      36      65                        
CTRIA3       143       1      35      65      58                        
CTRIA3       144       1      35      58      81                        
CTRIA3       145       1      34      35      81                        
CTRIA3       146       1      34      81      80                        
CTRIA3       147       1      33      34      80                        
CTRIA3       148       1      32      33      80                        
CTRIA3       149       1      32      80      61                        
CTRIA3       150       1      31      32      61                        
CTRIA3       151       1      31      61      78                        
CTRIA3       152       1      31      78      66                        
CTRIA3       153       1      30      31      66                        
CTRIA3       154       1      29      30      66                        
CTRIA3       155       1      29      66      78                        
CTRIA3       156       1      29      78      77                        
CTRIA3       157       1      28      29      77                        
CTRIA3       158       1      28      77      76                        
CTRIA3       159       1      27      28      76                        
CTRIA3       160       1      26      27      76                        
CTRIA3       161       1      26      76      49                        
CTRIA3       162       1      25      26      49                        
CTRIA3       163       1      25      49      50                        
CTRIA3       164       1      25      50      67                        
CTRIA3       165       1      24      25      67                        
CTRIA3       166       1      52      45      46                        
CTRIA3       167       1      73      52      46                        
CTRIA3       168       1      73      46      47                        
CTRIA3       169       1      50      73      47                        
CTRIA3       170       1      67      50      47                        
CTRIA3       171       1      24      67      47                        
CTRIA3       172       2      88      89      90                        
CTRIA3       173       2      86      87      88                        
CTRIA3       174       2      86      88      90                        
CTRIA3       175       2      85      86      90                        
CTRIA3       176       2      76      77      84                        
CTRIA3       177       2      76      84      91                        
CTRIA3       178       2      49      76      91                        
CTRIA3       179       2      49      91      92                        
CTRIA3       180       2      50      49      92                        
CTRIA3       181       2      73      50      92                        
CTRIA3       182       2      91      84      85                        
CTRIA3       183       2      91      85      90                        
CTRIA3       184       2      92      91      90                        
CTRIA3       185       2      92      90      89                        
CTRIA3       186       2      73      92      89                        
CTRIA3       187       2      52      73      89                        
CTRIA3       188       2      53      52      89                        
CTRIA3       189       2      53      89      88                        
CTRIA3       190       2      70      53      88                        
CTRIA3       191       2      70      88      87                        
CTRIA3       192       2      55      70      87                        
CTRIA3       193       2      80      81      96                        
CTRIA3       194       2      61      80      96                        
CTRIA3       195       2      61      96      97                        
CTRIA3       196       2      78      61      97                        
CTRIA3       197       2      84      77      78                        
CTRIA3       198       2      84      78      97                        
CTRIA3       199       2      96      81      93                        
CTRIA3       200       2      96      93      95                        
CTRIA3       201       2      97      96      95                        
CTRIA3       202       2      85      84      97                        
CTRIA3       203       2      85      97      95                        
CTRIA3       204       2      95      93      94                        
CTRIA3       205       2      86      85      95                        
CTRIA3       206       2      86      95      94                        
CTRIA3       207       2      87      86      94                        
CTRIA3       208       2      56      55      87                        
CTRIA3       209       2      56      87      94                        
CTRIA3       210       2      93      81      58                        
CTRIA3       211       2      93      58      83                        
CTRIA3       212       2      94      93      83                        
CTRIA3       213       2      56      94      83                        
ENDDATA
