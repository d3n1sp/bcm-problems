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
$   Date       : Sat May 20 06:29:58 2006
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.85425 0.85425 0.623 0.623 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID          24       0      0.  1.7085      0.       0        
GRID          25       0      0. 1.46443      0.       0        
GRID          26       0      0. 1.22036      0.       0        
GRID          27       0      0. 0.97629      0.       0        
GRID          28       0      0. 0.73221      0.       0        
GRID          29       0      0. 0.48814      0.       0        
GRID          30       0      0. 0.24407      0.       0        
GRID          31       0      0.      0.      0.       0        
GRID          32       0 0.24407      0.      0.       0        
GRID          33       0 0.48814      0.      0.       0        
GRID          34       0 0.73221      0.      0.       0        
GRID          35       0 0.97629      0.      0.       0        
GRID          36       0 1.22036      0.      0.       0        
GRID          37       0 1.46443      0.      0.       0        
GRID          38       0  1.7085      0.      0.       0        
GRID          39       0  1.7085 0.24407      0.       0        
GRID          40       0  1.7085 0.48814      0.       0        
GRID          41       0  1.7085 0.73221      0.       0        
GRID          42       0  1.7085 0.97629      0.       0        
GRID          43       0  1.7085 1.22036      0.       0        
GRID          44       0  1.7085 1.46443      0.       0        
GRID          45       0  1.7085  1.7085      0.       0        
GRID          46       0 1.46443  1.7085      0.       0        
GRID          47       0 1.22036  1.7085      0.       0        
GRID          48       0 0.97629  1.7085      0.       0        
GRID          49       0 0.73221  1.7085      0.       0        
GRID          50       0 0.48814  1.7085      0.       0        
GRID          51       0 0.24407  1.7085      0.       0        
GRID          52       0 0.23125 0.85425      0.       0        
GRID          53       0 0.27867 1.09266      0.       0        
GRID          54       0 0.41372 1.29478      0.       0        
GRID          55       0 0.61584 1.42983      0.       0        
GRID          57       0 1.09266 1.42983      0.       0        
GRID          60       0 1.47725 0.85425      0.       0        
GRID          61       0 1.42983 0.61584      0.       0        
GRID          62       0 1.29478 0.41372      0.       0        
GRID          63       0 1.09266 0.27867      0.       0        
GRID          66       0 0.41372 0.41372      0.       0        
GRID          68       0 1.38087 0.22219      0.       0        
GRID          69       0 1.50454 0.39674      0.       0        
GRID          70       0 1.48631 1.38087      0.       0        
GRID          71       0 1.31176 1.50454      0.       0        
GRID          72       0 0.22219 0.32763      0.       0        
GRID          73       0 0.39674 0.20396      0.       0        
GRID          74       0 0.32763 1.48631      0.       0        
GRID          75       0 0.20396 1.31176      0.       0        
GRID          77       0 1.42983 1.09266      0.       0        
GRID          78       0 1.29478 1.29478      0.       0        
GRID          80       0 0.85425 1.47725      0.       0        
GRID          85       0 0.27867 0.61584      0.       0        
GRID          87       0 0.61584 0.27867      0.       0        
GRID          88       0 0.85425 0.23125      0.       0        
GRID          92       0 0.47214 0.69598      0.       0        
GRID          93       0  0.7248 0.80063      0.       0        
GRID          94       0  0.9837 0.90787      0.       0        
GRID          95       0 1.23636 1.01252      0.       0        
GRID          96       0 1.03396 1.16356      0.       0        
GRID          97       0 0.82321 1.28861      0.       0        
GRID          98       0 0.77389 1.04827      0.       0        
GRID          99       0 0.50846 0.94589      0.       0        
GRID         100       0 0.56906 1.18334      0.       0        
GRID         101       0 1.13944 0.52516      0.       0        
GRID         102       0 1.20004 0.76261      0.       0        
GRID         103       0 0.93461 0.66023      0.       0        
GRID         104       0 0.88529 0.41989      0.       0        
GRID         105       0 0.67454 0.54494      0.       0        
CTRIA3       124       1      87      33      34                        
CTRIA3       125       1      88      87      34                        
CTRIA3       126       1      88      34      35                        
CTRIA3       127       1      63      88      35                        
CTRIA3       128       1      63      35      36                        
CTRIA3       129       1      36      37      68                        
CTRIA3       130       1      63      36      68                        
CTRIA3       131       1      62      63      68                        
CTRIA3       132       1      39      40      69                        
CTRIA3       133       1      37      38      39                        
CTRIA3       134       1      68      37      39                        
CTRIA3       135       1      68      39      69                        
CTRIA3       136       1      62      68      69                        
CTRIA3       137       1      61      62      69                        
CTRIA3       138       1      61      69      40                        
CTRIA3       139       1      61      40      41                        
CTRIA3       140       1      60      61      41                        
CTRIA3       141       1      60      41      42                        
CTRIA3       142       1      77      60      42                        
CTRIA3       143       1      77      42      43                        
CTRIA3       144       1      43      44      70                        
CTRIA3       145       1      77      43      70                        
CTRIA3       146       1      78      77      70                        
CTRIA3       147       1      46      47      71                        
CTRIA3       148       1      44      45      46                        
CTRIA3       149       1      70      44      46                        
CTRIA3       150       1      70      46      71                        
CTRIA3       151       1      78      70      71                        
CTRIA3       152       1      57      78      71                        
CTRIA3       153       1      57      71      47                        
CTRIA3       154       1      66      72      73                        
CTRIA3       155       1      87      66      73                        
CTRIA3       156       1      33      87      73                        
CTRIA3       157       1      32      33      73                        
CTRIA3       158       1      30      31      32                        
CTRIA3       159       1      32      73      72                        
CTRIA3       160       1      30      32      72                        
CTRIA3       161       1      29      30      72                        
CTRIA3       162       1      72      66      85                        
CTRIA3       163       1      29      72      85                        
CTRIA3       164       1      28      29      85                        
CTRIA3       165       1      28      85      52                        
CTRIA3       166       1      27      28      52                        
CTRIA3       167       1      27      52      53                        
CTRIA3       168       1      26      27      53                        
CTRIA3       169       1      57      47      48                        
CTRIA3       170       1      80      57      48                        
CTRIA3       171       1      80      48      49                        
CTRIA3       172       1      55      80      49                        
CTRIA3       173       1      55      49      50                        
CTRIA3       174       1      50      51      74                        
CTRIA3       175       1      55      50      74                        
CTRIA3       176       1      54      55      74                        
CTRIA3       177       1      54      74      75                        
CTRIA3       178       1      53      54      75                        
CTRIA3       179       1      26      53      75                        
CTRIA3       180       1      25      26      75                        
CTRIA3       181       1      25      75      74                        
CTRIA3       182       1      25      74      51                        
CTRIA3       183       1      24      25      51                        
CTRIA3       184       2      96      97      98                        
CTRIA3       185       2      94      95      96                        
CTRIA3       186       2      94      96      98                        
CTRIA3       187       2      93      94      98                        
CTRIA3       188       2      52      85      92                        
CTRIA3       189       2      52      92      99                        
CTRIA3       190       2      53      52      99                        
CTRIA3       191       2      53      99     100                        
CTRIA3       192       2      54      53     100                        
CTRIA3       193       2      55      54     100                        
CTRIA3       194       2      99      92      93                        
CTRIA3       195       2      99      93      98                        
CTRIA3       196       2     100      99      98                        
CTRIA3       197       2     100      98      97                        
CTRIA3       198       2      55     100      97                        
CTRIA3       199       2      80      55      97                        
CTRIA3       200       2      57      80      97                        
CTRIA3       201       2      57      97      96                        
CTRIA3       202       2      78      57      96                        
CTRIA3       203       2      78      96      95                        
CTRIA3       204       2      77      78      95                        
CTRIA3       205       2      88      63     104                        
CTRIA3       206       2      87      88     104                        
CTRIA3       207       2      87     104     105                        
CTRIA3       208       2      66      87     105                        
CTRIA3       209       2      92      85      66                        
CTRIA3       210       2      92      66     105                        
CTRIA3       211       2     104      63     101                        
CTRIA3       212       2     104     101     103                        
CTRIA3       213       2     105     104     103                        
CTRIA3       214       2      93      92     105                        
CTRIA3       215       2      93     105     103                        
CTRIA3       216       2     103     101     102                        
CTRIA3       217       2      94      93     103                        
CTRIA3       218       2      94     103     102                        
CTRIA3       219       2      95      94     102                        
CTRIA3       220       2      60      77      95                        
CTRIA3       221       2      60      95     102                        
CTRIA3       222       2     101      63      62                        
CTRIA3       223       2     101      62      61                        
CTRIA3       224       2     102     101      61                        
CTRIA3       225       2      60     102      61                        
ENDDATA
