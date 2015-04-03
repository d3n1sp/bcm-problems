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
$   Date       : Sat Nov 19 14:28:53 2005
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
GRID          61       0      0.     0.5 0.86603       0        
GRID          62       0      0. 0.86603     0.5       0        
GRID          63       0      0.      1.      0.       0        
GRID          64       0    -0.5 0.86603      0.       0        
GRID          65       0-0.86603     0.5      0.       0        
GRID          66       0-0.45825  0.6717 0.58208       0        
GRID          67       0-0.75582 0.47388 0.45185       0        
GRID          71       0      0.      0.     -1.       0        
GRID          72       0      0.    -0.5-0.86603       0        
GRID          73       0      0.-0.86603    -0.5       0        
GRID          74       0      0.     -1.      0.       0        
GRID          76       0-0.86603    -0.5      0.       0        
GRID          77       0-0.45825 -0.6717-0.58208       0        
GRID          78       0-0.75582-0.47388-0.45185       0        
GRID          81       0    -0.5-0.86603      0.       0        
GRID          83       0      0.-0.86603     0.5       0        
GRID          84       0      0.    -0.5 0.86603       0        
GRID          85       0      0.      0.      1.       0        
GRID          86       0    -0.5      0. 0.86603       0        
GRID          87       0-0.86603      0.     0.5       0        
GRID          88       0-0.45825-0.58208  0.6717       0        
GRID          89       0-0.75582-0.45185 0.47388       0        
GRID          94       0      0. 0.86603    -0.5       0        
GRID          95       0      0.     0.5-0.86603       0        
GRID          97       0    -0.5      0.-0.86603       0        
GRID          98       0-0.86603      0.    -0.5       0        
GRID          99       0-0.45825 0.58208 -0.6717       0        
GRID         100       0-0.75582 0.45185-0.47388       0        
GRID         105       0     0.5 0.86603      0.       0        
GRID         109       0     0.5      0.-0.86603       0        
GRID         110       0 0.45825 0.58208 -0.6717       0        
GRID         111       0 0.75582 0.45185-0.47388       0        
GRID         117       0 0.86603    -0.5      0.       0        
GRID         119       0 0.86603      0.     0.5       0        
GRID         120       0     0.5      0. 0.86603       0        
GRID         121       0 0.45825-0.58208  0.6717       0        
GRID         122       0 0.75582-0.45185 0.47388       0        
GRID         128       0 0.86603      0.    -0.5       0        
GRID         131       0     0.5-0.86603      0.       0        
GRID         132       0 0.45825 -0.6717-0.58208       0        
GRID         133       0 0.75582-0.47388-0.45185       0        
GRID         140       0      1.      0.      0.       0        
GRID         141       0 0.86603     0.5      0.       0        
GRID         143       0 0.45825 0.58208  0.6717       0        
GRID         144       0 0.75582 0.45185 0.47388       0        
CQUAD4        25       1      62      63      64      66                
CTRIA3        26       1      86      85      61                        
CQUAD4        27       1      86      61      62      66                
CQUAD4        28       1      66      64      65      67                
CQUAD4        29       1      87      86      66      67                
CQUAD4        30       1      57      87      67      65                
CQUAD4        31       1      73      74      81      77                
CTRIA3        32       1      97      71      72                        
CQUAD4        33       1      97      72      73      77                
CQUAD4        34       1      77      81      76      78                
CQUAD4        35       1      98      97      77      78                
CQUAD4        36       1      57      98      78      76                
CQUAD4        37       1      84      85      86      88                
CTRIA3        38       1      81      74      83                        
CQUAD4        39       1      81      83      84      88                
CQUAD4        40       1      88      86      87      89                
CQUAD4        41       1      76      81      88      89                
CQUAD4        42       1      57      76      89      87                
CQUAD4        43       1      95      71      97      99                
CTRIA3        44       1      64      63      94                        
CQUAD4        45       1      64      94      95      99                
CQUAD4        46       1      99      97      98     100                
CQUAD4        47       1      65      64      99     100                
CQUAD4        48       1      57      65     100      98                
CQUAD4        49       1     141     140     128     111                
CQUAD4        50       1     111     128     109     110                
CQUAD4        51       1     105     141     111     110                
CTRIA3        52       1      94      63     105                        
CQUAD4        53       1      95      94     105     110                
CQUAD4        54       1      71      95     110     109                
CQUAD4        55       1     117     140     119     122                
CQUAD4        56       1     122     119     120     121                
CQUAD4        57       1     131     117     122     121                
CTRIA3        58       1      83      74     131                        
CQUAD4        59       1      84      83     131     121                
CQUAD4        60       1      85      84     121     120                
CQUAD4        61       1     128     140     117     133                
CQUAD4        62       1     133     117     131     132                
CQUAD4        63       1     109     128     133     132                
CTRIA3        64       1      72      71     109                        
CQUAD4        65       1      73      72     109     132                
CQUAD4        66       1      74      73     132     131                
CQUAD4        67       1     119     140     141     144                
CQUAD4        68       1     144     141     105     143                
CQUAD4        69       1     120     119     144     143                
CQUAD4        70       1      61      85     120     143                
CQUAD4        71       1      62      61     143     105                
CTRIA3        72       1      63      62     105                        
ENDDATA
