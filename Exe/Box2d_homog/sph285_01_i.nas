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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box2d_homog\Solid\model.mod
$   Date       : Thu Jun 14 14:53:47 2007
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
$ FEMAP Load Set 111 : Inclusion BSOURCE 0.4806 0.4806 0.285 0.285 0.
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : mat
MAT1           1                              0.      0.      0.        
GRID           9       0 0.24372 0.71748      0.       0        
GRID          10       0  0.1456  0.4806      0.       0        
GRID          14       0 0.12186 0.83934      0.       0        
GRID          15       0      0.  0.9612      0.       0        
GRID          16       0      0.  0.7209      0.       0        
GRID          18       00.089153 0.68028      0.       0        
GRID          22       0 0.27907 0.68213      0.       0        
GRID          25       0  0.4806  0.7656      0.       0        
GRID          28       0  0.2403  0.9612      0.       0        
GRID          29       0  0.4806  0.9612      0.       0        
GRID          30       0 0.28092 0.87205      0.       0        
GRID          35       0 0.71748 0.71748      0.       0        
GRID          36       0  0.8156  0.4806      0.       0        
GRID          37       0  0.7656  0.4806      0.       0        
GRID          40       0  0.9612  0.7209      0.       0        
GRID          42       0 0.87205 0.68028      0.       0        
GRID          48       0  0.4806  0.8156      0.       0        
GRID          50       0 0.83934 0.83934      0.       0        
GRID          51       0  0.9612  0.9612      0.       0        
GRID          52       0  0.7209  0.9612      0.       0        
GRID          54       0 0.68028 0.87205      0.       0        
GRID          56       0 0.68213 0.68213      0.       0        
GRID          59       0 0.24372 0.24372      0.       0        
GRID          61       0  0.1956  0.4806      0.       0        
GRID          62       0 0.12186 0.12186      0.       0        
GRID          63       0      0.      0.      0.       0        
GRID          64       0      0.  0.2403      0.       0        
GRID          65       0      0.  0.4806      0.       0        
GRID          66       00.089153 0.28092      0.       0        
GRID          73       0  0.4806  0.1956      0.       0        
GRID          76       0  0.2403      0.      0.       0        
GRID          78       0 0.280920.089153      0.       0        
GRID          79       0  0.4806  0.4806      0.       0        
GRID          80       0 0.27907 0.27907      0.       0        
GRID          82       0 0.68213 0.27907      0.       0        
GRID          86       0 0.83934 0.12186      0.       0        
GRID          88       0  0.9612  0.2403      0.       0        
GRID          89       0  0.9612  0.4806      0.       0        
GRID          90       0 0.87205 0.28092      0.       0        
GRID          95       0 0.71748 0.24372      0.       0        
GRID          96       0  0.4806  0.1456      0.       0        
GRID          99       0  0.9612      0.      0.       0        
GRID         100       0  0.7209      0.      0.       0        
GRID         101       0  0.4806      0.      0.       0        
GRID         102       0 0.680280.089153      0.       0        
CQUAD4         4       1      22       9      10      61                
CQUAD4         5       1      14      15      16      18                
CTRIA3         6       1      10       9      18                        
CQUAD4         7       1      10      18      16      65                
CTRIA3         8       1       9      14      18                        
CTRIA3         9       2      79      22      61                        
CQUAD4        10       1       9      22      25      48                
CQUAD4        11       1      15      14      30      28                
CTRIA3        12       1       9      48      30                        
CQUAD4        13       1      30      48      29      28                
CTRIA3        14       1      14       9      30                        
CTRIA3        15       2      22      79      25                        
CQUAD4        16       1      35      56      37      36                
CQUAD4        17       1      51      50      42      40                
CTRIA3        18       1      35      36      42                        
CQUAD4        19       1      42      36      89      40                
CTRIA3        20       1      50      35      42                        
CTRIA3        21       2      56      79      37                        
CQUAD4        22       1      56      35      48      25                
CQUAD4        23       1      50      51      52      54                
CTRIA3        24       1      48      35      54                        
CQUAD4        25       1      48      54      52      29                
CTRIA3        26       1      35      50      54                        
CTRIA3        27       2      79      56      25                        
CQUAD4        28       1      59      80      61      10                
CQUAD4        29       1      63      62      66      64                
CTRIA3        30       1      59      10      66                        
CQUAD4        31       1      66      10      65      64                
CTRIA3        32       1      62      59      66                        
CTRIA3        33       2      80      79      61                        
CQUAD4        34       1      80      59      96      73                
CQUAD4        35       1      62      63      76      78                
CTRIA3        36       1      96      59      78                        
CQUAD4        37       1      96      78      76     101                
CTRIA3        38       1      59      62      78                        
CTRIA3        39       2      79      80      73                        
CQUAD4        40       1      82      95      36      37                
CQUAD4        41       1      86      99      88      90                
CTRIA3        42       1      36      95      90                        
CQUAD4        43       1      36      90      88      89                
CTRIA3        44       1      95      86      90                        
CTRIA3        45       2      79      82      37                        
CQUAD4        46       1      95      82      73      96                
CQUAD4        47       1      99      86     102     100                
CTRIA3        48       1      95      96     102                        
CQUAD4        49       1     102      96     101     100                
CTRIA3        50       1      86      95     102                        
CTRIA3        51       2      82      79      73                        
ENDDATA
