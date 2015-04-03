ID C:\Users\Dima\Lame3d2\IDE,FEMAP
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
$   From Model : C:\Users\Dima\Lame3d2\IDENT_12septmber2006\Ident\Exe\box2d_fulleren\Solid\model.MOD
$   Date       : Thu Nov 09 12:41:54 2006
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
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID         224       0      0.      0.      0.       0        
GRID         225       00.093795      0.      0.       0        
GRID         226       0 0.18759      0.      0.       0        
GRID         227       0 0.18759 0.12506      0.       0        
GRID         231       0 0.18759  0.6253      0.       0        
GRID         232       00.093795  0.6253      0.       0        
GRID         233       0      0.  0.6253      0.       0        
GRID         234       0      0. 0.50024      0.       0        
GRID         235       0      0. 0.37518      0.       0        
GRID         236       0      0. 0.25012      0.       0        
GRID         237       0      0. 0.12506      0.       0        
GRID         238       00.093795 0.12506      0.       0        
GRID         239       00.093795 0.25012      0.       0        
GRID         240       00.093795 0.37518      0.       0        
GRID         241       00.093795 0.50024      0.       0        
GRID         242       0  0.6253  0.6253      0.       0        
GRID         243       0 0.53151  0.6253      0.       0        
GRID         244       0 0.43771  0.6253      0.       0        
GRID         248       0 0.43771 0.12506      0.       0        
GRID         249       0 0.43771      0.      0.       0        
GRID         250       0 0.53151      0.      0.       0        
GRID         251       0  0.6253      0.      0.       0        
GRID         252       0  0.6253 0.12506      0.       0        
GRID         253       0  0.6253 0.25012      0.       0        
GRID         254       0  0.6253 0.37518      0.       0        
GRID         255       0  0.6253 0.50024      0.       0        
GRID         256       0 0.53151 0.50024      0.       0        
GRID         257       0  0.5315 0.37518      0.       0        
GRID         258       0  0.5315 0.25012      0.       0        
GRID         259       0 0.53151 0.12506      0.       0        
GRID         261       0 0.18759 0.50024      0.       0        
GRID         262       0 0.18759 0.37518      0.       0        
GRID         263       0 0.18759 0.25012      0.       0        
GRID         266       0 0.31265      0.      0.       0        
GRID         269       0 0.43771 0.25012      0.       0        
GRID         270       0 0.43771 0.37518      0.       0        
GRID         271       0 0.43771 0.50024      0.       0        
GRID         273       0 0.31265  0.6253      0.       0        
GRID         274       0 0.31265 0.50024      0.       0        
GRID         275       0 0.31265 0.37518      0.       0        
GRID         276       0 0.31265 0.25012      0.       0        
GRID         277       0 0.31265 0.12506      0.       0        
CQUAD4       370       1     224     225     238     237                
CQUAD4       371       1     225     226     227     238                
CQUAD4       372       1     237     238     239     236                
CQUAD4       373       1     238     227     263     239                
CQUAD4       374       1     236     239     240     235                
CQUAD4       375       1     239     263     262     240                
CQUAD4       376       1     235     240     241     234                
CQUAD4       377       1     240     262     261     241                
CQUAD4       378       1     234     241     232     233                
CQUAD4       379       1     241     261     231     232                
CQUAD4       380       1     242     243     256     255                
CQUAD4       381       1     243     244     271     256                
CQUAD4       382       1     255     256     257     254                
CQUAD4       383       1     256     271     270     257                
CQUAD4       384       1     254     257     258     253                
CQUAD4       385       1     257     270     269     258                
CQUAD4       386       1     253     258     259     252                
CQUAD4       387       1     258     269     248     259                
CQUAD4       388       1     252     259     250     251                
CQUAD4       389       1     259     248     249     250                
CQUAD4       390       2     231     261     274     273                
CQUAD4       391       2     261     262     275     274                
CQUAD4       392       2     262     263     276     275                
CQUAD4       393       2     263     227     277     276                
CQUAD4       394       2     227     226     266     277                
CQUAD4       395       2     273     274     271     244                
CQUAD4       396       2     274     275     270     271                
CQUAD4       397       2     275     276     269     270                
CQUAD4       398       2     276     277     248     269                
CQUAD4       399       2     277     266     249     248                
ENDDATA
