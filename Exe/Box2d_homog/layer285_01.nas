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
$   From Model : C:\Users\Dima\Lame3d2\IDENT_9november2006\Ident\Exe\box2d_fulleren\Solid\model.MOD
$   Date       : Fri Nov 10 20:27:48 2006
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
GRID         392       0      0.      0.      0.       0        
GRID         394       0 0.43254  0.3204      0.       0        
GRID         396       0 0.43254  0.9612      0.       0        
GRID         397       0      0.  0.9612      0.       0        
GRID         398       0      0.  0.6408      0.       0        
GRID         399       0      0.  0.3204      0.       0        
GRID         400       0  0.9612  0.9612      0.       0        
GRID         402       0 0.52866  0.6408      0.       0        
GRID         403       0 0.52866  0.3204      0.       0        
GRID         404       0 0.52866      0.      0.       0        
GRID         405       0  0.9612      0.      0.       0        
GRID         406       0  0.9612  0.3204      0.       0        
GRID         407       0  0.9612  0.6408      0.       0        
GRID         409       0 0.43254  0.6408      0.       0        
GRID         411       0 0.43254      0.      0.       0        
GRID         415       0 0.52866  0.9612      0.       0        
CTRIA3       547       1     392     394     399                        
CTRIA3       548       1     392     411     394                        
CTRIA3       549       1     399     394     398                        
CTRIA3       550       1     394     409     398                        
CTRIA3       551       1     398     396     397                        
CTRIA3       552       1     398     409     396                        
CTRIA3       553       1     400     402     407                        
CTRIA3       554       1     400     415     402                        
CTRIA3       555       1     407     402     406                        
CTRIA3       556       1     402     403     406                        
CTRIA3       557       1     406     404     405                        
CTRIA3       558       1     406     403     404                        
CTRIA3       559       2     396     402     415                        
CTRIA3       560       2     396     409     402                        
CTRIA3       561       2     409     394     402                        
CTRIA3       562       2     394     403     402                        
CTRIA3       563       2     394     404     403                        
CTRIA3       564       2     394     411     404                        
ENDDATA
