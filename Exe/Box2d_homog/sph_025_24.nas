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
$   From Model : C:\Users\Dima\Lame3d2\Mpls_all\Exe\Box2d_homog\Solid\model3.MOD
$   Date       : Thu Dec 22 18:36:56 2005
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
$ FEMAP Property 1 : Matrix
PSHELL         1       1      0.       1               1              0.
$ FEMAP Property 2 : Fulleren
PSHELL         2       1      0.       1               1              0.
$ FEMAP Material 1 : Untitled
MAT1           1     7.8     2.3     0.3      0.      0.      0.        
GRID        2579       0 0.32322 0.67678      0.       0        
GRID        2585       0      0.      1.      0.       0        
GRID        2840       0    0.25     0.5      0.       0        
GRID        2841       0      0.     0.5      0.       0        
GRID        2849       0 0.67678 0.67678      0.       0        
GRID        2850       0      1.      1.      0.       0        
GRID        2851       0     0.5    0.75      0.       0        
GRID        2852       0     0.5      1.      0.       0        
GRID        2855       0     0.5     0.5      0.       0        
GRID        2856       0 0.32322 0.32322      0.       0        
GRID        2857       0      0.      0.      0.       0        
GRID        2863       0 0.67678 0.32322      0.       0        
GRID        2864       0      1.      0.      0.       0        
GRID        2865       0     0.5    0.25      0.       0        
GRID        2866       0     0.5      0.      0.       0        
GRID        2867       0    0.75     0.5      0.       0        
GRID        2868       0      1.     0.5      0.       0        
CTRIA3       267       2    2840    2855    2579                        
CTRIA3       268       2    2851    2579    2855                        
CTRIA3       303       1    2585    2579    2852                        
CTRIA3       304       1    2851    2852    2579                        
CTRIA3       305       1    2840    2579    2841                        
CTRIA3       306       1    2841    2579    2585                        
CTRIA3       307       2    2855    2867    2849                        
CTRIA3       308       2    2851    2855    2849                        
CTRIA3       309       1    2849    2850    2852                        
CTRIA3       310       1    2852    2851    2849                        
CTRIA3       311       1    2849    2867    2868                        
CTRIA3       312       1    2849    2868    2850                        
CTRIA3       313       2    2855    2840    2856                        
CTRIA3       314       2    2865    2855    2856                        
CTRIA3       315       1    2856    2857    2866                        
CTRIA3       316       1    2866    2865    2856                        
CTRIA3       317       1    2856    2840    2841                        
CTRIA3       318       1    2856    2841    2857                        
CTRIA3       319       2    2867    2855    2863                        
CTRIA3       320       2    2865    2863    2855                        
CTRIA3       321       1    2864    2863    2866                        
CTRIA3       322       1    2865    2866    2863                        
CTRIA3       323       1    2867    2863    2868                        
CTRIA3       324       1    2868    2863    2864                        
ENDDATA
