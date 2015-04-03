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
$   Date       : Fri Nov 17 15:01:36 2006
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
GRID         350       0      0.      0.      0.       0        
GRID         351       0  0.1578      0.      0.       0        
GRID         356       0  0.3156  0.6312      0.       0        
GRID         358       0  0.1578   0.789      0.       0        
GRID         359       0      0.   0.789      0.       0        
GRID         360       0      0.  0.6312      0.       0        
GRID         361       0      0.  0.4734      0.       0        
GRID         362       0      0.  0.3156      0.       0        
GRID         363       0      0.  0.1578      0.       0        
GRID         364       0  0.1578  0.1578      0.       0        
GRID         365       0  0.1578  0.3156      0.       0        
GRID         366       0  0.1578  0.4734      0.       0        
GRID         367       0  0.1578  0.6312      0.       0        
GRID         368       0   0.789   0.789      0.       0        
GRID         369       0  0.6312   0.789      0.       0        
GRID         370       0  0.4734   0.789      0.       0        
GRID         372       0  0.4734  0.4734      0.       0        
GRID         373       0  0.4734  0.3156      0.       0        
GRID         374       0  0.4734  0.1578      0.       0        
GRID         375       0  0.4734      0.      0.       0        
GRID         376       0  0.6312      0.      0.       0        
GRID         377       0   0.789      0.      0.       0        
GRID         378       0   0.789  0.1578      0.       0        
GRID         379       0   0.789  0.3156      0.       0        
GRID         380       0   0.789  0.4734      0.       0        
GRID         381       0   0.789  0.6312      0.       0        
GRID         382       0  0.6312  0.6312      0.       0        
GRID         383       0  0.6312  0.4734      0.       0        
GRID         384       0  0.6312  0.3156      0.       0        
GRID         385       0  0.6312  0.1578      0.       0        
GRID         386       0  0.3156   0.789      0.       0        
GRID         388       0  0.3156  0.4734      0.       0        
GRID         389       0  0.3156  0.3156      0.       0        
GRID         390       0  0.3156  0.1578      0.       0        
GRID         391       0  0.3156      0.      0.       0        
GRID         396       0  0.4734  0.6312      0.       0        
CQUAD4       430       1     350     351     364     363                
CQUAD4       431       1     351     391     390     364                
CQUAD4       432       1     363     364     365     362                
CQUAD4       433       1     364     390     389     365                
CQUAD4       434       1     362     365     366     361                
CQUAD4       435       1     365     389     388     366                
CQUAD4       436       1     361     366     367     360                
CQUAD4       437       1     366     388     356     367                
CQUAD4       438       1     360     367     358     359                
CQUAD4       439       1     367     356     386     358                
CQUAD4       440       1     368     369     382     381                
CQUAD4       441       1     369     370     396     382                
CQUAD4       442       1     381     382     383     380                
CQUAD4       443       1     382     396     372     383                
CQUAD4       444       1     380     383     384     379                
CQUAD4       445       1     383     372     373     384                
CQUAD4       446       1     379     384     385     378                
CQUAD4       447       1     384     373     374     385                
CQUAD4       448       1     378     385     376     377                
CQUAD4       449       1     385     374     375     376                
CQUAD4       450       2     386     356     396     370                
CQUAD4       451       2     356     388     372     396                
CQUAD4       452       2     388     389     373     372                
CQUAD4       453       2     389     390     374     373                
CQUAD4       454       2     390     391     375     374                
ENDDATA
