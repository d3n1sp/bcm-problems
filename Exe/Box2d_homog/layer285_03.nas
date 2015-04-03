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
$   Date       : Fri Nov 17 15:07:18 2006
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
GRID         398       0      0.      0.      0.       0        
GRID         399       0 0.12024      0.      0.       0        
GRID         400       0 0.24049      0.      0.       0        
GRID         401       0 0.24049 0.13742      0.       0        
GRID         402       0 0.24049 0.27484      0.       0        
GRID         404       0 0.24049 0.54968      0.       0        
GRID         406       0 0.12024  0.6871      0.       0        
GRID         407       0      0.  0.6871      0.       0        
GRID         408       0      0. 0.54968      0.       0        
GRID         409       0      0. 0.41226      0.       0        
GRID         410       0      0. 0.27484      0.       0        
GRID         411       0      0. 0.13742      0.       0        
GRID         412       0 0.12024 0.13742      0.       0        
GRID         413       0 0.12024 0.27484      0.       0        
GRID         414       0 0.12024 0.41226      0.       0        
GRID         415       0 0.12024 0.54968      0.       0        
GRID         416       0  0.6871  0.6871      0.       0        
GRID         417       0 0.56686  0.6871      0.       0        
GRID         418       0 0.44661  0.6871      0.       0        
GRID         422       0 0.44661 0.13742      0.       0        
GRID         423       0 0.44661      0.      0.       0        
GRID         424       0 0.56686      0.      0.       0        
GRID         425       0  0.6871      0.      0.       0        
GRID         426       0  0.6871 0.13742      0.       0        
GRID         427       0  0.6871 0.27484      0.       0        
GRID         428       0  0.6871 0.41226      0.       0        
GRID         429       0  0.6871 0.54968      0.       0        
GRID         430       0 0.56686 0.54968      0.       0        
GRID         431       0 0.56686 0.41226      0.       0        
GRID         432       0 0.56686 0.27484      0.       0        
GRID         433       0 0.56686 0.13742      0.       0        
GRID         434       0 0.24049  0.6871      0.       0        
GRID         436       0 0.24049 0.41226      0.       0        
GRID         440       0 0.34355      0.      0.       0        
GRID         443       0 0.44661 0.27484      0.       0        
GRID         444       0 0.44661 0.41226      0.       0        
GRID         445       0 0.44661 0.54968      0.       0        
GRID         447       0 0.34355  0.6871      0.       0        
GRID         448       0 0.34355 0.54968      0.       0        
GRID         449       0 0.34355 0.41226      0.       0        
GRID         450       0 0.34355 0.27484      0.       0        
GRID         451       0 0.34355 0.13742      0.       0        
CQUAD4       455       1     398     399     412     411                
CQUAD4       456       1     399     400     401     412                
CQUAD4       457       1     411     412     413     410                
CQUAD4       458       1     412     401     402     413                
CQUAD4       459       1     410     413     414     409                
CQUAD4       460       1     413     402     436     414                
CQUAD4       461       1     409     414     415     408                
CQUAD4       462       1     414     436     404     415                
CQUAD4       463       1     408     415     406     407                
CQUAD4       464       1     415     404     434     406                
CQUAD4       465       1     416     417     430     429                
CQUAD4       466       1     417     418     445     430                
CQUAD4       467       1     429     430     431     428                
CQUAD4       468       1     430     445     444     431                
CQUAD4       469       1     428     431     432     427                
CQUAD4       470       1     431     444     443     432                
CQUAD4       471       1     427     432     433     426                
CQUAD4       472       1     432     443     422     433                
CQUAD4       473       1     426     433     424     425                
CQUAD4       474       1     433     422     423     424                
CQUAD4       475       2     434     404     448     447                
CQUAD4       476       2     404     436     449     448                
CQUAD4       477       2     436     402     450     449                
CQUAD4       478       2     402     401     451     450                
CQUAD4       479       2     401     400     440     451                
CQUAD4       480       2     447     448     445     418                
CQUAD4       481       2     448     449     444     445                
CQUAD4       482       2     449     450     443     444                
CQUAD4       483       2     450     451     422     443                
CQUAD4       484       2     451     440     423     422                
ENDDATA
