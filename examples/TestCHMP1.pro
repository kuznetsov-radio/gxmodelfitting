pro TestCHMP1
 RefDir='C:\MCloud\CoronalMW\AR-SRH\Data\SRH\20230226_034000\'
 ModelFileName='C:\MCloud\CoronalMW\AR-SRH\Data\models\model20230226_0340.sav'
 EBTELfileName='C:\MCloud\CoronalMW\AR-SRH\Data\ebtel\ebtelDEM.sav'
 LibFileName='C:\MCloud\CoronalMW\AR-SRH\RenderX\x64\Release\RenderX_64.dll'
 OutDir='C:\MCloud\CoronalMW\AR-SRH\Data\Images\chmp\'
 alist=[1.00]; [1.00, 1.25]
 blist=[1.25]; [1.25, 1.50]
 xc=-25.0
 yc=505.0
 dx=2.0
 dy=2.0
 Nx=150
 Ny=150
 
 MultiScanAB, RefDir, ModelFileName, EBTELfileName, LibFileName, OutDir, $
              alist, blist, xc, yc, dx, dy, Nx, Ny, ObsDateTime='20230226_0340'
end