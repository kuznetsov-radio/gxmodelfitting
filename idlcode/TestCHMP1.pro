pro TestCHMP1
 RefDir='Data\ref\'
 ModelFileName='Data\models\model20230226_0340.sav'
 EBTELfileName='Data\ebtel\ebtelDEM.sav'
 LibFileName='RenderGRFF_64.dll'
 OutDir='Data\output\'
 alist=[1.00, 1.25]
 blist=[1.25, 1.50]
 xc=-25.0
 yc=505.0
 dx=2.0
 dy=2.0
 Nx=150
 Ny=150
 
 MultiScanAB, RefDir, ModelFileName, EBTELfileName, LibFileName, OutDir, $
              alist, blist, xc, yc, dx, dy, Nx, Ny, ObsDateTime='20230226_0340'
end