pro TestCHMP1_RATAN
 RefDir='Data\ref_RATAN\'
 ModelFileName='Data\models\model20180930_0900.sav'
 EBTELfileName='Data\ebtel\ebtelDEM.sav'
 LibFileName='RenderGRFF_64.dll'
 OutDir='Data\output\'
 alist=[1.0, 2.0]
 blist=[1.0, 2.0]
 xc=90.0
 yc=-260.0
 dx=2.0
 dy=2.0
 Nx=120
 Ny=120
 
 MultiScanAB, RefDir, ModelFileName, EBTELfileName, LibFileName, OutDir, $
              alist, blist, xc, yc, dx, dy, Nx, Ny, ObsDateTime='20180930_0900'
end