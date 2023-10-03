pro TestCHMP2
 RefFileName='Data\ref\SRH0612_CH00_5.80GHz.sav'
 ModelFileName='Data\models\model20230226_0340.sav'
 EBTELfileName='Data\ebtel\ebtelDEM.sav'
 LibFileName='RenderGRFF_64.dll'
 OutDir='Data\output\'

 xc=-25.0
 yc=505.0
 dx=2.0
 dy=2.0
 Nx=150
 Ny=150
 
 a_start=2.0
 b_start=1.5
 da=0.25
 db=0.25
 
 SearchForLocalMinimumAB, RefFileName, ModelFileName, EBTELfileName, LibFileName, OutDir, $
                          xc, yc, dx, dy, Nx, Ny, $
                          a_start, b_start, da, db, $
                          ObsDateTime='20230226_0340', ObsFreq='05.80GHz'
end