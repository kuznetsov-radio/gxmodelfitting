pro TestCHMP2
 RefFileName='C:\MCloud\CoronalMW\AR-SRH\Data\SRH\20220109_044500\SRH0612_CH00_5.80GHz.sav'
 ModelFileName='C:\MCloud\CoronalMW\AR-SRH\Data\models\model20220109_0445.sav'
 EBTELfileName='C:\MCloud\CoronalMW\AR-SRH\Data\ebtel\ebtelDEM.sav'
 LibFileName='C:\MCloud\CoronalMW\AR-SRH\RenderX\x64\Release\RenderX_64.dll'
 OutDir='C:\MCloud\CoronalMW\AR-SRH\Data\Images\chmp\'

 xc=55.0
 yc=-450.0
 dx=2.0
 dy=2.0
 Nx=150
 Ny=150
 
 a_start=2.0
 b_start=2.5
 da=0.25
 db=0.25
 
 SearchForLocalMinimumAB, RefFileName, ModelFileName, EBTELfileName, LibFileName, OutDir, $
                          xc, yc, dx, dy, Nx, Ny, $
                          a_start, b_start, da, db, $
                          ObsDateTime='20220109_0445', ObsFreq='05.80GHz'
end