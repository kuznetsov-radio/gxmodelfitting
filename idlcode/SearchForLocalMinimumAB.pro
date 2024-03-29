pro FindMinMetricLocation, mtr_arr, k 
 mtr_arrX=mtr_arr
 u=where(finite(mtr_arrX) and (mtr_arrX lt 0), k)
 if k gt 0 then mtr_arrX[u]=!values.d_NaN
 mtrMin=min(mtr_arrX, k, /NaN)
end  

pro ExpandArrays1, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db
 FindMinMetricLocation, mtr_arr, k
 if a_arr[k] eq min(a_arr) then begin
  N_a+=1
  a_arr1D=[a_arr1D[0]-da, a_arr1D]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for j=0, N_b-1 do begin
   a_arrX[*, j]=a_arr1D
   b_arrX[*, j]=b_arr1D[j]
   mtr_arrX[*, j]=[-1d0, reform(mtr_arr[*, j])]
   Q0_arrX[*, j]=[-1d0, reform(Q0_arr[*, j])]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif 
  
 FindMinMetricLocation, mtr_arr, k
 if a_arr[k] eq max(a_arr) then begin
  N_a+=1
  a_arr1D=[a_arr1D, a_arr1D[N_a-2]+da]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for j=0, N_b-1 do begin
   a_arrX[*, j]=a_arr1D
   b_arrX[*, j]=b_arr1D[j]
   mtr_arrX[*, j]=[reform(mtr_arr[*, j]), -1d0]
   Q0_arrX[*, j]=[reform(Q0_arr[*, j]), -1d0]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif   
  
 FindMinMetricLocation, mtr_arr, k
 if b_arr[k] eq min(b_arr) then begin
  N_b+=1
  b_arr1D=[b_arr1D[0]-db, b_arr1D]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for i=0, N_a-1 do begin
   a_arrX[i, *]=a_arr1D[i]
   b_arrX[i, *]=b_arr1D
   mtr_arrX[i, *]=[-1d0, reform(mtr_arr[i, *])]
   Q0_arrX[i, *]=[-1d0, reform(Q0_arr[i, *])]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif 
  
 FindMinMetricLocation, mtr_arr, k
 if b_arr[k] eq max(b_arr) then begin
  N_b+=1
  b_arr1D=[b_arr1D, b_arr1D[N_b-2]+db]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for i=0, N_a-1 do begin
   a_arrX[i, *]=a_arr1D[i]
   b_arrX[i, *]=b_arr1D
   mtr_arrX[i, *]=[reform(mtr_arr[i, *]), -1d0]
   Q0_arrX[i, *]=[reform(Q0_arr[i, *]), -1d0]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif   
 
 u=where((abs(a_arr) ge 10.0) or (abs(b_arr) ge 10.0), k)
 if k gt 0 then begin
  mtr_arr[u]=!values.d_NaN
  Q0_arr[u]=!values.d_NaN
 endif
end

pro ExpandArrays2, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db, threshold_metric
 FindMinMetricLocation, mtr_arr, k
 mtr_min=mtr_arr[k]
 u=where(finite(mtr_arr) and (mtr_arr gt 0) and (mtr_arr lt (mtr_min*threshold_metric)))
 
 if min(a_arr[u]) eq min(a_arr) then begin
  N_a+=1
  a_arr1D=[a_arr1D[0]-da, a_arr1D]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for j=0, N_b-1 do begin
   a_arrX[*, j]=a_arr1D
   b_arrX[*, j]=b_arr1D[j]
   mtr_arrX[*, j]=[-1d0, reform(mtr_arr[*, j])]
   Q0_arrX[*, j]=[-1d0, reform(Q0_arr[*, j])]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif 
  
 FindMinMetricLocation, mtr_arr, k
 mtr_min=mtr_arr[k]
 u=where(finite(mtr_arr) and (mtr_arr gt 0) and (mtr_arr lt (mtr_min*threshold_metric)))

 if max(a_arr[u]) eq max(a_arr) then begin
  N_a+=1
  a_arr1D=[a_arr1D, a_arr1D[N_a-2]+da]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for j=0, N_b-1 do begin
   a_arrX[*, j]=a_arr1D
   b_arrX[*, j]=b_arr1D[j]
   mtr_arrX[*, j]=[reform(mtr_arr[*, j]), -1d0]
   Q0_arrX[*, j]=[reform(Q0_arr[*, j]), -1d0]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif   
  
 FindMinMetricLocation, mtr_arr, k
 mtr_min=mtr_arr[k]
 u=where(finite(mtr_arr) and (mtr_arr gt 0) and (mtr_arr lt (mtr_min*threshold_metric)))

 if min(b_arr[u]) eq min(b_arr) then begin
  N_b+=1
  b_arr1D=[b_arr1D[0]-db, b_arr1D]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for i=0, N_a-1 do begin
   a_arrX[i, *]=a_arr1D[i]
   b_arrX[i, *]=b_arr1D
   mtr_arrX[i, *]=[-1d0, reform(mtr_arr[i, *])]
   Q0_arrX[i, *]=[-1d0, reform(Q0_arr[i, *])]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif 
  
 FindMinMetricLocation, mtr_arr, k
 mtr_min=mtr_arr[k]
 u=where(finite(mtr_arr) and (mtr_arr gt 0) and (mtr_arr lt (mtr_min*threshold_metric)))

 if max(b_arr[u]) eq max(b_arr) then begin
  N_b+=1
  b_arr1D=[b_arr1D, b_arr1D[N_b-2]+db]
  a_arrX=dblarr(N_a, N_b)
  b_arrX=dblarr(N_a, N_b)
  mtr_arrX=dblarr(N_a, N_b)
  Q0_arrX=dblarr(N_a, N_b)
  for i=0, N_a-1 do begin
   a_arrX[i, *]=a_arr1D[i]
   b_arrX[i, *]=b_arr1D
   mtr_arrX[i, *]=[reform(mtr_arr[i, *]), -1d0]
   Q0_arrX[i, *]=[reform(Q0_arr[i, *]), -1d0]
  endfor
  a_arr=a_arrX
  b_arr=b_arrX
  mtr_arr=mtr_arrX
  Q0_arr=Q0_arrX
 endif   
 
 u=where((abs(a_arr) ge 10.0) or (abs(b_arr) ge 10.0), k)
 if k gt 0 then begin
  mtr_arr[u]=!values.d_NaN
  Q0_arr[u]=!values.d_NaN
 endif
end

pro SaveLocalResults, OutDir, metric, threshold, iso, ObsDateTime, ObsFreq, a, b, $
                      bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                      modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
                      freqList, allQ, allMetrics, modImageConvArr, obsImageArr, modelFileName, EBTELfileName
 fname=OutDir+'fit_'+metric+'_thr'+string(threshold, format='(F5.3)')+$
       (iso ? '_I' : '_M')+ObsDateTime+ObsFreq+$
       '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav'
        
 save, a, b, bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
       modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
       freqList, allQ, allMetrics, modImageConvArr, obsImageArr, modelFileName, EBTELfileName, $
       filename=fname, /compress
end 

function LoadLocalResults, OutDir, metric, threshold, iso, ObsDateTime, ObsFreq, a, b, $
         bestQarr, chiArr, rhoArr, etaArr
 fname=OutDir+'fit_'+metric+'_thr'+string(threshold, format='(F5.3)')+$
       (iso ? '_I' : '_M')+ObsDateTime+ObsFreq+$
       '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav'

 if file_exist(fname) then begin
  o=obj_new('IDL_Savefile', fname)
  o->restore, 'bestQarr'
  o->restore, 'chiArr'
  o->restore, 'rhoArr'
  o->restore, 'etaArr' 
  obj_destroy, o
  res=1
 endif else res=0   
 
 return, res
end

pro SearchForLocalMinimumAB, RefFileName, ModelFileName, EBTELfileName, LibFileName, OutDir, $
                             xc, yc, dx, dy, Nx, Ny, $
                             a_start, b_start, da, db, $
                             Q0start=Q0start, metric=metric, $
                             threshold_img=threshold_img, threshold_metric=threshold_metric, $
                             MultiThermal=MultiThermal, ObsDateTime=ObsDateTime, ObsFreq=ObsFreq
;This program searches for the parameters of the coronal heating model (a, b, Q0) that provide the best agreement 
;between the model and observed radio maps. The search provides a local minimum of the selected model-to-observations
;comparison metric. 
;The program also determines the region of "good agreement" in the (a, b) space, where the model-to-observations 
;comparison metric is below a certain threshold (relative to the minimum one). Potentially, if the metric has 
;several local minima, this step may alter the best-fit parameters found at the previous step.
;
;Input parameters:
; RefFileName - name of the .sav file that contains the observed radio maps (at a single frequency).
; The file should contain a 'ref' map object with three maps:
; I_obs=ref.getmap(0) - the observed radio map, with the tags I_obs.freq specifying the emission frequency in GHz,
;                       and I_obs.id specifying the map title,
; sigma=ref.getmap(1) - the corresponding instrumental noise (with the same dimensions as I_obs),
; beam =ref.getmap(2) - the instrument beam (point-spread function), with the tags beam.a_beam and beam.b_beam
;                       specifying the beam half-widths at 1/e level in two ortogonal directions, in arcseconds.
; Other required tags of these maps are standard for the SSW map structure.
;
; ModelFileName - name of the .sav file that contains the GX Simulator model.
;
; EBTELfileName - name of the .sav file that contains the EBTEL table.
;
; LibFileName - name of the appropriate executable library that computes the radio emission
; (see https://github.com/kuznetsov-radio/gximagecomputing).
;
; OutDir - the directory where the results will be stored.
;
; xc, yc - x and y coordinates of the model map center, in arcseconds.
;
; dx, dy - x and y resolutions (pixel sizes) of the model map, in arcseconds.
;
; Nx, Ny - x and y sizes of the model map, in pixels.
;
; a_start, b_start - the starting point of the metric minimization algorithm in the (a, b) space.
;
; da, db - the grid sizes of the metric minimization algorithm in the a and b directions, respectively.
;
;The optional parameters include:
; Q0start - the starting point of the metric minimization algorithm in the Q0 space. This parameter is applied
; at the point (a_start, b_start) only; at other (a, b) points, the best-fit value of Q0 at the nearest adjacent
; point is used.
; Default: the best-fit Q0 parameters for AR 12924 (extrapolated to the specified a_start and b_start) will be used.
;
; metric - the metric to minimize. It can be one of the following three options:
; 'rho': rho^2=mean(((I_obs-I_mod)/I_obs)^2),
; 'chi': chi^2=mean(((I_obs-I_mod)/sigma)^2),
; 'eta': eta^2=mean(((I_obs-I_mod)/mean(I_obs))^2).
; All averagings are performed over the above-mentioned sub-region determined by the threshold.
; Default: 'eta'.
;
; threshold_img - the threshold value to compute the image mask.
; Comparison of the model and observed radio maps is performed in the area where
; (I_obs gt threshold_img*max(I_obs)) || (I_mod gt threshold_img*max(I_mod))
; Default: 0.1.
;
; threshold_metric - the threshold value that determines the area of "good agreement" between the model and
; observations in the (a, b) space. The algorithm will scan the (a, b) space with the da and db resolutions, until
; it finds all points where the condition 
; metric(a, b) < threshold_metric*min(metric) 
; is satisfied.
; Default: 2.
;
; MultiThermal - if set, the formulae for the free-free and gyroresonance emissions from multithermal plasma
; (described by the DEM and DDM, respectively) are used, see Fleishman, Kuznetsov & Landi (2021).
; If not set, the isothermal approximation with the plasma density and temperature derived from the DDM or DEM
; (from the DDM, if both are specified) is used.
; Default: 0 (isothermal approximation is assumed).
;
; ObsDateTime - an additional string added to the names of the resulting files.
; Default: ''
;
; ObsFreq - an additional string added to the names of the resulting files.
; Default: ''
;
;Results:
; The output of the program is similar to that of the MultiScanAB.pro, with the difference that only one frequency
; is considered. The program creates in the OutDir directory a .sav file with the name starting with 'Summary' and 
; including the used metric, threshold, indicator of the multithermal approach, and (optionally) the 
; ObsDateTime and ObsFreq strings. This .sav file contains the following fields:
;  alist, blist - arrays of a and b parameters covering the range searched by the program during the metric 
;                 minimization process, with da and db steps.
;  freqList - 1-element array containing the emission frequency, in GHz.
;  bestQ - 3D array (N_a*N_b*1, where N_a and N_b are the sizes of the alist and blist arrays, respectively) of the 
;          obtained best-fit heating rates Q0 at different values of a and b.
;  Iobs, Imod - 3D arrays (N_a*N_b*1) of the total observed and model radio fluxes at different values of a and b. 
;               The fluxes correspond to the obtained best-fit Q0 values.
;  CC - 3D array (N_a*N_b*1) of the correlation coefficients of the observed and model radio maps at different 
;       values of a and b. The coefficients correspond to the obtained best-fit Q0 values.
;  shiftX, shiftY - 3D arrays (N_a*N_b*1) of the shifts (in arcseconds) applied to the observed radio maps
;                   to obtain the best correlation with the model maps, at different values of a, b, and frequency.
;                   The shifts correspond to the obtained best-fit Q0 values.
;  rho, chi, eta - 3D arrays (N_a*N_b*1) of the obtained rho^2, chi^2, and eta^2 metrics at different values 
;                  of a and b. Note that only one of those metrics (defined by the 'metric' keyword) is actually
;                  minimized; two other metrics correspond to the obtained best-fit Q0 values.
; Note, that if the algorithm has failed to find the best-fit heating rate Q0 at a certain combination of a and b 
; (e.g., the used metric has no minimum within the valid Q0 range, or has more than one local minimum), the 
; corresponding bestQ, Imod, CC, rho, chi, eta, shiftX, and shiftY are set to NaN.
; If the data for a certain combination of a and b are missing (because the search for the best-fit Q0 is performed 
; only within a subset of the rectangular area determined by the regular grids alist and blist), the 
; corresponding bestQ, Imod, CC, rho, chi, and eta are set to -1. Therefore, when  analyzing the results, only the 
; points where bestQ is finite and bestQ>0 should be considered. If the Summary*.sav file exists, it will be 
; overwritten.

; The program saves the temporary progress: for each (a, b) combination it creates in the OutDir directory a .sav 
; file with the name starting with 'fit' and including the used metric, threshold, indicator of the multithermal 
; approach, a and b values, and (optionally) the ObsDateTime and ObsFreq strings.
; These .sav files contain the following fields:
;  freqList - 1-element array containing the emission frequency, in GHz.
;  bestQarr - 1-element array containing the obtained best-fit heating rate Q0.
;  rhoArr, chiArr, etaArr - 1-element arrays containing the obtained rho^2, chi^2, and eta^2 metrics.
;                           Note that only one of those metrics (defined by the 'metric' keyword) is actually
;                           minimized; two other metrics correspond to the obtained best-fit Q0 values.
;  modImageArr - map object containing the best-fit model radio map (corresponding to the best-fit heating rate Q0). 
;                The map is not convolved with the instrument beam.
;  modImageConvArr - map object containing the above-mentioned best-fit model radio map convolved with the 
;                    instrument beam.
;  obsImageArr - map object containing the observed radio map rebinned and shifted to match the best-fit model map.
; If the algorithm failed to find the best-fit heating rate (e.g., the used metric has no minimum within the 
; valid Q0 range, or has more than one local minimum), the corresponding BestQarr, rhoArr, chiArr, and etaArr are 
; set to NaN, and the corresponding image maps contain all zeros.
; Note: the program does not overwrite the existing fit*.sav files. If the program is interrupted, on the next launch
; it will compute the results only for those (a, b) values that have not been processed before.
           
 tstart0=systime(1)    

 LoadObservations, RefFileName, obsImaps, obsSImaps, obsInfo

 model=LoadGXmodel(ModelFileName)
 
 ebtel=LoadEBTEL(EBTELfileName)
 
 if ~exist(Q0start) then Q0start=exp(-10.1-0.193*a_start+2.17*b_start)
 
 if ~exist(threshold_img) then threshold_img=0.1d0
 
 if ~exist(threshold_metric) then threshold_metric=2.0
 
 if ~exist(metric) then metric='eta'
 if (metric ne 'eta') && (metric ne 'chi') && (metric ne 'rho') then metric='eta'
 
 if ~exist(MultiThermal) then iso=1 else iso=MultiThermal eq 0

 ObsID=exist(ObsDateTime) ? ObsDateTime : ''
 if exist(ObsDateTime) then ObsDateTime='_'+ObsDateTime else ObsDateTime=''
 
 if exist(ObsFreq) then ObsFreq='_'+ObsFreq else ObsFreq=''

 simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, ObsInfo.freq)   
 
 print, 'Searching for a local minimum'
 
 a_arr1D=double(a_start)
 b_arr1D=double(b_start)
 a_arr=double(a_start)
 b_arr=double(b_start)
 N_a=1
 N_b=1
 
 if ~LoadLocalResults(OutDir, metric, threshold_img, iso, ObsDateTime, ObsFreq, a_start, b_start, $
                      bestQarr, chiArr, rhoArr, etaArr) then begin
  FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, a_start, b_start, Q0start, iso, $
                bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                freqList, allQ, allMetrics, modImageConvArr, obsImageArr, thr=threshold_img, metric=metric
  SaveLocalResults, OutDir, metric, threshold_img, iso, ObsDateTime, ObsFreq, a_start, b_start, $
                    bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                    modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
                    freqList, allQ, allMetrics, modImageConvArr, obsImageArr, modelFileName, EBTELfileName
 endif   
 
 case metric of                
  'chi': mtr_arr=reform(chiArr)
  'rho': mtr_arr=reform(rhoArr)              
  'eta': mtr_arr=reform(etaArr)
 endcase
 
 Q0_arr=reform(bestQarr)
 
 done=~finite(mtr_arr)
 
 while ~done do begin
  ExpandArrays1, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db
  
  FindMinMetricLocation, mtr_arr, k
  idx=array_indices(mtr_arr, k)
  i0=idx[0]
  j0=idx[1]
  
  for i=i0-1, i0+1 do for j=j0-1, j0+1 do if finite(mtr_arr[i, j]) && (mtr_arr[i, j] lt 0) then begin
   if ~LoadLocalResults(OutDir, metric, threshold_img, iso, ObsDateTime, ObsFreq, a_arr[i, j], b_arr[i, j], $
                        bestQarr, chiArr, rhoArr, etaArr) then begin
    FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $
                  a_arr[i, j], b_arr[i, j], Q0_arr[i0, j0], iso, $
                  bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                  IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                  freqList, allQ, allMetrics, modImageConvArr, obsImageArr, thr=threshold_img, metric=metric
    SaveLocalResults, OutDir, metric, threshold_img, iso, ObsDateTime, ObsFreq, a_arr[i, j], b_arr[i, j], $
                      bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                      modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
                      freqList, allQ, allMetrics, modImageConvArr, obsImageArr, modelFileName, EBTELfileName
   endif       
   
   case metric of                      
    'chi': mtr_arr[i, j]=reform(chiArr)
    'rho': mtr_arr[i, j]=reform(rhoArr)    
    'eta': mtr_arr[i, j]=reform(etaArr)
   endcase
   
   Q0_arr[i, j]=reform(bestQarr)
  endif
  
  done=1
  for i=i0-1, i0+1 do for j=j0-1, j0+1 do $
   if ((i ne i0) || (j ne j0)) && finite(mtr_arr[i, j]) && (mtr_arr[i, j] lt mtr_arr[i0, j0]) then done=0 
 endwhile
 
 print, 'Exploring the area within the threshold'
 
 done=0
 
 while ~done do begin
  done=1
  
  ExpandArrays2, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db, threshold_metric
  
  FindMinMetricLocation, mtr_arr, k
  mtr_min=mtr_arr[k]
  u=where(finite(mtr_arr) and (mtr_arr gt 0) and (mtr_arr lt (mtr_min*threshold_metric)), n)
  
  for l=0, n-1 do begin
   idx=array_indices(mtr_arr, u[l])
   i0=idx[0]
   j0=idx[1]
   
   for i=i0-1, i0+1 do for j=j0-1, j0+1 do if finite(mtr_arr[i, j]) && (mtr_arr[i, j] lt 0) then begin
    if ~LoadLocalResults(OutDir, metric, threshold_img, iso, ObsDateTime, ObsFreq, a_arr[i, j], b_arr[i, j], $
                         bestQarr, chiArr, rhoArr, etaArr) then begin
     FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $
                   a_arr[i, j], b_arr[i, j], Q0_arr[i0, j0], iso, $
                   bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                   IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                   freqList, allQ, allMetrics, modImageConvArr, obsImageArr, thr=threshold_img, metric=metric
     SaveLocalResults, OutDir, metric, threshold_img, iso, ObsDateTime, ObsFreq, a_arr[i, j], b_arr[i, j], $
                       bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                       modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
                       freqList, allQ, allMetrics, modImageConvArr, obsImageArr, modelFileName, EBTELfileName
    endif                                 

    case metric of                      
     'chi': mtr_arr[i, j]=reform(chiArr)
     'rho': mtr_arr[i, j]=reform(rhoArr)    
     'eta': mtr_arr[i, j]=reform(etaArr)
    endcase
    
    Q0_arr[i, j]=reform(bestQarr)
    done=0
   endif
  endfor
 endwhile
 
 print, 'Creating the summary file'
 
 freqList=obsInfo.freq
 N_freq=n_elements(freqList)
 alist=a_arr1D
 blist=b_arr1D
 
 bestQ=dblarr(N_a, N_b, N_freq)
 chi=dblarr(N_a, N_b, N_freq)
 chiVar=dblarr(N_a, N_b, N_freq)
 rho=dblarr(N_a, N_b, N_freq)
 rhoVar=dblarr(N_a, N_b, N_freq)
 eta=dblarr(N_a, N_b, N_freq)
 etaVar=dblarr(N_a, N_b, N_freq)
 Iobs=dblarr(N_a, N_b, N_freq)
 Imod=dblarr(N_a, N_b, N_freq)
 CC=dblarr(N_a, N_b, N_freq)    
 shiftX=dblarr(N_a, N_b, N_freq)
 shiftY=dblarr(N_a, N_b, N_freq)
 
 bestQ[*]=-1
 chi[*]=-1
 chiVar[*]=-1
 rho[*]=-1
 rhoVar[*]=-1 
 eta[*]=-1
 etaVar[*]=-1 
 Iobs[*]=-1
 Imod[*]=-1
 CC[*]=-1
 shiftX[*]=!values.d_NaN
 shiftY[*]=!values.d_NaN
 
 for i=0, N_a-1 do for j=0, N_b-1 do begin
  a=a_arr1D[i]
  b=b_arr1D[j]
  
  if (abs(a) ge 10.0) || (abs(b) ge 10.0) then begin
   chi[i, j, *]=!values.d_NaN 
   chiVar[i, j, *]=!values.d_NaN
   rho[i, j, *]=!values.d_NaN 
   rhoVar[i, j, *]=!values.d_NaN      
   eta[i, j, *]=!values.d_NaN 
   etaVar[i, j, *]=!values.d_NaN   
   bestQ[i, j, *]=!values.d_NaN 
   Iobs[i, j, *]=!values.d_NaN 
   Imod[i, j, *]=!values.d_NaN 
   CC[i, j, *]=!values.d_NaN 
  endif else begin 
   fname=OutDir+'fit_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
         (iso ? '_I' : '_M')+ObsDateTime+ObsFreq+$
         '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav'
  
   if file_exist(fname) then begin
    o=obj_new('IDL_Savefile', fname)
    o->restore, 'chiArr'
    o->restore, 'chiVarArr'
    o->restore, 'rhoArr'
    o->restore, 'rhoVarArr'        
    o->restore, 'etaArr'
    o->restore, 'etaVarArr'
    o->restore, 'bestQarr'
    o->restore, 'IobsArr'
    o->restore, 'ImodArr'
    o->restore, 'CCarr'
    o->restore, 'obsImageArr'
    obj_destroy, o      
  
    for k=0, N_freq-1 do begin
     chi[i, j, k]=chiArr[k]
     chiVar[i, j, k]=chiVarArr[k]
     rho[i, j, k]=rhoArr[k]
     rhoVar[i, j, k]=rhoVarArr[k]         
     eta[i, j, k]=etaArr[k]
     etaVar[i, j, k]=etaVarArr[k]  
     bestQ[i, j, k]=bestQarr[k]
     Iobs[i, j, k]=IobsArr[k]
     Imod[i, j, k]=ImodArr[k]
     CC[i, j, k]=CCarr[k]
     
     m=obsimagearr.getmap(k)
     if tag_exist(m, 'shiftX') then shiftX[i, j, k]=m.shiftX
     if tag_exist(m, 'shiftY') then shiftY[i, j, k]=m.shiftY
    endfor
   endif
  endelse  
 endfor 
 
 fname1=OutDir+'Summary_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
        (iso ? '_I': '_M')+ObsDateTime+ObsFreq+'.sav'
 
 save, alist, blist, freqList, bestQ, Iobs, Imod, CC, chi, chiVar, rho, rhoVar, eta, etaVar, $
       modelFileName, EBTELfileName, ObsID, shiftX, shiftY, filename=fname1
 
 print, 'Done'
 
 print, 'Elapsed time: ', systime(1)-tstart0, ' s'
end