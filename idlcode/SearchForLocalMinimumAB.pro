pro FindMinMetricLocation, mtr_arr, k 
 mtr_arrX=mtr_arr
 u=where(finite(mtr_arrX) and (mtr_arrX lt 0), k)
 if k gt 0 then mtr_arrX[u]=!values.d_NaN
 mtrMin=min(mtr_arrX, k, /NaN)
end  

pro ExpandArrays1, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db, a_range, b_range
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
 
 u=where((a_arr le a_range[0]) or (a_arr ge a_range[1]) or (b_arr le b_range[0]) or (b_arr ge b_range[1]), k)
 if k gt 0 then begin
  mtr_arr[u]=!values.d_NaN
  Q0_arr[u]=!values.d_NaN
 endif
end

pro ExpandArrays2, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db, threshold_metric, a_range, b_range
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
 
 u=where((a_arr le a_range[0]) or (a_arr ge a_range[1]) or (b_arr le b_range[0]) or (b_arr ge b_range[1]), k)
 if k gt 0 then begin
  mtr_arr[u]=!values.d_NaN
  Q0_arr[u]=!values.d_NaN
 endif
end

pro SaveLocalResults, OutDir, ObsDateTime, ObsFreq, $
                      LibFileName, modelFileName, EBTELfileName, DEM_on, DDM_on, $
                      sxArr, syArr, beamArr, $
                      a, b, Qstart, Qstep, iso, threshold_img, threshold_metric, metric, MultiFreq_on, fixed_shifts, $
                      freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ 
                      ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                      obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                      modFlagArr, allQ, allMetrics
 fname=OutDir+'fit_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
       (iso ? '_I' : '_M')+ObsDateTime+ObsFreq+$
       '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav'
        
 save, LibFileName, modelFileName, EBTELfileName, DEM_on, DDM_on, $
       sxArr, syArr, beamArr, $
       a, b, Qstart, Qstep, iso, threshold_img, threshold_metric, metric, MultiFreq_on, fixed_shifts, $
       freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ 
       ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
       obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
       modFlagArr, allQ, allMetrics, $
       filename=fname, /compress
end 

function InSav, o, name
 s=o->Names()
 res=0
 for i=0, n_elements(s)-1 do if strcmp(s[i], name, /fold_case) then res=1
 return, res
end

function LoadLocalResults, OutDir, metric, threshold, iso, ObsDateTime, ObsFreq, a, b, $
         bestQarr, chiArr, rhoArr, etaArr
 fname=OutDir+'fit_'+metric+'_thr'+string(threshold, format='(F5.3)')+$
       (iso ? '_I' : '_M')+ObsDateTime+ObsFreq+$
       '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav'

 if file_exist(fname) then begin
  o=obj_new('IDL_Savefile', fname)
  o->restore, 'bestQarr'
  if InSav(o, 'chiArr') then o->restore, 'chiArr' else chiArr=dblarr(n_elements(bestQarr))
  if InSav(o, 'rhoArr') then o->restore, 'rhoArr' else rhoArr=dblarr(n_elements(bestQarr))
  if InSav(o, 'etaArr') then o->restore, 'etaArr' else etaArr=dblarr(n_elements(bestQarr))
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
                             MultiThermal=MultiThermal, ObsDateTime=ObsDateTime, ObsFreq=ObsFreq, DEM=DEM, DDM=DDM, $
                             a_range=a_range, b_range=b_range, noArea=noArea, Qstep=Qstep, xy_shift=xy_shift, loud=loud
;This program searches for the parameters of the coronal heating model (a, b, Q0) that provide the best agreement 
;between the model and observed radio maps. The search provides a local minimum of the selected model-to-observations
;comparison metric. 
;The program also determines (optionally) the region of "good agreement" in the (a, b) space, 
;where the model-to-observations comparison metric is below a certain threshold (relative to the minimum one). 
;Potentially, if the metric has several local minima, this step may alter the best-fit parameters found at the 
;previous step.
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
; DEM, DDM - these keywords are only applicable if the chosen EBTELfileName .sav file contains both the DEM and DDM tables. 
; In this case, if the /DEM keyword is set, the code loads the DEM table only (the DDM table is ignored). 
; Similarly, if the /DDM keyword is set, the code loads the DDM table only (the DEM table is ignored). 
; If both /DEM and /DDM keywords (or none of them) are set, the code loads both tables.
; If the chosen file contains only one EBTEL table (either DEM or DDM), the code loads that table; 
; the /DEM and /DDM keywords are ignored.
;  
; a_range, b_range - the allowed ranges of the a and b indices, 2-element arrays. The search for the best fit is
; performed in the range of a_range[0] < a < a_range[1], b_range[0] < b < b_range[1].
; Default: a_range=[-10, 10], b_range=[-10, 10].
; Note that you cannot extend the a and b ranges beyond the default values, due to the limitations related to the
; file-naming convention.
; 
; noArea - if set, the area of good agreement (within the threshold_metric threshold) is not computed, i.e., the code
; stops after finding the local minimum.
; Default: 0 (the area of good agreement is computed). 
; 
; Qstep - the initial relative step over Q0 to search for the optimal heating rate value (must be >1).
; Default: the golden ratio value (1.6180339). 
; 
; xy_shift - shifts applied to the observed microwave map, a 2-element vector in the form of xy_shift=[dx, dy], in arcseconds.
; If this parameter is not specified (by default), the shifts are computed automatically each time (i.e., for each frequency
; and a, b, and Q0 values) to provide the maximum correlation between the observed and model images. 
; 
; loud - if set, the code displays more detailed information when it fails to find a solution (e.g., when the minimization
;  procedure goes beyond the EBTEL table).
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
;  ItotalObs, ItotalMod - 3D arrays (N_a*N_b*1) of the total observed and model radio fluxes at different values of a and b. 
;                         The fluxes correspond to the obtained best-fit Q0 values.
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
; corresponding bestQ, ItotalMod, CC, rho, chi, eta, shiftX, and shiftY are set to NaN.
; If the data for a certain combination of a and b are missing (because the search for the best-fit Q0 is performed 
; only within a subset of the rectangular area determined by the regular grids alist and blist), the 
; corresponding bestQ, ItotalMod, CC, rho, chi, and eta are set to -1. Therefore, when analyzing the results, only the 
; points where, e.g., bestQ is finite and bestQ>0 should be considered. If the Summary*.sav file exists, it will be 
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
 sxArr=obsInfo.sx
 syArr=obsInfo.sy
 beamArr=obj_new('map')
 for i=0, obsInfo.Nfreq-1 do begin
  beam=(obsInfo.psf)[i]
  m=make_map(beam, xc=0, yc=0, dx=obsInfo.psf_dx[i], dy=obsInfo.psf_dy[i])
  beamArr->setmap, i, m
 endfor 

 model=LoadGXmodel(ModelFileName)
 
 ebtel=LoadEBTEL(EBTELfileName, DEM=DEM, DDM=DDM)
 DEM_on=ebtel.DEM_on
 DDM_on=ebtel.DDM_on 
 
 if ~exist(Q0start) then Q0start=exp(-10.1-0.193*a_start+2.17*b_start)
 Qstart=Q0start
 
 if ~exist(threshold_img) then threshold_img=0.1d0
 
 if ~exist(threshold_metric) then threshold_metric=2.0
 
 if ~exist(metric) then metric='eta'
 if (metric ne 'eta') && (metric ne 'chi') && (metric ne 'rho') then metric='eta'
 
 if ~exist(MultiThermal) then iso=1 else iso=MultiThermal eq 0

 ObsID=exist(ObsDateTime) ? ObsDateTime : ''
 if exist(ObsDateTime) then ObsDateTime1='_'+ObsDateTime else ObsDateTime1=''
 
 if exist(ObsFreq) then ObsFreq1='_'+ObsFreq else ObsFreq1=''
 
 if exist(a_range) then a_range=[a_range[0]>(-9.999), a_range[1]<9.999] else a_range=[-9.999, 9.999]
 if exist(b_range) then b_range=[b_range[0]>(-9.999), b_range[1]<9.999] else b_range=[-9.999, 9.999]
 
 if ~exist(Qstep) then Qstep=(1d0+sqrt(5d0))/2
 
 if exist(xy_shift) then fixed_shifts=1 else begin
  xy_shift=0
  fixed_shifts=0
 endelse  
 
 MultiFreq_on=1

 simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, ObsInfo.freq)   
 
 print, 'Searching for a local minimum'
 
 a_arr1D=double(a_start)
 b_arr1D=double(b_start)
 a_arr=double(a_start)
 b_arr=double(b_start)
 N_a=1
 N_b=1
 
 if ~LoadLocalResults(OutDir, metric, threshold_img, iso, ObsDateTime1, ObsFreq1, a_start, b_start, $
                      bestQarr, chiArr, rhoArr, etaArr) then begin
  FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ 
                a_start, b_start, Qstart, Qstep, iso, threshold_img, metric, MultiFreq_on, fixed_shifts, xy_shift, $         
                freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $        
                ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                modFlagArr, allQ, allMetrics, loud=loud
                 
  SaveLocalResults, OutDir, ObsDateTime1, ObsFreq1, $
                    LibFileName, modelFileName, EBTELfileName, DEM_on, DDM_on, $
                    sxArr, syArr, beamArr, $
                    a_start, b_start, Qstart, Qstep, iso, threshold_img, threshold_metric, metric, MultiFreq_on, fixed_shifts, $
                    freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ 
                    ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                    obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                    modFlagArr, allQ, allMetrics
 endif   
 
 case metric of                
  'chi': mtr_arr=reform(chiArr)
  'rho': mtr_arr=reform(rhoArr)              
  'eta': mtr_arr=reform(etaArr)
 endcase
 
 Q0_arr=reform(bestQarr)
 
 done=~finite(mtr_arr)
 
 while ~done do begin
  ExpandArrays1, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db, a_range, b_range
  
  FindMinMetricLocation, mtr_arr, k
  idx=array_indices(mtr_arr, k)
  i0=idx[0]
  j0=idx[1]
  
  for i=i0-1, i0+1 do for j=j0-1, j0+1 do if finite(mtr_arr[i, j]) && (mtr_arr[i, j] lt 0) then begin
   if ~LoadLocalResults(OutDir, metric, threshold_img, iso, ObsDateTime1, ObsFreq1, a_arr[i, j], b_arr[i, j], $
                        bestQarr, chiArr, rhoArr, etaArr) then begin
    FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ 
                  a_arr[i, j], b_arr[i, j], Q0_arr[i0, j0], Qstep, iso, threshold_img, metric, $
                  MultiFreq_on, fixed_shifts, xy_shift, $         
                  freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $        
                  ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                  obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                  modFlagArr, allQ, allMetrics, loud=loud
                  
    SaveLocalResults, OutDir, ObsDateTime1, ObsFreq1, $
                      LibFileName, modelFileName, EBTELfileName, DEM_on, DDM_on, $
                      sxArr, syArr, beamArr, $
                      a_arr[i, j], b_arr[i, j], Q0_arr[i0, j0], Qstep, iso, threshold_img, threshold_metric, $
                      metric, MultiFreq_on, fixed_shifts, $
                      freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ 
                      ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                      obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                      modFlagArr, allQ, allMetrics
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
 
 if ~keyword_set(noArea) then begin
  print, 'Exploring the area within the threshold'
 
  done=0
 
  while ~done do begin
   done=1
  
   ExpandArrays2, mtr_arr, Q0_arr, a_arr, b_arr, a_arr1D, b_arr1D, N_a, N_b, da, db, threshold_metric, a_range, b_range
  
   FindMinMetricLocation, mtr_arr, k
   mtr_min=mtr_arr[k]
   u=where(finite(mtr_arr) and (mtr_arr gt 0) and (mtr_arr lt (mtr_min*threshold_metric)), n)
  
   for l=0, n-1 do begin
    idx=array_indices(mtr_arr, u[l])
    i0=idx[0]
    j0=idx[1]
   
    for i=i0-1, i0+1 do for j=j0-1, j0+1 do if finite(mtr_arr[i, j]) && (mtr_arr[i, j] lt 0) then begin
     if ~LoadLocalResults(OutDir, metric, threshold_img, iso, ObsDateTime1, ObsFreq1, a_arr[i, j], b_arr[i, j], $
                          bestQarr, chiArr, rhoArr, etaArr) then begin
      FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ 
                    a_arr[i, j], b_arr[i, j], Q0_arr[i0, j0], Qstep, iso, threshold_img, metric, $
                    MultiFreq_on, fixed_shifts, xy_shift, $         
                    freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $        
                    ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                    obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                    modFlagArr, allQ, allMetrics, loud=loud
      SaveLocalResults, OutDir, ObsDateTime1, ObsFreq1, $
                        LibFileName, modelFileName, EBTELfileName, DEM_on, DDM_on, $
                        sxArr, syArr, beamArr, $
                        a_arr[i, j], b_arr[i, j], Q0_arr[i0, j0], Qstep, iso, threshold_img, threshold_metric, $
                        metric, MultiFreq_on, fixed_shifts, $
                        freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ 
                        ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                        obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                        modFlagArr, allQ, allMetrics
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
 endif 
 
 print, 'Creating the summary file'
 
 freqList=obsInfo.freq
 N_freq=n_elements(freqList)
 alist=a_arr1D
 blist=b_arr1D
 
 bestQ=dblarr(N_a, N_b, N_freq)
 chi=dblarr(N_a, N_b, N_freq)
 rho=dblarr(N_a, N_b, N_freq)
 eta=dblarr(N_a, N_b, N_freq)
 ItotalObs=dblarr(N_a, N_b, N_freq)
 ItotalMod=dblarr(N_a, N_b, N_freq)
 CC=dblarr(N_a, N_b, N_freq)    
 shiftX=dblarr(N_a, N_b, N_freq)
 shiftY=dblarr(N_a, N_b, N_freq)
 
 bestQ[*]=-1
 chi[*]=-1
 rho[*]=-1
 eta[*]=-1
 ItotalObs[*]=-1
 ItotalMod[*]=-1
 CC[*]=-1
 shiftX[*]=!values.d_NaN
 shiftY[*]=!values.d_NaN
 
 for i=0, N_a-1 do for j=0, N_b-1 do begin
  a=a_arr1D[i]
  b=b_arr1D[j]
  
  if (a le a_range[0]) || (a ge a_range[1]) || (b le b_range[0]) || (b ge b_range[1]) then begin
   chi[i, j, *]=!values.d_NaN 
   rho[i, j, *]=!values.d_NaN 
   eta[i, j, *]=!values.d_NaN 
   bestQ[i, j, *]=!values.d_NaN 
   ItotalObs[i, j, *]=!values.d_NaN 
   ItotalMod[i, j, *]=!values.d_NaN 
   CC[i, j, *]=!values.d_NaN 
  endif else begin 
   fname=OutDir+'fit_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
         (iso ? '_I' : '_M')+ObsDateTime1+ObsFreq1+$
         '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav'
  
   if file_exist(fname) then begin
    o=obj_new('IDL_Savefile', fname)
    if InSav(o, 'chiArr') then o->restore, 'chiArr' else chiArr=dblarr(N_freq)
    if InSav(o, 'rhoArr') then o->restore, 'rhoArr' else rhoArr=dblarr(N_freq)
    if InSav(o, 'etaArr') then o->restore, 'etaArr' else etaArr=dblarr(N_freq)
    o->restore, 'bestQarr'
    o->restore, 'ItotalObsArr'
    o->restore, 'ItotalModArr'
    o->restore, 'CCarr'
    o->restore, 'obsImageArr'
    obj_destroy, o      
  
    for k=0, N_freq-1 do begin
     chi[i, j, k]=chiArr[k]
     rho[i, j, k]=rhoArr[k]
     eta[i, j, k]=etaArr[k]
     bestQ[i, j, k]=bestQarr[k]
     ItotalObs[i, j, k]=ItotalObsArr[k]
     ItotalMod[i, j, k]=ItotalModArr[k]
     CC[i, j, k]=CCarr[k]
     
     m=obsImageArr.getmap(k)
     if tag_exist(m, 'shiftX') then shiftX[i, j, k]=m.shiftX
     if tag_exist(m, 'shiftY') then shiftY[i, j, k]=m.shiftY
    endfor
   endif
  endelse  
 endfor 
 
 fname1=OutDir+'Summary_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
        (iso ? '_I': '_M')+ObsDateTime1+ObsFreq1+'.sav'
 
 save, alist, blist, freqList, bestQ, ItotalObs, ItotalMod, CC, chi, rho, eta, metric, iso, threshold_img, threshold_metric, $
       modelFileName, EBTELfileName, DEM_on, DDM_on, MultiFreq_on, ObsID, shiftX, shiftY, fixed_shifts, $
       filename=fname1
 
 print, 'Done'
 
 print, 'Elapsed time: ', systime(1)-tstart0, ' s'
end