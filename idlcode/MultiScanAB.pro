pro MultiScanAB, RefDir, ModelFileName, EBTELfileName, LibFileName, OutDir, $
                 alist, blist, xc, yc, dx, dy, Nx, Ny, $
                 RefFiles=RefFiles, Q0start=Q0start, threshold=threshold, metric=metric, $
                 MultiThermal=MultiThermal, ObsDateTime=ObsDateTime
;This program searches for the heating rate value Q0 that provides the best agreement between the model and
;observed radio maps, for the specified parameters a and b of the coronal heating model.
;
;Input parameters:
; RefDir - the directory where the observed radio maps are stored.
; If the parameter RefFiles is omitted, the program loads all *.sav files in the RefDir directory.
; Otherwise, the program loads the file(s) specified by RefDir+RefFiles.
; Each .sav file should contain a 'ref' map object with three maps:
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
; alist, blist - a and b parameters of the coronal heating models (scalars or 1D arrays).
;
; xc, yc - x and y coordinates of the model map center, in arcseconds.
;
; dx, dy - x and y resolutions (pixel sizes) of the model map, in arcseconds.
;
; Nx, Ny - x and y sizes of the model map, in pixels.
;
;The optional parameters include:
; RefFiles - if specified, the program loads the observed radio maps from the files given by RefDir+RefFiles.
; Default: all *.sav files in the RefDir directory.
;
; Q0start - the initial value of Q0. It can be either:
; a scalar value (applied to all a an b), or
; a 2D N_a*N_b array, where N_a and N_b are the sizes of alist and blist arrays.
; Default: the best-fit Q0 parameters for AR 12924 (extrapolated to the specified a and b) will be used.
;
; threshold - the threshold value to compute the image mask.
; Comparison of the model and observed radio maps is performed in the area where
; (I_obs gt threshold*max(I_obs)) || (I_mod gt threshold*max(I_mod))
; Default: 0.1.
;
; metric - the metric to minimize. It can be one of the following three options:
; 'rho': rho^2=mean(((I_obs-I_mod)/I_obs)^2),
; 'chi': chi^2=mean(((I_obs-I_mod)/sigma)^2),
; 'eta': eta^2=mean(((I_obs-I_mod)/mean(I_obs))^2).
; All averagings are performed over the above-mentioned sub-region determined by the threshold.
; Default: 'eta'.
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
;Results:
; As the result, for each (a, b) combination the program creates in the OutDir directory a .sav file
; with the name starting with 'fit' and including the used metric, threshold, indicator of the multithermal 
; approach, a and b values, and (optionally) 
; the ObsDateTime string.
; These .sav files contain the following fields:
;  freqList - array of the emission frequencies, in GHz.
;  bestQarr - array of the obtained best-fit heating rates Q0 at different frequencies.
;  rhoArr / chiArr / etaArr - array of the obtained best (minimum) rho^2 / chi^2 / eta^2 metrics at 
;                             different frequencies.
;  modImageArr - (multi-frequency) map object containing the best-fit model radio maps (corresponding to the
;                best-fit heating rates Q0) at different frequencies. The maps are not convolved with the 
;                instrument beam.
;  modImageConvArr - (multi-frequency) map object containing the above-mentioned best-fit model radio maps
;                    convolved with the instrument beam.
;  obsImageArr - (multi-frequency) map object containing the observed radio maps rebinned and shifted to
;                match the best-fit model maps at the corresponding frequencies.
; If the algorithm failed to find the best-fit heating rate at a certain frequency (e.g., the used metric has 
; no minimum within the valid Q0 range, or has more than one local minimum), the corresponding BestQarr and
; rhoArr / chiArr / etaArr are set to NaN, and the corresponding image maps contain all zeros.
; Note: the program does not overwrite the existing fit*.sav files. If the program is interrupted, on the next launch
; it will compute the results only for those (a, b) values that have not been processed before.
; 
; Also, the program creates in the OutDir directory a .sav file with the name starting with 'Summary' and 
; including the used metric, threshold, indicator of the multithermal approach, and (optionally) the 
; ObsDateTime string. This .sav file contains the following fields:
;  alist, blist - the input alist and blist parameters.
;  freqList - array of the emission frequencies, in GHz.
;  bestQ - 3D array (N_a*N_b*N_freq, where N_a, N_b, and N_freq are the sizes of the alist, blist, and freqList
;          arrays, respectively) of the obtained best-fit heating rates Q0 at different values of a, b, and frequency.
;  Iobs, Imod - 3D arrays (N_a*N_b*N_freq) of the total observed and model radio fluxes at different values of
;               a, b, and frequency. The fluxes correspond to the obtained best-fit Q0 values.
;  CC - 3D array (N_a*N_b*N_freq) of the correlation coefficients of the observed and model radio maps at 
;       different values of a, b, and frequency. The coefficients correspond to the obtained best-fit Q0 values.
;  shiftX, shiftY - 3D arrays (N_a*N_b*N_freq) of the shifts (in arcseconds) applied to the observed radio maps
;                   to obtain the best correlation with the model maps, at different values of a, b, and frequency.
;                   The shifts correspond to the obtained best-fit Q0 values.
;  rho / chi / eta - 3D arrays (N_a*N_b*N_freq) of the obtained best (minimum) rho^2 / chi^2 / eta^2 metrics 
;                    at different values of a, b, and frequency.
;  rhoVar / chiVar / etaVar - 3D arrays (N_a*N_b*N_freq) of the shifted metrics defined as:
;                             rhoVar=variance((I_obs-I_mod)/I_obs),
;                             chiVar=variance((I_obs-I_mod)/sigma),
;                             etaVar=variance((I_obs-I_mod)/mean(I_obs)).
;                             The shifted metrics (at different values of a, b, and frequency) correspond to the
;                             obtained best-fit Q0 values, i.e., to the minima of the non-shifted rho / chi / eta
;                             metrics. The shifted metrics are not used for finding the best-fit heating rate.
; If the Summary*.sav file exists, it will be overwritten.

 if exist(RefFiles) then ObsFileNames=RefDir+RefFiles else ObsFileNames=file_search(RefDir+'*.sav')
 LoadObservations, ObsFileNames, obsImaps, obsSImaps, obsInfo

 model=LoadGXmodel(ModelFileName)
 
 ebtel=LoadEBTEL(EBTELfileName)
 
 N_a=n_elements(alist)
 N_b=n_elements(blist)
 
 if exist(Q0start) then begin
  s=size(Q0start, /dimensions)
  if ~((n_elements(s) eq 2) && (s[0] eq N_a) && (s[1] eq N_b)) then begin
   Q0=Q0start[0]
   Q0start=dblarr(N_a, N_b)
   Q0start[*]=Q0
  endif
 endif else begin
  Q0start=dblarr(N_a, N_b)
  for i=0, N_a-1 do for j=0, N_b-1 do Q0start[i, j]=exp(-10.1-0.193*alist[i]+2.17*blist[j])
 endelse
 
 if ~exist(threshold) then threshold=0.1d0
 
 if ~exist(metric) then metric='eta'
 if (metric ne 'eta') && (metric ne 'chi') && (metric ne 'rho') then metric='eta'
 
 if ~exist(MultiThermal) then iso=1 else iso=MultiThermal eq 0

 ObsID=exist(ObsDateTime) ? ObsDateTime : ''
 
 if exist(ObsDateTime) then ObsDateTime='_'+ObsDateTime else ObsDateTime=''

 simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, ObsInfo.freq)       
 
 tstart0=systime(1)
 
 for i=0, N_a-1 do for j=0, N_b-1 do begin
  a=alist[i]
  b=blist[j]
  Qstart=Q0start[i, j]
  
  print, 'Computing the best fit Q for a=', a, ', b=', b
  
  fname=OutDir+'fit_'+metric+'_thr'+string(threshold, format='(F5.3)')+$
        (iso ? '_I' : '_M')+ObsDateTime+$
        '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav' 
  
  if ~file_exist(fname) then begin
   tstart1=systime(1)
   
   FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, a, b, Qstart, iso, $
                 bestQarr, etaArr, etaVarArr, IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                 freqList, allQ, allEta, modImageConvArr, obsImageArr, thr=threshold, metric=metric 
                 
   if metric eq 'eta' then begin              
    save, a, b, bestQarr, etaArr, etaVarArr, modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
          freqList, allQ, allEta, modImageConvArr, obsImageArr, modelFileName, EBTELfileName, $
          filename=fname, /compress
   endif else if metric eq 'chi' then begin
    chiArr=etaArr
    chiVarArr=etaVarArr
    allChi=allEta       
    save, a, b, bestQarr, chiArr, chiVarArr, modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
          freqList, allQ, allChi, modImageConvArr, obsImageArr, modelFileName, EBTELfileName, $
          filename=fname, /compress
   endif else if metric eq 'rho' then begin
    rhoArr=etaArr
    rhoVarArr=etaVarArr
    allRho=allEta       
    save, a, b, bestQarr, rhoArr, rhoVarArr, modImageArr, modFlagArr, IobsArr, ImodArr, CCarr, $
          freqList, allQ, allRho, modImageConvArr, obsImageArr, modelFileName, EBTELfileName, $
          filename=fname, /compress    
   endif       
   
   print, 'The best fit Q for a=', a, ', b=', b, ' computed in ', systime(1)-tstart1, ' s'
   print, 'Total elapsed time: ', systime(1)-tstart0, ' s'
  endif
 endfor
 
 print, 'Creating the summary file'
 
 freqList=obsInfo.freq
 N_freq=n_elements(freqList)
 
 bestQ=dblarr(N_a, N_b, N_freq)
 eta=dblarr(N_a, N_b, N_freq)
 etaVar=dblarr(N_a, N_b, N_freq)
 Iobs=dblarr(N_a, N_b, N_freq)
 Imod=dblarr(N_a, N_b, N_freq)
 CC=dblarr(N_a, N_b, N_freq)    
 shiftX=dblarr(N_a, N_b, N_freq)
 shiftY=dblarr(N_a, N_b, N_freq)
 
 for i=0, N_a-1 do for j=0, N_b-1 do begin
  a=alist[i]
  b=blist[j]
  
  fname=OutDir+'fit_'+metric+'_thr'+string(threshold, format='(F5.3)')+$
        (iso ? '_I': '_M')+ObsDateTime+$
        '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav' 
  o=obj_new('IDL_Savefile', fname)
  o->restore, (metric eq 'eta') ? 'etaArr' : ((metric eq 'chi') ? 'chiArr' : 'rhoArr')
  o->restore, (metric eq 'eta') ? 'etaVarArr' : ((metric eq 'chi') ? 'chiVarArr' : 'rhoVarArr')
  o->restore, 'bestQarr'
  o->restore, 'IobsArr'
  o->restore, 'ImodArr'
  o->restore, 'CCarr'
  obj_destroy, o      
  
  for k=0, N_freq-1 do begin
   eta[i, j, k]=(metric eq 'eta') ? etaArr[k] : ((metric eq 'chi') ? chiArr[k] : rhoArr[k]) 
   etaVar[i, j, k]=(metric eq 'eta') ? etaVarArr[k] : ((metric eq 'chi') ? chiVarArr[k] : rhoVarArr[k])  
   bestQ[i, j, k]=bestQarr[k]
   Iobs[i, j, k]=IobsArr[k]
   Imod[i, j, k]=ImodArr[k]
   CC[i, j, k]=CCarr[k]
   
   m=obsimagearr.getmap(k)
   if tag_exist(m, 'shiftX') then shiftX[i, j, k]=m.shiftX
   if tag_exist(m, 'shiftY') then shiftY[i, j, k]=m.shiftY
  endfor
 endfor 
 
 fname1=OutDir+'Summary_'+metric+'_thr'+string(threshold, format='(F5.3)')+$
        (iso ? '_I': '_M')+ObsDateTime+'.sav'
 
 if metric eq 'eta' then begin
  save, alist, blist, freqList, bestQ, Iobs, Imod, CC, eta, etaVar, modelFileName, EBTELfileName, ObsID, $
        shiftX, shiftY, filename=fname1
 endif else if metric eq 'chi' then begin
  chi=eta
  chiVar=etaVar
  save, alist, blist, freqList, bestQ, Iobs, Imod, CC, chi, chiVar, modelFileName, EBTELfileName, ObsID, $
        shiftX, shiftY, filename=fname1
 endif else begin
  rho=eta
  rhoVar=etaVar
  save, alist, blist, freqList, bestQ, Iobs, Imod, CC, rho, rhoVar, modelFileName, EBTELfileName, ObsID, $
        shiftX, shiftY, filename=fname1  
 endelse
 
 print, 'Done'
end