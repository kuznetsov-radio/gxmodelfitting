pro MultiScanAB, RefDir, ModelFileName, EBTELfileName, LibFileName, OutDir, $
                 alist, blist, xc, yc, dx, dy, Nx, Ny, $
                 RefFiles=RefFiles, Q0start=Q0start, threshold=threshold, metric=metric, $
                 MultiThermal=MultiThermal, ObsDateTime=ObsDateTime, noMultiFreq=noMultiFreq, DEM=DEM, DDM=DDM, $
                 Qstep=Qstep, xy_shift=xy_shift, loud=loud, SHtable=SHtable, Nthreads=Nthreads, $
                 analyticalNT=analyticalNT, EMthreshold=EMthreshold
;This program searches for the heating rate value Q0 that provides the best agreement between the model and
;observed radio maps, for the specified parameters a and b of the coronal heating model.
;
;Input parameters:
; RefDir - the directory where the observed radio maps/profiles are stored.
; If the parameter RefFiles is omitted, the program loads all *.sav files in the RefDir directory.
; Otherwise, the program loads the file(s) specified by RefDir+RefFiles.
; For general 2D radio maps, each .sav file should contain a 'ref' map object with three maps:
; I_obs=ref.getmap(0) - the observed radio map (in terms of brightness temperature in K), 
;                       with the tags I_obs.freq specifying the emission frequency in GHz,
;                       and I_obs.id specifying the map title,
; sigma=ref.getmap(1) - the corresponding instrumental noise (with the same dimensions as I_obs),
; beam =ref.getmap(2) - the instrument beam (point-spread function), with the tags beam.a_beam and beam.b_beam
;                       specifying the beam half-widths at 1/e level in two ortogonal directions, in arcseconds.
; Other required tags of these maps are standard for the SSW map structure.
; For 1D intensity profiles observed by RATAN-600, each .sav file should contain two fields:
; instrument='RATAN' - the label to identify the data format,
; ref={flux, x, freq, time, rot} - the structure specifying the data, i.e.,
; flux - the intensity profile, 1D array, in units of sfu/arcsec,
; x - the corresponding coordinates, 1D array, in arcsec,
; freq - the emission frequency, in GHz,
; time - the observation time, string (the same as in SSW map structures),
; rot - the RATAN positional angle, in degrees.
; Note that 1D profiles and 2D maps cannot be mixed together.
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
; a scalar value (applied to all a and b), or
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
; noMultiFreq - if not set (by default), the code optimizes the computations by computing the metrics at all
;  specified frequencies simultaneously. Although the minimization is performed frequency-by-frequency, the Q0 grid
;  and the corresponding metrics computed during the minimization at lower frequencies are then used as pre-computed
;  data at higher frequencies.
;  If set, all frequencies are processed independently; this can be slower, but sometimes more reliable.
;
; DEM, DDM - these keywords are only applicable if the chosen EBTELfileName .sav file contains both the DEM and DDM 
;  tables. 
;  In this case, if the /DEM keyword is set, the code loads the DEM table only (the DDM table is ignored). 
;  Similarly, if the /DDM keyword is set, the code loads the DDM table only (the DEM table is ignored). 
;  If both /DEM and /DDM keywords (or none of them) are set, the code loads both tables.
;  If the chosen file contains only one EBTEL table (either DEM or DDM), the code loads that table; 
;  the /DEM and /DDM keywords are ignored.
;  
; Qstep - the initial relative step over Q0 to search for the optimal heating rate value (must be >1).
; Default: the golden ratio value (1.6180339). 
; 
; xy_shift - shift applied to the observed microwave maps/profiles, a 2-element vector in the form of xy_shift=[dx, dy] 
;  (for 2D maps), or a scalar value (for 1D profiles), in arcseconds.
;  If this parameter is not specified (by default), the shift is computed automatically each time 
;  (i.e., for each frequency and a, b, and Q0 values) to provide the maximum correlation between the observed 
;  and model images/profiles. 
; 
; loud - if set, the code displays more detailed information when it fails to find a solution (e.g., when the minimization
;  procedure goes beyond the EBTEL table).
;
; SHtable - a 7*7 table specifying the selective heating coefficients applied to the field lines with different
;  footpoint combinations. Default: no selective heating (all elements of the table equal 1).
;  
; Nthreads - number of processor threads used for computing the model microwave images. Cannot exceed
;            the number of available processors. Default: a system-defined value (typically, the number 
;            of available processors).
;            
; analyticalNT - defines how the code computes the emission if the parameters of a closed field line (Q, L) 
;                fall beyond the used EBTEL table. If not set (by default): the isothermal barometric formula
;                with the base density of 1e8 cm^{-3} and the temperature of 1 MK is used. If set:
;                the plasma density and temperature are computed using approximate analytical formulae
;                for continuosly heated coronal loops. Note that for the open magnetic lines, the isothermal
;                barometric formula is always used.
;                
; EMthreshold - the maximum allowed EBTEL miss ratio. The EBTEL miss ratio is computed as the ratio of the 
;               number of the closed field lines where the parameter Q falls beyond the used EBTEL table
;               (is too large or too small) to the total number of closed field lines. If the EBTEL miss ratio
;               exceeds the threshold, the corresponding (Q0, a, b) combination is considered falling beyond
;               an acceptable range. Default: 0.1. The threshold is not applicable if there are no closed field 
;               lines.
;
;Results:
; As the result, for each (a, b) combination the program creates in the OutDir directory a .sav file
; with the name starting with 'fit' and including the used metric, threshold, indicator of the multithermal 
; approach, a and b values, and (optionally) 
; the ObsDateTime string.
; These .sav files contain the following fields:
;  freqList - array of the emission frequencies, in GHz.
;  bestQarr - array of the obtained best-fit heating rates Q0 at different frequencies.
;  rhoArr, chiArr, etaArr - arrays of the obtained rho^2, chi^2, and eta^2 metrics at different frequencies.
;                           Note that only one of those metrics (defined by the 'metric' keyword) is actually
;                           minimized; two other metrics correspond to the obtained best-fit Q0 values. 
;  modImageArr - (multi-frequency) map object containing the best-fit model radio maps (corresponding to the
;                best-fit heating rates Q0) at different frequencies. The maps are not convolved with the 
;                instrument beam.
;  modImageConvArr - if the input is in the form of 2D maps, this is a (multi-frequency) map object containing 
;                    the above-mentioned best-fit model radio maps convolved with the instrument beam.
;  obsImageArr - if the input is in the form of 2D maps, this is a (multi-frequency) map object containing the observed 
;                radio maps rebinned and shifted to match the best-fit model maps at the corresponding frequencies.
; If the input is in the form of 1D profiles, the fields modImageConvArr and obsImageArr are lists of structures
; (in the same format as described above) specifying the model 1D scans and the observed 1D scans rebinned and
; shifted to match the models, respectively.
; If the algorithm failed to find the best-fit heating rate at a certain frequency (e.g., the used metric has 
; no minimum within the valid Q0 range, or has more than one local minimum), the corresponding bestQarr, rhoArr, 
; chiArr, and etaArr are set to NaN, and the corresponding image maps contain all zeros.
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
;  ItotalObs, ItotalMod - 3D arrays (N_a*N_b*N_freq) of the total observed and model radio fluxes at different values of
;                         a, b, and frequency. The fluxes correspond to the obtained best-fit Q0 values.
;  CC - 3D array (N_a*N_b*N_freq) of the correlation coefficients of the observed and model radio maps at 
;       different values of a, b, and frequency. The coefficients correspond to the obtained best-fit Q0 values.
;  shiftX, shiftY - 3D arrays (N_a*N_b*N_freq) of the shifts (in arcseconds) applied to the observed radio maps
;                   to obtain the best correlation with the model maps, at different values of a, b, and frequency.
;                   The shifts correspond to the obtained best-fit Q0 values. If the input is in the form of 1D
;                   profiles, shiftY is always zero.
;  rho, chi, eta - 3D arrays (N_a*N_b*N_freq) of the obtained rho^2, chi^2, and eta^2 metrics at different values of 
;                  a, b, and frequency. Note that only one of those metrics (defined by the 'metric' keyword) is 
;                  actually minimized; two other metrics correspond to the obtained best-fit Q0 values. 
; If the Summary*.sav file exists, it will be overwritten.

 if RefDir ne '' then RefDir=file_dirname(RefDir+path_sep()+'*', /mark_directory)
 if OutDir ne '' then OutDir=file_dirname(OutDir+path_sep()+'*', /mark_directory)

 if exist(RefFiles) then ObsFileNames=RefDir+RefFiles else ObsFileNames=file_search(RefDir+'*.sav')
 LoadObservations, ObsFileNames, obsImaps, obsSImaps, obsInfo
 instrument=obsInfo.id
 if obsInfo.id eq 'RATAN' then begin
  sxArr=0
  syArr=0
  beamArr=0
 endif else begin
  sxArr=obsInfo.sx
  syArr=obsInfo.sy
  beamArr=obj_new('map')
  for i=0, obsInfo.Nfreq-1 do begin
   beam=(obsInfo.psf)[i]
   m=make_map(beam, xc=0, yc=0, dx=obsInfo.psf_dx[i], dy=obsInfo.psf_dy[i])
   beamArr->setmap, i, m
  endfor
 endelse

 model=LoadGXmodel(ModelFileName, /noVoxelID)
 
 ebtel=LoadEBTEL(EBTELfileName, DEM=DEM, DDM=DDM)
 DEM_on=ebtel.DEM_on
 DDM_on=ebtel.DDM_on
 
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
 threshold_img=threshold
 
 if ~exist(metric) then metric='eta'
 if (metric ne 'eta') && (metric ne 'chi') && (metric ne 'rho') then metric='eta'
 
 if ~exist(MultiThermal) then iso=1 else iso=(MultiThermal eq 0)
 
 if ~exist(noMultiFreq) then MultiFreq_on=1 else MultiFreq_on=(noMultiFreq eq 0)

 ObsID=exist(ObsDateTime) ? ObsDateTime : ''
 
 if exist(ObsDateTime) then ObsDateTime='_'+ObsDateTime else ObsDateTime=''
 
 if ~exist(Qstep) then Qstep=(1d0+sqrt(5d0))/2
 
 if exist(xy_shift) then fixed_shifts=1 else begin
  xy_shift=0
  fixed_shifts=0
 endelse

 if exist(SHtable) then begin
  SHtable=double(SHtable)
  s=size(SHtable)
  if (s[0] ne 2) || (s[1] ne 7) || (s[2] ne 7) then begin
   print, 'Incorrect size of the selective heating table; returning to default.'
   SHtable=dblarr(7, 7)
   SHtable[*]=1
  endif
 endif
 
 if ~exist(EMthreshold) then EMthreshold=0.1

 simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, obsInfo.freq, $
                          rot=(obsInfo.id eq 'RATAN') ? obsInfo.rot[0] : 0d0, Nthreads=Nthreads)       
 
 tstart0=systime(1)
 
 for i=0, N_a-1 do for j=0, N_b-1 do begin
  a=alist[i]
  b=blist[j]
  Qstart=Q0start[i, j]
  
  print, 'Computing the best fit Q for a=', a, ', b=', b
  
  fname=OutDir+'fit_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
        (iso ? '_I' : '_M')+ObsDateTime+$
        '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav' 
  
  if ~file_exist(fname) then begin
   tstart1=systime(1)
   
   FindBestFitQ, LibFileName, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ 
                 a, b, Qstart, Qstep, iso, threshold_img, metric, MultiFreq_on, fixed_shifts, xy_shift, $         
                 freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $        
                 ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                 obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                 modFlagArr, allQ, allMetrics, loud=loud, SHtable=SHtable, Nthreads=Nthreads, $
                 analyticalNT=analyticalNT, EMthreshold=EMthreshold
         
   save, LibFileName, modelFileName, EBTELfileName, DEM_on, DDM_on, $
         sxArr, syArr, beamArr, $
         a, b, Qstart, Qstep, iso, threshold_img, metric, MultiFreq_on, fixed_shifts, $
         freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ 
         ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
         obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
         modFlagArr, allQ, allMetrics, instrument, $
         filename=fname, /compress
   
   print, 'The best fit Q for a=', a, ', b=', b, ' computed in ', systime(1)-tstart1, ' s'
   print, 'Total elapsed time: ', systime(1)-tstart0, ' s'
  endif
 endfor
 
 print, 'Creating the summary file'
 
 freqList=obsInfo.freq
 N_freq=n_elements(freqList)
 
 bestQ=dblarr(N_a, N_b, N_freq)
 chi=dblarr(N_a, N_b, N_freq)
 rho=dblarr(N_a, N_b, N_freq)
 eta=dblarr(N_a, N_b, N_freq)
 ItotalObs=dblarr(N_a, N_b, N_freq)
 ItotalMod=dblarr(N_a, N_b, N_freq)
 CC=dblarr(N_a, N_b, N_freq)    
 shiftX=dblarr(N_a, N_b, N_freq)
 shiftY=dblarr(N_a, N_b, N_freq)
 
 for i=0, N_a-1 do for j=0, N_b-1 do begin
  a=alist[i]
  b=blist[j]
  
  fname=OutDir+'fit_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
        (iso ? '_I': '_M')+ObsDateTime+$
        '_a'+string(a, format='(F+6.3)')+'_b'+string(b, format='(F+6.3)')+'.sav' 
  o=obj_new('IDL_Savefile', fname)
  o->restore, 'chiArr'
  o->restore, 'rhoArr'
  o->restore, 'etaArr'
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
   
   if obsInfo.id eq 'RATAN' then begin
    m=obsImageArr[k]
    shiftX[i, j, k]=m.shiftX
   endif else begin
    m=obsImageArr.getmap(k)
    if tag_exist(m, 'shiftX') then shiftX[i, j, k]=m.shiftX
    if tag_exist(m, 'shiftY') then shiftY[i, j, k]=m.shiftY
   endelse 
  endfor
 endfor 
 
 fname1=OutDir+'Summary_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
        (iso ? '_I': '_M')+ObsDateTime+'.sav'
 
 save, alist, blist, freqList, bestQ, ItotalObs, ItotalMod, CC, chi, rho, eta, metric, iso, threshold_img, $
       modelFileName, EBTELfileName, DEM_on, DDM_on, MultiFreq_on, ObsID, shiftX, shiftY, fixed_shifts, instrument, $
       filename=fname1
 
 print, 'Done'
end