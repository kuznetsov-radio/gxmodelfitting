function CheckSummary, OutDir, a_min, b_min, Q0, $
                       metric=metric, threshold_img=threshold_img, MultiThermal=MultiThermal, $
                       ObsDateTime=ObsDateTime, ObsFreq=ObsFreq
 if ~exist(threshold_img) then threshold_img=0.1d0

 if ~exist(metric) then metric='eta'
 if (metric ne 'eta') && (metric ne 'chi') && (metric ne 'rho') then metric='eta'
 
 if ~exist(MultiThermal) then iso=1 else iso=MultiThermal eq 0

 if exist(ObsDateTime) then ObsDateTime1='_'+ObsDateTime else ObsDateTime1=''
 
 if exist(ObsFreq) then ObsFreq1='_'+ObsFreq else ObsFreq1=''

 fname=OutDir+'Summary_'+metric+'_thr'+string(threshold_img, format='(F5.3)')+$
       (iso ? '_I': '_M')+ObsDateTime1+ObsFreq1+'.sav'
 print, 'Searching for the file: ', fname      

 if file_exist(fname) then begin
  restore, fname
  
  case metric of
   'chi': mtr=chi
   'rho': mtr=rho
   'eta': mtr=eta
  endcase
  
  u=where(mtr lt 0, g)
  if g gt 0 then mtr[u]=!values.d_NaN
  
  mmin=1d10
  for i=0, n_elements(alist)-1 do for j=0, n_elements(blist)-1 do if finite(mtr[i, j]) && (mtr[i, j] lt mmin) then begin
   mmin=mtr[i, j]
   a_min=alist[i]
   b_min=blist[j]
   Q0=bestQ[i, j]
  endif
  
  print, '... found'
  print, 'Minimum metric: ', mmin, ' at a=', a_min, ', b=', b_min, ' Q0=', Q0
  return, 1
 endif else begin
  print, '... not found'
  return, 0
 endelse
end