pro FindShift1D, obsscan, modscan, dx
 d_dx=1d0
 
 dx=0d0
 done=0
 
 while ~done do begin
  cc_l=c_correlate(modscan.flux, interpol(obsscan.flux, obsscan.x, modscan.x+dx-d_dx, /spline), 0, /covariance)
  cc_c=c_correlate(modscan.flux, interpol(obsscan.flux, obsscan.x, modscan.x+dx, /spline), 0, /covariance)
  cc_r=c_correlate(modscan.flux, interpol(obsscan.flux, obsscan.x, modscan.x+dx+d_dx, /spline), 0, /covariance)
  
  a=max([cc_l, cc_c, cc_r], k)
  if k eq 1 then done=1 else if k eq 0 then dx-=d_dx else dx+=d_dx
  
  if min(modscan.x+dx) lt min(obsscan.x) then begin
   dx+=d_dx
   done=1
  endif
  if max(modscan.x+dx) gt max(obsscan.x) then begin
   dx-=d_dx
   done=1
  endif
 endwhile
end