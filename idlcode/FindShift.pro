function correlateMaps, map1, map2
 m1=mean(map1.data)
 m2=mean(map2.data)
 return, total((map1.data-m1)*(map2.data-m2))/sqrt(total((map1.data-m1)^2))/sqrt(total((map2.data-m2)^2))
end

pro FindShift, obsmap, modmap, dx, dy
 n=301
 d_dx=1d0
 d_dy=1d0
 
 iarr=indgen(n)-(n-1)/2
 jarr=indgen(n)-(n-1)/2
 
 sarr=dblarr(n, n)
 sarr[*]=!values.d_nan
 
 i=(n-1)/2
 j=(n-1)/2
 
 ExtractSubmap, obsmap, modmap, d_dx*iarr[i], d_dy*jarr[j], '*', newmap
 sarr[i, j]=correlateMaps(modmap, newmap)
 
 done=0
 
 while ~done do begin
  for ii=-1, 1 do for jj=-1, 1 do if ~finite(sarr[i+ii, j+jj]) then begin
   ExtractSubmap, obsmap, modmap, d_dx*iarr[i+ii], d_dy*jarr[j+jj], '*', newmap
   sarr[i+ii, j+jj]=correlateMaps(modmap, newmap)
  endif
  
  smax=-1d100
  for ii=-1, 1 do for jj=-1, 1 do if sarr[i+ii, j+jj] gt smax then begin
   smax=sarr[i+ii, j+jj]
   ii0=ii
   jj0=jj
  endif
  
  if (ii0 eq 0) && (jj0 eq 0) then done=1 else begin
   i=i+ii0
   j=j+jj0
  endelse
 endwhile
 
 dx=d_dx*iarr[i]
 dy=d_dy*jarr[j]
end       
