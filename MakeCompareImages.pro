pro MakeCompareImages, model, obsImaps, obsSImaps, obsInfo, a, b, bestQarr, thr, modImagesI, modImagesV, $
                       rho, chi, avg_subrange=avg_subrange
 forward_function GetSmoothedMax 

 Nfreq=obsInfo.Nfreq
 freqlist=obsInfo.freq
 
 rho=dblarr(Nfreq)
 chi=dblarr(Nfreq)
 
 modImagesI=obj_new('map')
 modimagesV=obj_new('map')
 
 if ~keyword_set(avg_subrange) then for i=0, Nfreq-1 do begin
  ComputeGXmaps, model, freqlist[i], bestQarr[i], a, b, modImaps, modVmaps
  
  modImagesI->setmap, i, modImaps.getmap(0)
  modImagesV->setmap, i, modVmaps.getmap(0)
  
  obj_destroy, modImaps
  obj_destroy, modVmaps
 endfor else begin
  ComputeGXmaps, model, freqlist, mean(bestQarr[avg_subrange]), a, b, modImaps, modVmaps
  
  for i=0, Nfreq-1 do begin
   modImagesI->setmap, i, modImaps.getmap(i)
   modImagesV->setmap, i, modVmaps.getmap(i)
  endfor
  
  obj_destroy, modImaps
  obj_destroy, modVmaps
 endelse
 
 for i=0, Nfreq-1 do begin
  modI=modImagesI.getmap(i)
  _obsI=obsImaps.getmap(i)
  _obsSigma=obsSImaps.getmap(i)
    
  MakeLocalBeam, obsInfo, i, modI.dx, modI.dy, beam
  modI.data=convol(modI.data, beam, /edge_zero)
    
  FindShift, _obsI, modI, dx, dy
  ExtractSubmap, _obsI, modI, dx, dy, 'SRH I', obsI
  ExtractSubmap, _obsSigma, modI, dx, dy, 'SRH I sigma', obsSigma
    
  obsMax=GetSmoothedMax(obsI, obsInfo.sx[i], obsInfo.sy[i])
  modMax=max(modI.data)
    
  u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)))
      
  rho[i]=variance(modI.data[u]/obsI.data[u]-1)
  chi[i]=variance((modI.data[u]-obsI.data[u])/obsSigma.data[u])         
 endfor
end