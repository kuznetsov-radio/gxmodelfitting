pro FindBestFitQ, model, obsImaps, obsInfo, a, b, Qstart, thr, bestQarr, Qlist, modImageList, $
                  reprocess=reprocess, badf=badf
 forward_function GetSmoothedMax, interpol_localized

 Nfreq=obsInfo.Nfreq
 freqlist=obsInfo.freq

 if ~keyword_set(badf) then badf=intarr(Nfreq)
 
 if keyword_set(reprocess) then begin
  fdone=intarr(n_elements(Qlist))
  Iobs=dblarr(n_elements(Qlist), Nfreq)
  Imod=dblarr(n_elements(Qlist), Nfreq)  
 endif else begin
  Qlist=[Qstart]
  fdone=[0]
  modImageList=objarr(1)
  Iobs=dblarr(1, Nfreq)
  Imod=dblarr(1, Nfreq)
 endelse 
 
 bestQarr=dblarr(Nfreq)
 
 alldone=0
 
 while ~alldone do begin
  for i=0, n_elements(Qlist)-1 do if ~fdone[i] then begin
   if modImageList[i] then modImaps=modImageList[i] else begin
    ComputeGXmaps, model, freqlist, Qlist[i], a, b, modImaps, modVmaps
    obj_destroy, modVmaps
   endelse 
   
   for j=0, Nfreq-1 do if ~badf[j] then begin
    modI=modImaps.getmap(j)
    _obsI=obsImaps.getmap(j)
    
    MakeLocalBeam, obsInfo, j, modI.dx, modI.dy, beam
    modI.data=convol(modI.data, beam, /edge_zero)
    
    FindShift, _obsI, modI, dx, dy
    ExtractSubmap, _obsI, modI, dx, dy, 'SRH I', obsI
    
    obsMax=GetSmoothedMax(obsI, obsInfo.sx[j], obsInfo.sy[j])
    modMax=max(modI.data)
    
    u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)))
      
    Iobs[i, j]=total(obsI.data[u])*obsI.dx*obsI.dy
    Imod[i, j]=total(modI.data[u])*modI.dx*modI.dy
   endif
   
   modImageList[i]=modImaps
   fdone[i]=1             
  endif
  
  NQ=n_elements(Qlist)
  if NQ eq 1 then begin
   aw=total(sgn(Iobs-Imod))
   if aw eq 0 then aw=1
  endif else begin
   Qm=exp(alog(min(Qlist/1.001))+(alog(max(Qlist*1.001))-alog(min(Qlist/1.001)))*dindgen(10000)/9999)
   
   lmin=10000
   rmin=10000
   
   for j=0, Nfreq-1 do if ~badf[j] then begin
    Imod2=exp(interpol_localized(alog(Imod[*, j]), alog(Qlist), alog(Qm)))
    Iobs2=exp(interpol_localized(alog(Iobs[*, j]), alog(Qlist), alog(Qm)))
 
    res=min(abs(Iobs2-Imod2), k)
    bestQarr[j]=Qm[k]
    
    u=where(Qlist lt bestQarr[j], kl)
    u=where(Qlist gt bestQarr[j], kr)
    lmin=lmin<kl
    rmin=rmin<kr
   endif
 
   if (lmin ge 2) && (rmin ge 2) then aw=0 $
   else aw=(lmin lt rmin) ? -1 : 1
  endelse
  
  if aw lt 0 then begin
   Qlist=[min(Qlist)/2, Qlist]
   fdone=[0, fdone]
   modImageList=[obj_new(), modImageList]
   Iobs2=dblarr(NQ+1, Nfreq)
   Iobs2[1 : NQ, *]=Iobs
   Iobs=Iobs2
   Imod2=dblarr(NQ+1, Nfreq)
   Imod2[1 : NQ, *]=Imod
   Imod=Imod2
  endif else if aw gt 0 then begin
   Qlist=[Qlist, max(Qlist)*2]
   fdone=[fdone, 0]
   modImageList=[modImageList, obj_new()]
   Iobs2=dblarr(NQ+1, Nfreq)
   Iobs2[0 : NQ-1, *]=Iobs
   Iobs=Iobs2
   Imod2=dblarr(NQ+1, Nfreq)
   Imod2[0 : NQ-1, *]=Imod
   Imod=Imod2
  endif else alldone=1
 endwhile
end