pro FindBestFitQ, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, a, b, Qstart, iso, $
                  bestQarr, chiArr, chiVarArr, IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                  freqList, allQ, allChi, modImageConvArr, obsImageArr, thr=thr, metric=metric
 forward_function GetSmoothedMax, DefineCoronaParms, ReserveOutputSpace, ConvertToMaps
                       
 if ~exist(thr) then thr=0.1d0 ;default map threshold
 if ~exist(metric) then metric='eta' ;default metric to use
 
 acc=1d-2 ;desired accuracy
 Tbase=1d6 ;analytical corona temperature
 nbase=1d8 ;analytical corona base density
 
 Nfreq=obsInfo.Nfreq
 freqList=obsInfo.freq
 
 Qgrid=[Qstart]
 fdone=[0]
 modImages=objarr(1)
 modImagesConv=objarr(1)
 obsImages=objarr(1)
 flags=lonarr(1, 6)
 chi=dblarr(1, Nfreq)
 chiVar=dblarr(1, Nfreq)
 Iobs=dblarr(1, Nfreq)
 Imod=dblarr(1, Nfreq)
 CC=dblarr(1, Nfreq)
 
 bestQarr=dblarr(Nfreq)
 chiArr=dblarr(Nfreq)
 chiVarArr=dblarr(Nfreq)
 IobsArr=dblarr(Nfreq)
 ImodArr=dblarr(Nfreq)
 CCarr=dblarr(Nfreq) 
 modImageArr=obj_new('map')
 modImageConvArr=obj_new('map')
 obsImageArr=obj_new('map')
 modFlagArr=lonarr(Nfreq, 6)
 
 bestQarr[*]=!values.d_NaN
 chiArr[*]=!values.d_NaN
 chiVarArr[*]=!values.d_NaN
 IobsArr[*]=!values.d_NaN
 ImodArr[*]=!values.d_NaN
 CCarr[*]=!values.d_NaN 
 
 G=(1d0+sqrt(5d0))/2
 
 ;---------------------------------------------------------
 
 done=0
 
 while ~done do begin
  NQ=n_elements(Qgrid)
  
  for i=0, NQ-1 do if ~fdone[i] then begin
   print, 'a=', string(a, format='(F0)'), ' b=', string(b, format='(F0)')
   print, 'Array of Q0: ', Qgrid
   print, 'Computing images for the item #', string(i, format='(I0)')  
   
   coronaparms=DefineCoronaParams(Tbase, nbase, Qgrid[i], a, b, force_isothermal=iso)
   outspace=ReserveOutputSpace(simbox)
    
   r=call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace)
                    
   ConvertToMaps, outspace, simbox, model, modImaps, modVmaps
   obj_destroy, modVmaps
   modImages[i]=modImaps
   flags[i, *]=outspace.flagsCorona
   
   modX=obj_new('map')
   obsX=obj_new('map')
   
   for j=0, Nfreq-1 do begin
    modI=modImaps.getmap(j)
    _obsI=obsImaps.getmap(j)
    _obsSigma=obsSImaps.getmap(j)
    
    MakeLocalBeam, obsInfo, j, modI.dx, modI.dy, beam
    FixLocalBeam, beam, simbox.Nx, simbox.Ny
    modI.data=convol_fft(modI.data, beam)
    modX->setmap, j, modI
    
    FindShift, _obsI, modI, dx, dy
    ExtractSubmap, _obsI, modI, dx, dy, 'SRH I', obsI
    ExtractSubmap, _obsSigma, modI, dx, dy, 'SRH I sigma', obsSigma
    obsX->setmap, j, obsI
    
    obsMax=GetSmoothedMax(obsI, obsInfo.sx[j], obsInfo.sy[j])
    modMax=max(modI.data)    
    
    u=where(obsI.data gt (obsMax*thr), ms)
    maskObs=1d0*ms/n_elements(obsI.data)
    
    u=where(modI.data gt (modMax*thr), ms)
    maskMod=1d0*ms/n_elements(modI.data)    
    
    u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)), ms)
    Iobs[i, j]=total(obsI.data[u])*obsI.dx*obsI.dy
    Imod[i, j]=total(modI.data[u])*modI.dx*modI.dy
    CC[i, j]=c_correlate(obsI.data[u], modI.data[u], 0)
    
    case metric of
     'chi': begin
             chi[i, j]=mean(((modI.data[u]-obsI.data[u])/obsSigma.data[u])^2)    ;\chi^2
             chiVar[i, j]=variance((modI.data[u]-obsI.data[u])/obsSigma.data[u]) ;\chi^2 corrected
            end
     'eta': begin     
             chi[i, j]=mean(((modI.data[u]-obsI.data[u])/mean(obsI.data[u]))^2)    ;\eta^2
             chiVar[i, j]=variance((modI.data[u]-obsI.data[u])/mean(obsI.data[u])) ;\eta^2 corrected
            end
     'rho': begin
             chi[i, j]=mean((modI.data[u]/obsI.data[u]-1)^2)    ;\rho^2
             chiVar[i, j]=variance(modI.data[u]/obsI.data[u]-1) ;\rho^2 corrected
            end
    endcase            
    
    if (maskMod gt 0.99) || ((maskMod/maskObs) gt 4) then begin
     chi[i, j]=!values.d_NaN
     chiVar[i, j]=!values.d_NaN
    endif
   endfor
   
   modImagesConv[i]=modX
   obsImages[i]=obsX
 
   fdone[i]=1
  endif
  
  if NQ eq 1 then begin
   aw=total(sgn(Iobs-Imod))
   if aw eq 0 then aw=1 
  endif else begin
   lmins=0
   rmins=0
   
   for j=0, Nfreq-1 do begin
    u=where(~finite(chi[*, j]), k)
    if k eq 0 then begin
     chi_min=min(chi[*, j], k)
     if k eq 0 then lmins+=1
     if k eq (NQ-1) then rmins+=1
    endif
   endfor
   
   if (1d0*flags[0, 5]/flags[0, 2]) gt 0.1 then lmins=0
   if (1d0*flags[NQ-1, 5]/flags[NQ-1, 2]) gt 0.1 then rmins=0
   
   if (lmins eq 0) && (rmins eq 0) then aw=0 else begin
    aw=rmins-lmins
    if aw eq 0 then aw=1
   endelse
  endelse
  
  if aw lt 0 then begin
   Qgrid=[min(Qgrid)/G, Qgrid]
   fdone=[0, fdone]
   modImages=[obj_new(), modImages]
   modImagesConv=[obj_new(), modImagesConv]
   obsImages=[obj_new(), obsImages]
      
   flags2=lonarr(NQ+1, 6)
   flags2[1 : NQ, *]=flags
   flags=flags2
 
   chi2=dblarr(NQ+1, Nfreq)
   chi2[1 : NQ, *]=chi
   chi=chi2
   
   chiVar2=dblarr(NQ+1, Nfreq)
   chiVar2[1 : NQ, *]=chiVar
   chiVar=chiVar2   
      
   Iobs2=dblarr(NQ+1, Nfreq)
   Iobs2[1 : NQ, *]=Iobs
   Iobs=Iobs2
   
   Imod2=dblarr(NQ+1, Nfreq)
   Imod2[1 : NQ, *]=Imod
   Imod=Imod2
   
   CC2=dblarr(NQ+1, Nfreq)
   CC2[1 : NQ, *]=CC
   CC=CC2      
  endif else if aw gt 0 then begin
   Qgrid=[Qgrid, max(Qgrid)*G]
   fdone=[fdone, 0]
   modImages=[modImages, obj_new()]
   modImagesConv=[modImagesConv, obj_new()]
   obsImages=[obsImages, obj_new()]
   
   flags2=lonarr(NQ+1, 6)
   flags2[0 : NQ-1, *]=flags
   flags=flags2
   
   chi2=dblarr(NQ+1, Nfreq)
   chi2[0 : NQ-1, *]=chi
   chi=chi2
   
   chiVar2=dblarr(NQ+1, Nfreq)
   chiVar2[0 : NQ-1, *]=chiVar
   chiVar=chiVar2   
   
   Iobs2=dblarr(NQ+1, Nfreq)
   Iobs2[0 : NQ-1, *]=Iobs
   Iobs=Iobs2   
   
   Imod2=dblarr(NQ+1, Nfreq)
   Imod2[0 : NQ-1, *]=Imod
   Imod=Imod2      
   
   CC2=dblarr(NQ+1, Nfreq)
   CC2[0 : NQ-1, *]=CC
   CC=CC2        
  endif else done=1
 endwhile
 
 ;-----------------------------------------------------------------
 
 badf=intarr(Nfreq)
 
 for j=0, Nfreq-1 do begin
  u=where(~finite(chi[*, j]), k)
  
  if k gt 0 then badf[j]=1 else begin
   chi_min=min(chi[*, j], k)
   
   if (k eq 0) || (k eq (n_elements(Qgrid)-1)) then badf[j]=1 else begin
    nmin=0
    for i=1, n_elements(Qgrid)-2 do if (chi[i, j] lt chi[i-1, j]) && (chi[i, j] lt chi[i+1, j]) then nmin+=1
    if nmin ne 1 then badf[j]=1
   endelse 
  endelse
 endfor
 
 u=where(~badf, k)
 if k gt 0 then print, 'Good frequencies: ', freqList[u]
 
 u=where(badf, k)
 if k gt 0 then print, 'Bad frequencies: ', freqList[u]
 
 for j=0, Nfreq-1 do if ~badf[j] then begin
  done=0
  step=1
  
  while ~done do begin
   NQ=n_elements(Qgrid) 
    
   chi_b=min(chi[*, j], ib)
   chi_a=chi[ib-1, j]
   chi_c=chi[ib+1, j]
   Qa=Qgrid[ib-1]
   Qb=Qgrid[ib]
   Qc=Qgrid[ib+1]

   forward_function IRatio
   
   if ((Qc-Qa)/(Qc+Qa)) lt acc then done=1 else begin
    QxG=((Qc-Qb) gt (Qb-Qa)) ? Qb+(Qc-Qb)*(1d0-1d0/G) : $
                               Qb-(Qb-Qa)*(1d0-1d0/G)
    QxB=Qb-0.5*((Qb-Qa)^2*(chi_b-chi_c)-(Qb-Qc)^2*(chi_b-chi_a))/$
               ((Qb-Qa)*(chi_b-chi_c)-(Qb-Qc)*(chi_b-chi_a))   
             
    dNewG=(QxG gt Qb) ? (Qc-Qb+QxG-Qa)/2 : (Qb-Qa+Qc-QxG)/2
    dNewB=(QxB gt Qb) ? (Qc-Qb+QxB-Qa)/2 : (Qb-Qa+Qc-QxB)/2

    if QxB gt Qb then begin
     rB_hit=IRatio(Qc-QxB, QxB-Qb)
     rB_miss=IRatio(QxB-Qb, Qb-Qa)
    endif else begin
     rB_hit=IRatio(Qb-QxB, QxB-Qa)
     rb_miss=IRatio(Qc-Qb, Qb-QxB)
    endelse
             
    if (QxB ge Qc) || (QxB le Qa) || $
       (dNewB gt dNewG) || (rB_hit gt 10) || (rB_miss gt 10) then begin
     Qx=QxG
     st=' (golden step)'
    endif else begin
     Qx=QxB
     st=' (Brent step)'
    endelse           
                              
    print, 'a=', string(a, format='(F0)'), ' b=', string(b, format='(F0)'), $
           ' ', string(freqlist[j], format='(F0)'), ' GHz, minimization step #', string(step, format='(I0)')
    print, 'Bracketed range: ', Qa, Qb, Qc, ' (', (Qc-Qa)/(Qc+Qa)*100, '%)'
    print, 'Image metrics:   ', chi_a, chi_b, chi_c
    print, 'Computing images for Q0=', Qx, st
   
    coronaparms=DefineCoronaParams(Tbase, nbase, Qx, a, b, force_isothermal=iso)
    outspace=ReserveOutputSpace(simbox)
    
    r=call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace)
                    
    ConvertToMaps, outspace, simbox, model, modImaps, modVmaps
    obj_destroy, modVmaps
    
    chi_x=dblarr(Nfreq)
    chiVar_x=dblarr(Nfreq)
    Iobs_x=dblarr(Nfreq)
    Imod_x=dblarr(Nfreq)
    CC_x=dblarr(Nfreq)
    
    modX=obj_new('map')
    obsX=obj_new('map')
    
    for k=0, Nfreq-1 do begin
     modI=modImaps.getmap(k)
     _obsI=obsImaps.getmap(k)
     _obsSigma=obsSImaps.getmap(k)
    
     MakeLocalBeam, obsInfo, k, modI.dx, modI.dy, beam
     FixLocalBeam, beam, simbox.Nx, simbox.Ny
     modI.data=convol_fft(modI.data, beam)
     modX->setmap, k, modI
    
     FindShift, _obsI, modI, dx, dy
     ExtractSubmap, _obsI, modI, dx, dy, 'SRH I', obsI
     ExtractSubmap, _obsSigma, modI, dx, dy, 'SRH I sigma', obsSigma
     obsX->setmap, k, obsI
    
     obsMax=GetSmoothedMax(obsI, obsInfo.sx[k], obsInfo.sy[k])
     modMax=max(modI.data)    
    
     u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)), ms)
     Iobs_x[k]=total(obsI.data[u])*obsI.dx*obsI.dy
     Imod_x[k]=total(modI.data[u])*modI.dx*modI.dy
     CC_x[k]=c_correlate(obsI.data[u], modI.data[u], 0)
          
     case metric of
      'chi': begin
              chi_x[k]=mean(((modI.data[u]-obsI.data[u])/obsSigma.data[u])^2)    ;\chi^2
              chiVar_x[k]=variance((modI.data[u]-obsI.data[u])/obsSigma.data[u]) ;\chi^2 corrected
             end
      'eta': begin     
              chi_x[k]=mean(((modI.data[u]-obsI.data[u])/mean(obsI.data[u]))^2)    ;\eta^2
              chiVar_x[k]=variance((modI.data[u]-obsI.data[u])/mean(obsI.data[u])) ;\eta^2 corrected
             end
      'rho': begin
              chi_x[k]=mean((modI.data[u]/obsI.data[u]-1)^2)    ;\rho^2
              chiVar_x[k]=variance(modI.data[u]/obsI.data[u]-1) ;\rho^2 corrected
             end
     endcase  
    endfor 
    
    l=(Qx gt Qb) ? ib : ib-1
    
    Qgrid=[Qgrid[0 : l], Qx, Qgrid[l+1 : NQ-1]]
    modImages=[modImages[0 : l], modImaps, modImages[l+1 : NQ-1]]
    modImagesConv=[modImagesConv[0 : l], modX, modImagesConv[l+1 : NQ-1]]
    obsImages=[obsImages[0 : l], obsX, obsImages[l+1 : NQ-1]]
    
    flags2=lonarr(NQ+1, 6)
    flags2[0 : l, *]=flags[0 : l, *]
    flags2[l+1, *]=outspace.flagsCorona
    flags2[l+2 : NQ, *]=flags[l+1 : NQ-1, *]
    flags=flags2
    
    chi2=dblarr(NQ+1, Nfreq)
    chi2[0 : l, *]=chi[0 : l, *]
    chi2[l+1, *]=chi_x
    chi2[l+2 : NQ, *]=chi[l+1 : NQ-1, *]
    chi=chi2
    
    chiVar2=dblarr(NQ+1, Nfreq)
    chiVar2[0 : l, *]=chiVar[0 : l, *]
    chiVar2[l+1, *]=chiVar_x
    chiVar2[l+2 : NQ, *]=chiVar[l+1 : NQ-1, *]
    chiVar=chiVar2    
    
    Iobs2=dblarr(NQ+1, Nfreq)
    Iobs2[0 : l, *]=Iobs[0 : l, *]
    Iobs2[l+1, *]=Iobs_x
    Iobs2[l+2 : NQ, *]=Iobs[l+1 : NQ-1, *]
    Iobs=Iobs2
    
    Imod2=dblarr(NQ+1, Nfreq)
    Imod2[0 : l, *]=Imod[0 : l, *]
    Imod2[l+1, *]=Imod_x
    Imod2[l+2 : NQ, *]=Imod[l+1 : NQ-1, *]
    Imod=Imod2
    
    CC2=dblarr(NQ+1, Nfreq)
    CC2[0 : l, *]=CC[0 : l, *]
    CC2[l+1, *]=CC_x
    CC2[l+2 : NQ, *]=CC[l+1 : NQ-1, *]
    CC=CC2                    
   endelse
   
   step+=1
  endwhile
  
  chiArr[j]=min(chi[*, j], ib)
  chiVarArr[j]=chiVar[ib, j]
  IobsArr[j]=Iobs[ib, j]
  ImodArr[j]=Imod[ib, j]
  CCarr[j]=CC[ib, j]  
  bestQarr[j]=Qgrid[ib]
  modFlagArr[j, *]=flags[ib, *]
  
  bmap=modImages[ib]
  m=bmap.getmap(j)
  modImageArr->setmap, j, m
  bmap=modImagesConv[ib]
  m=bmap.getmap(j)
  modImageConvArr->setmap, j, m
  bmap=obsImages[ib]
  m=bmap.getmap(j)
  obsImageArr->setmap, j, m  
 endif else begin
  bmap=modImages[0]
  m=bmap.getmap(j)
  m.data[*]=0
  modImageArr->setmap, j, m
  modImageConvArr->setmap, j, m
  
  _obsI=obsImaps.getmap(j)
  ExtractSubmap, _obsI, m, 0, 0, 'SRH I', obsI
  obsImageArr->setmap, j, obsI
 endelse
 
 allQ=Qgrid
 allChi=chi
end