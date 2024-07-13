pro FindBestFitQmf, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, a, b, Qstart, iso, $
                    bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                    IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                    freqList, allQ, allMetrics, modImageConvArr, obsImageArr, thr=thr, metric=metric, $
                    Qstep=Qstep, loud=loud
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
 rho=dblarr(1, Nfreq)
 rhoVar=dblarr(1, Nfreq)
 eta=dblarr(1, Nfreq)
 etaVar=dblarr(1, Nfreq)  
 Iobs=dblarr(1, Nfreq)
 Imod=dblarr(1, Nfreq)
 CC=dblarr(1, Nfreq)
 
 bestQarr=dblarr(Nfreq)
 chiArr=dblarr(Nfreq)
 chiVarArr=dblarr(Nfreq)
 rhoArr=dblarr(Nfreq)
 rhoVarArr=dblarr(Nfreq) 
 etaArr=dblarr(Nfreq)
 etaVarArr=dblarr(Nfreq)  
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
 rhoArr[*]=!values.d_NaN
 rhoVarArr[*]=!values.d_NaN
 etaArr[*]=!values.d_NaN
 etaVarArr[*]=!values.d_NaN  
 IobsArr[*]=!values.d_NaN
 ImodArr[*]=!values.d_NaN
 CCarr[*]=!values.d_NaN 
 
 G=(1d0+sqrt(5d0))/2
 if ~exist(Qstep) then Qstep=G
 
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
    ExtractSubmap, _obsI, modI, dx, dy, _obsI.id, obsI
    obsI=create_struct('shiftX', dx, $
                       'shiftY', dy, $
                       obsI)
    ExtractSubmap, _obsSigma, modI, dx, dy, _obsI.id+' sigma', obsSigma
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
    
    chi[i, j]=mean(((modI.data[u]-obsI.data[u])/obsSigma.data[u])^2)    ;\chi^2
    chiVar[i, j]=variance((modI.data[u]-obsI.data[u])/obsSigma.data[u]) ;\chi^2 corrected
    rho[i, j]=mean((modI.data[u]/obsI.data[u]-1)^2)    ;\rho^2
    rhoVar[i, j]=variance(modI.data[u]/obsI.data[u]-1) ;\rho^2 corrected    
    eta[i, j]=mean(((modI.data[u]-obsI.data[u])/mean(obsI.data[u]))^2)    ;\eta^2
    etaVar[i, j]=variance((modI.data[u]-obsI.data[u])/mean(obsI.data[u])) ;\eta^2 corrected
    
    if (maskMod gt 0.99) || ((maskMod/maskObs) gt 4) then begin
     if exist(loud) then print, '*** GR contribution is too low at ', freqList[j], ' GHz ***'
     chi[i, j]=!values.d_NaN
     chiVar[i, j]=!values.d_NaN
     rho[i, j]=!values.d_NaN
     rhoVar[i, j]=!values.d_NaN
     eta[i, j]=!values.d_NaN
     etaVar[i, j]=!values.d_NaN          
    endif
   endfor
   
   case metric of
    'chi': mtr=chi
    'rho': mtr=rho    
    'eta': mtr=eta   
   endcase  
   
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
    u=where(~finite(mtr[*, j]), k)
    if k eq 0 then begin
     mtr_min=min(mtr[*, j], k)
     if k eq 0 then lmins+=1
     if k eq (NQ-1) then rmins+=1
    endif
   endfor
   
   if (1d0*flags[0, 5]/flags[0, 2]) gt 0.1 then begin
    if exist(loud) then print, '*** Out of EBTEL table at left boundary, Q0=', Qgrid[0], ' ***'
    lmins=0
   endif
   if (1d0*flags[NQ-1, 5]/flags[NQ-1, 2]) gt 0.1 then begin
    if exist(loud) then print, '*** Out of EBTEL table at right boundary, Q0=', Qgrid[NQ-1], ' ***'
    rmins=0
   endif
   
   if (lmins eq 0) && (rmins eq 0) then aw=0 else begin
    aw=rmins-lmins
    if aw eq 0 then aw=1
   endelse
  endelse
  
  if aw lt 0 then begin
   Qgrid=[min(Qgrid)/Qstep, Qgrid]
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
   
   rho2=dblarr(NQ+1, Nfreq)
   rho2[1 : NQ, *]=rho
   rho=rho2
   
   rhoVar2=dblarr(NQ+1, Nfreq)
   rhoVar2[1 : NQ, *]=rhoVar
   rhoVar=rhoVar2
   
   eta2=dblarr(NQ+1, Nfreq)
   eta2[1 : NQ, *]=eta
   eta=eta2
   
   etaVar2=dblarr(NQ+1, Nfreq)
   etaVar2[1 : NQ, *]=etaVar
   etaVar=etaVar2            
      
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
   Qgrid=[Qgrid, max(Qgrid)*Qstep]
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
   
   rho2=dblarr(NQ+1, Nfreq)
   rho2[0 : NQ-1, *]=rho
   rho=rho2
   
   rhoVar2=dblarr(NQ+1, Nfreq)
   rhoVar2[0 : NQ-1, *]=rhoVar
   rhoVar=rhoVar2      
   
   eta2=dblarr(NQ+1, Nfreq)
   eta2[0 : NQ-1, *]=eta
   eta=eta2
   
   etaVar2=dblarr(NQ+1, Nfreq)
   etaVar2[0 : NQ-1, *]=etaVar
   etaVar=etaVar2      
   
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
  u=where(~finite(mtr[*, j]), k)
  
  if k gt 0 then badf[j]=1 else begin
   mtr_min=min(mtr[*, j], k)
   
   if (k eq 0) || (k eq (n_elements(Qgrid)-1)) then badf[j]=1 else begin
    nmin=0
    for i=1, n_elements(Qgrid)-2 do if (mtr[i, j] lt mtr[i-1, j]) && (mtr[i, j] lt mtr[i+1, j]) then nmin+=1
    if nmin ne 1 then begin
     if exist(loud) then print, '*** More than one local minimum at ', freqList[j], ' GHz ***'
     badf[j]=1
    endif
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
    
   mtr_b=min(mtr[*, j], ib)
   mtr_a=mtr[ib-1, j]
   mtr_c=mtr[ib+1, j]
   Qa=Qgrid[ib-1]
   Qb=Qgrid[ib]
   Qc=Qgrid[ib+1]

   forward_function IRatio
   
   if ((Qc-Qa)/(Qc+Qa)) lt acc then done=1 else begin
    QxG=((Qc-Qb) gt (Qb-Qa)) ? Qb+(Qc-Qb)*(1d0-1d0/G) : $
                               Qb-(Qb-Qa)*(1d0-1d0/G)
    QxB=Qb-0.5*((Qb-Qa)^2*(mtr_b-mtr_c)-(Qb-Qc)^2*(mtr_b-mtr_a))/$
               ((Qb-Qa)*(mtr_b-mtr_c)-(Qb-Qc)*(mtr_b-mtr_a))   
             
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
    print, 'Image metrics:   ', mtr_a, mtr_b, mtr_c
    print, 'Computing images for Q0=', Qx, st
   
    coronaparms=DefineCoronaParams(Tbase, nbase, Qx, a, b, force_isothermal=iso)
    outspace=ReserveOutputSpace(simbox)
    
    r=call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace)
                    
    ConvertToMaps, outspace, simbox, model, modImaps, modVmaps
    obj_destroy, modVmaps
    
    chi_x=dblarr(Nfreq)
    chiVar_x=dblarr(Nfreq)
    rho_x=dblarr(Nfreq)
    rhoVar_x=dblarr(Nfreq)
    eta_x=dblarr(Nfreq)
    etaVar_x=dblarr(Nfreq)        
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
     ExtractSubmap, _obsI, modI, dx, dy, _obsI.id, obsI
     obsI=create_struct('shiftX', dx, $
                        'shiftY', dy, $
                        obsI)
     ExtractSubmap, _obsSigma, modI, dx, dy, _obsI.id+' sigma', obsSigma
     obsX->setmap, k, obsI
    
     obsMax=GetSmoothedMax(obsI, obsInfo.sx[k], obsInfo.sy[k])
     modMax=max(modI.data)    
    
     u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)), ms)
     Iobs_x[k]=total(obsI.data[u])*obsI.dx*obsI.dy
     Imod_x[k]=total(modI.data[u])*modI.dx*modI.dy
     CC_x[k]=c_correlate(obsI.data[u], modI.data[u], 0)
          
     chi_x[k]=mean(((modI.data[u]-obsI.data[u])/obsSigma.data[u])^2)    ;\chi^2
     chiVar_x[k]=variance((modI.data[u]-obsI.data[u])/obsSigma.data[u]) ;\chi^2 corrected
     eta_x[k]=mean(((modI.data[u]-obsI.data[u])/mean(obsI.data[u]))^2)    ;\eta^2
     etaVar_x[k]=variance((modI.data[u]-obsI.data[u])/mean(obsI.data[u])) ;\eta^2 corrected
     rho_x[k]=mean((modI.data[u]/obsI.data[u]-1)^2)    ;\rho^2
     rhoVar_x[k]=variance(modI.data[u]/obsI.data[u]-1) ;\rho^2 corrected  
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
    
    rho2=dblarr(NQ+1, Nfreq)
    rho2[0 : l, *]=rho[0 : l, *]
    rho2[l+1, *]=rho_x
    rho2[l+2 : NQ, *]=rho[l+1 : NQ-1, *]
    rho=rho2
    
    rhoVar2=dblarr(NQ+1, Nfreq)
    rhoVar2[0 : l, *]=rhoVar[0 : l, *]
    rhoVar2[l+1, *]=rhoVar_x
    rhoVar2[l+2 : NQ, *]=rhoVar[l+1 : NQ-1, *]
    rhoVar=rhoVar2        
    
    eta2=dblarr(NQ+1, Nfreq)
    eta2[0 : l, *]=eta[0 : l, *]
    eta2[l+1, *]=eta_x
    eta2[l+2 : NQ, *]=eta[l+1 : NQ-1, *]
    eta=eta2
    
    etaVar2=dblarr(NQ+1, Nfreq)
    etaVar2[0 : l, *]=etaVar[0 : l, *]
    etaVar2[l+1, *]=etaVar_x
    etaVar2[l+2 : NQ, *]=etaVar[l+1 : NQ-1, *]
    etaVar=etaVar2        
    
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
    
    case metric of
     'chi': mtr=chi
     'rho': mtr=rho    
     'eta': mtr=eta   
    endcase                          
   endelse
   
   step+=1
  endwhile
  
  mtr_min=min(mtr[*, j], ib)
  chiArr[j]=chi[ib, j]
  chiVarArr[j]=chiVar[ib, j]
  rhoArr[j]=rho[ib, j]
  rhoVarArr[j]=rhoVar[ib, j]
  etaArr[j]=eta[ib, j]
  etaVarArr[j]=etaVar[ib, j]    
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
  ExtractSubmap, _obsI, m, 0, 0, _obsI.id, obsI
  obsImageArr->setmap, j, obsI
 endelse
 
 allQ=Qgrid
 allMetrics=mtr
end

pro FindBestFitQ, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, a, b, Qstart, iso, $
                  bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                  IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                  freqList, allQ, allMetrics, modImageConvArr, obsImageArr, thr=thr, metric=metric, $
                  noMultiFreq=noMultiFreq, Qstep=Qstep, loud=loud
 forward_function MakeSimulationBox

 if ~keyword_set(noMultiFreq) then $
  FindBestFitQmf, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, a, b, Qstart, iso, $
                  bestQarr, chiArr, chiVarArr, rhoArr, rhoVarArr, etaArr, etaVarArr, $
                  IobsArr, ImodArr, CCarr, modImageArr, modFlagArr, $
                  freqList, allQ, allMetrics, modImageConvArr, obsImageArr, thr=thr, metric=metric, $
                  Qstep=Qstep, loud=loud $
 else begin
  Nfreq=obsInfo.Nfreq
  
  bestQarr=dblarr(Nfreq)
  chiArr=dblarr(Nfreq)
  chiVarArr=dblarr(Nfreq)
  rhoArr=dblarr(Nfreq)
  rhoVarArr=dblarr(Nfreq) 
  etaArr=dblarr(Nfreq)
  etaVarArr=dblarr(Nfreq)  
  IobsArr=dblarr(Nfreq)
  ImodArr=dblarr(Nfreq)
  CCarr=dblarr(Nfreq) 
  modImageArr=obj_new('map')
  modImageConvArr=obj_new('map')
  obsImageArr=obj_new('map')
  modFlagArr=lonarr(Nfreq, 6)
  
  for i=0, Nfreq-1 do begin
   simbox_loc=MakeSimulationBox(simbox.xc, simbox.yc, simbox.dx, simbox.dy, simbox.Nx, simbox.Ny, ObsInfo.freq[i]) 
   obsImaps_loc=obj_new('map')
   m=obsImaps.getmap(i)
   obsImaps_loc->setmap, 0, m
   obsSImaps_loc=obj_new('map')
   m=obsSImaps.getmap(i)
   obsSImaps_loc->setmap, 0, m
   psf_loc=obsInfo.psf
   psf_loc=psf_loc[i]
   psf_loc=list(psf_loc)
   obsInfo_loc={id: ' ', $
                Nfreq: 1, freq: obsInfo.freq[i], sx: obsInfo.sx[i], sy: obsInfo.sy[i], $
                psf_dx: obsInfo.psf_dx[i], psf_dy: obsInfo.psf_dy[i], psf: psf_loc}
                
   FindBestFitQmf, libname, model, ebtel, simbox_loc, obsImaps_loc, obsSImaps_loc, obsInfo_loc, a, b, Qstart, iso, $
                   bestQarr_loc, chiArr_loc, chiVarArr_loc, rhoArr_loc, rhoVarArr_loc, etaArr_loc, etaVarArr_loc, $
                   IobsArr_loc, ImodArr_loc, CCarr_loc, modImageArr_loc, modFlagArr_loc, $
                   freqList_loc, allQ_loc, allMetrics_loc, modImageConvArr_loc, obsImageArr_loc, thr=thr, $
                   metric=metric, Qstep=Qstep, loud=loud
                   
   bestQarr[i]=bestQarr_loc
   chiArr[i]=chiArr_loc
   chiVarArr[i]=chiVarArr_loc
   rhoArr[i]=rhoArr_loc
   rhoVarArr[i]=rhoVarArr_loc 
   etaArr[i]=etaArr_loc
   etaVarArr[i]=etaVarArr_loc  
   IobsArr[i]=IobsArr_loc
   ImodArr[i]=ImodArr_loc
   CCarr[i]=CCarr_loc 
   m=modImageArr_loc.getmap(0)
   modImageArr->setmap, i, m
   m=modImageConvArr_loc.getmap(0)
   modImageConvArr->setmap, i, m
   m=obsImageArr_loc.getmap(0)
   obsImageArr->setmap, i, m
   modFlagArr[i, *]=modFlagArr[0, *]
  endfor
  
  freqList=obsInfo.freq
  allQ=0
  allMetrics=0
 endelse
end                  