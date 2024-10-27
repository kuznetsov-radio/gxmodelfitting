pro FindBestFitQmf, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ ;input
                    a, b, Qstart, Qstep, iso, thr, metric, fixed_shifts, xy_shift, $ ;input
                    freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ ;output
                    ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ ;output
                    obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ ;output
                    modFlagArr, allQ, allMetrics, $ ;extra output
                    loud=loud
 forward_function GetSmoothedMax, DefineCoronaParms, ReserveOutputSpace, ConvertToMaps
                       
 acc=1d-2 ;desired accuracy
 Tbase=1d6 ;analytical corona temperature
 nbase=1d8 ;analytical corona base density
 
 Nfreq=obsInfo.Nfreq
 freqList=obsInfo.freq
 
 RATAN_on=obsInfo.id eq 'RATAN'
 
 Qgrid=[Qstart]
 fdone=[0]
 modImages=objarr(1)
 modImagesConv=objarr(1)
 obsImages=objarr(1)
 obsImagesSigma=objarr(1)
 flags=lonarr(1, 6)
 chi=dblarr(1, Nfreq)
 rho=dblarr(1, Nfreq)
 eta=dblarr(1, Nfreq)  
 ItotalObs=dblarr(1, Nfreq)
 ItotalMod=dblarr(1, Nfreq)
 ImaxObs=dblarr(1, Nfreq)
 ImaxMod=dblarr(1, Nfreq)
 IthrObs=dblarr(1, Nfreq)
 IthrMod=dblarr(1, Nfreq) 
 CC=dblarr(1, Nfreq)
 
 bestQarr=dblarr(Nfreq)
 chiArr=dblarr(Nfreq)
 rhoArr=dblarr(Nfreq)
 etaArr=dblarr(Nfreq)
 ItotalObsArr=dblarr(Nfreq)
 ItotalModArr=dblarr(Nfreq)
 ImaxObsArr=dblarr(Nfreq)
 ImaxModArr=dblarr(Nfreq)
 IthrObsArr=dblarr(Nfreq)
 IthrModArr=dblarr(Nfreq)  
 CCarr=dblarr(Nfreq) 
 modImageArr=obj_new('map')
 modImageConvArr=RATAN_on ? list() : obj_new('map')
 obsImageArr=RATAN_on ? list() : obj_new('map')
 obsImageSigmaArr=RATAN_on ? list() : obj_new('map')
 modFlagArr=lonarr(Nfreq, 6)
 
 bestQarr[*]=!values.d_NaN
 chiArr[*]=!values.d_NaN
 rhoArr[*]=!values.d_NaN
 etaArr[*]=!values.d_NaN
 ItotalObsArr[*]=!values.d_NaN
 ItotalModArr[*]=!values.d_NaN
 ImaxObsArr[*]=!values.d_NaN
 ImaxModArr[*]=!values.d_NaN
 IthrObsArr[*]=!values.d_NaN
 IthrModArr[*]=!values.d_NaN  
 CCarr[*]=!values.d_NaN 
 
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
                    
   ConvertToMaps, outspace, simbox, model, modImaps, modVmaps, flux=RATAN_on
   if ~RATAN_on then obj_destroy, modVmaps
   modImages[i]=modImaps
   flags[i, *]=outspace.flagsCorona
   
   if RATAN_on then begin
    modX=list()
    obsX=list()
    obsSX=list()
    
    modL=obj_new('map')
    modR=obj_new('map')
    for j=0, Nfreq-1 do begin
     modI=modImaps.getmap(j)
     modV=modVmaps.getmap(j)

     mR=modI
     mR.data=(modI.data+modV.data)/2
     modR->setmap, j, mR
     mL=modI
     mL.data=(modI.data-modV.data)/2
     modL->setmap, j, mL  
    endfor  
    
    asu_gxm_maplist2scans, modL, mScans_L, mXarc, out_index=out_index
    asu_gxm_maplist2scans, modR, mScans_R, mXarc, out_index=out_index
    mXarc-=out_index[0].xc
    mScans=mScans_L+mScans_R
    obj_destroy, modL
    obj_destroy, modR
    
    for j=0, Nfreq-1 do begin
     modI={flux: reform(mScans[*, j]), $
           x: mXarc, $
           freq: freqList[j]}
     modX.add, modI        
 
     _obsI=obsImaps[j]
     _obsSigma=obsSImaps[j]
     
     if fixed_shifts then dx=xy_shift else FindShift1D, _obsI, modI, dx
     ExtractSubProfile, _obsI, modI, dx, obsI
     ExtractSubProfile, _obsSigma, modI, dx, obsSigma
     obsX.add, obsI
     obsSX.add, obsSI
     
     obsMax=max(obsI.flux)
     modMax=max(modI.flux)    
     ImaxObs[i, j]=obsMax
     ImaxMod[i, j]=modMax
     IthrObs[i, j]=obsMax*thr
     IthrMod[i, j]=modMax*thr
     
     u=where(obsI.flux gt (obsMax*thr), ms)
     maskObs=1d0*ms/n_elements(obsI.flux)
    
     u=where(modI.flux gt (modMax*thr), ms)
     maskMod=1d0*ms/n_elements(modI.flux)    
    
     u=where((obsI.flux gt (obsMax*thr)) or (modI.flux gt (modMax*thr)), ms)
     ItotalObs[i, j]=total(obsI.flux[u])*(obsI.x[1]-obsI.x[0])
     ItotalMod[i, j]=total(modI.flux[u])*(modI.x[1]-modI.x[0])
     CC[i, j]=c_correlate(obsI.flux[u], modI.flux[u], 0)
    
     chi[i, j]=mean(((modI.flux[u]-obsI.flux[u])/obsSigma.flux[u])^2)    ;\chi^2
     rho[i, j]=mean((modI.flux[u]/obsI.flux[u]-1)^2)    ;\rho^2
     eta[i, j]=mean(((modI.flux[u]-obsI.flux[u])/mean(obsI.flux[u]))^2)    ;\eta^2     
    endfor
    
    obj_destroy, modVmaps
   endif else begin
    modX=obj_new('map')
    obsX=obj_new('map')
    obsSX=obj_new('map')
   
    for j=0, Nfreq-1 do begin
     modI=modImaps.getmap(j)
     _obsI=obsImaps.getmap(j)
     _obsSigma=obsSImaps.getmap(j)
    
     MakeLocalBeam, obsInfo, j, modI.dx, modI.dy, beam
     FixLocalBeam, beam, simbox.Nx, simbox.Ny
     modI.data=convol_fft(modI.data, beam)
     modX->setmap, j, modI
    
     if fixed_shifts then begin
      dx=xy_shift[0]
      dy=xy_shift[1]
     endif else FindShift, _obsI, modI, dx, dy
     ExtractSubmap, _obsI, modI, dx, dy, _obsI.id, obsI
     obsI=create_struct('shiftX', dx, $
                        'shiftY', dy, $
                        obsI)
     ExtractSubmap, _obsSigma, modI, dx, dy, _obsI.id+' sigma', obsSigma
     obsX->setmap, j, obsI
     obsSX->setmap, j, obsSigma
    
     obsMax=GetSmoothedMax(obsI, obsInfo.sx[j], obsInfo.sy[j])
     modMax=max(modI.data)    
     ImaxObs[i, j]=obsMax
     ImaxMod[i, j]=modMax
     IthrObs[i, j]=obsMax*thr
     IthrMod[i, j]=modMax*thr
    
     u=where(obsI.data gt (obsMax*thr), ms)
     maskObs=1d0*ms/n_elements(obsI.data)
    
     u=where(modI.data gt (modMax*thr), ms)
     maskMod=1d0*ms/n_elements(modI.data)    
    
     u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)), ms)
     ItotalObs[i, j]=total(obsI.data[u])*obsI.dx*obsI.dy
     ItotalMod[i, j]=total(modI.data[u])*modI.dx*modI.dy
     CC[i, j]=c_correlate(obsI.data[u], modI.data[u], 0)
    
     chi[i, j]=mean(((modI.data[u]-obsI.data[u])/obsSigma.data[u])^2)    ;\chi^2
     rho[i, j]=mean((modI.data[u]/obsI.data[u]-1)^2)    ;\rho^2
     eta[i, j]=mean(((modI.data[u]-obsI.data[u])/mean(obsI.data[u]))^2)    ;\eta^2
        
     if (maskMod gt 0.99) || ((maskMod/maskObs) gt 4) then begin
      if exist(loud) then print, '*** GR contribution is too low at ', freqList[j], ' GHz ***'
      chi[i, j]=!values.d_NaN
      rho[i, j]=!values.d_NaN
      eta[i, j]=!values.d_NaN
     endif
    endfor
   endelse
   
   case metric of
    'chi': mtr=chi
    'rho': mtr=rho    
    'eta': mtr=eta   
   endcase  
   
   modImagesConv[i]=modX
   obsImages[i]=obsX
   obsImagesSigma[i]=obsSX
 
   fdone[i]=1
  endif
  
  if NQ eq 1 then begin
   aw=total(sgn(ItotalObs-ItotalMod))
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
   obsImagesSigma=[obj_new(), obsImagesSigma]
      
   flags2=lonarr(NQ+1, 6)
   flags2[1 : NQ, *]=flags
   flags=flags2
 
   chi2=dblarr(NQ+1, Nfreq)
   chi2[1 : NQ, *]=chi
   chi=chi2
   
   rho2=dblarr(NQ+1, Nfreq)
   rho2[1 : NQ, *]=rho
   rho=rho2
   
   eta2=dblarr(NQ+1, Nfreq)
   eta2[1 : NQ, *]=eta
   eta=eta2
   
   ItotalObs2=dblarr(NQ+1, Nfreq)
   ItotalObs2[1 : NQ, *]=ItotalObs
   ItotalObs=ItotalObs2
   
   ItotalMod2=dblarr(NQ+1, Nfreq)
   ItotalMod2[1 : NQ, *]=ItotalMod
   ItotalMod=ItotalMod2
   
   ImaxObs2=dblarr(NQ+1, Nfreq)
   ImaxObs2[1 : NQ, *]=ImaxObs
   ImaxObs=ImaxObs2
   
   ImaxMod2=dblarr(NQ+1, Nfreq)
   ImaxMod2[1 : NQ, *]=ImaxMod
   ImaxMod=ImaxMod2   

   IthrObs2=dblarr(NQ+1, Nfreq)
   IthrObs2[1 : NQ, *]=IthrObs
   IthrObs=IthrObs2
   
   IthrMod2=dblarr(NQ+1, Nfreq)
   IthrMod2[1 : NQ, *]=IthrMod
   IthrMod=IthrMod2
   
   CC2=dblarr(NQ+1, Nfreq)
   CC2[1 : NQ, *]=CC
   CC=CC2      
  endif else if aw gt 0 then begin
   Qgrid=[Qgrid, max(Qgrid)*Qstep]
   fdone=[fdone, 0]
   modImages=[modImages, obj_new()]
   modImagesConv=[modImagesConv, obj_new()]
   obsImages=[obsImages, obj_new()]
   obsImagesSigma=[obsImagesSigma, obj_new()]
   
   flags2=lonarr(NQ+1, 6)
   flags2[0 : NQ-1, *]=flags
   flags=flags2
   
   chi2=dblarr(NQ+1, Nfreq)
   chi2[0 : NQ-1, *]=chi
   chi=chi2
   
   rho2=dblarr(NQ+1, Nfreq)
   rho2[0 : NQ-1, *]=rho
   rho=rho2
   
   eta2=dblarr(NQ+1, Nfreq)
   eta2[0 : NQ-1, *]=eta
   eta=eta2
   
   ItotalObs2=dblarr(NQ+1, Nfreq)
   ItotalObs2[0 : NQ-1, *]=ItotalObs
   ItotalObs=ItotalObs2   
   
   ItotalMod2=dblarr(NQ+1, Nfreq)
   ItotalMod2[0 : NQ-1, *]=ItotalMod
   ItotalMod=ItotalMod2
   
   ImaxObs2=dblarr(NQ+1, Nfreq)
   ImaxObs2[0 : NQ-1, *]=ImaxObs
   ImaxObs=ImaxObs2   
   
   ImaxMod2=dblarr(NQ+1, Nfreq)
   ImaxMod2[0 : NQ-1, *]=ImaxMod
   ImaxMod=ImaxMod2            

   IthrObs2=dblarr(NQ+1, Nfreq)
   IthrObs2[0 : NQ-1, *]=IthrObs
   IthrObs=IthrObs2   
   
   IthrMod2=dblarr(NQ+1, Nfreq)
   IthrMod2[0 : NQ-1, *]=IthrMod
   IthrMod=IthrMod2   
   
   CC2=dblarr(NQ+1, Nfreq)
   CC2[0 : NQ-1, *]=CC
   CC=CC2        
  endif else done=1
 endwhile
 
 ;-----------------------------------------------------------------
 
 G=(1d0+sqrt(5d0))/2
 
 badf=intarr(Nfreq)
 
 for j=0, Nfreq-1 do begin
  u=where(~finite(mtr[*, j]), k)
  
  if k gt 0 then badf[j]=1 else begin
   mtr_min=min(mtr[*, j], k)
   
   if (k eq 0) || (k eq (n_elements(Qgrid)-1)) then badf[j]=1 else begin
    nmin=0
    for i=1, n_elements(Qgrid)-2 do if (mtr[i, j] lt mtr[i-1, j]) && (mtr[i, j] lt mtr[i+1, j]) then nmin+=1
    if nmin ne 1 then begin
     if exist(loud) then print, (nmin gt 1) ? '*** More than one local minimum at ' : $
                                              '*** No local minima at ', freqList[j], ' GHz ***'
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
                    
    ConvertToMaps, outspace, simbox, model, modImaps, modVmaps, flux=RATAN_on
    if ~RATAN_on then obj_destroy, modVmaps
        
    chi_x=dblarr(Nfreq)
    rho_x=dblarr(Nfreq)
    eta_x=dblarr(Nfreq)
    ItotalObs_x=dblarr(Nfreq)
    ItotalMod_x=dblarr(Nfreq)
    ImaxObs_x=dblarr(Nfreq)
    ImaxMod_x=dblarr(Nfreq)
    IthrObs_x=dblarr(Nfreq)
    IthrMod_x=dblarr(Nfreq)        
    CC_x=dblarr(Nfreq)

    if RATAN_on then begin
     modX=list()
     obsX=list()
     obsSX=list()
    
     modL=obj_new('map')
     modR=obj_new('map')
     for k=0, Nfreq-1 do begin
      modI=modImaps.getmap(k)
      modV=modVmaps.getmap(k)

      mR=modI
      mR.data=(modI.data+modV.data)/2
      modR->setmap, k, mR
      mL=modI
      mL.data=(modI.data-modV.data)/2
      modL->setmap, k, mL  
     endfor  
    
     asu_gxm_maplist2scans, modL, mScans_L, mXarc, out_index=out_index
     asu_gxm_maplist2scans, modR, mScans_R, mXarc, out_index=out_index
     mXarc-=out_index[0].xc
     mScans=mScans_L+mScans_R
     obj_destroy, modL
     obj_destroy, modR
    
     for k=0, Nfreq-1 do begin
      modI={flux: reform(mScans[*, k]), $
            x: mXarc, $
            freq: freqList[k]}
      modX.add, modI        
 
      _obsI=obsImaps[k]
      _obsSigma=obsSImaps[k]
     
      if fixed_shifts then dx=xy_shift else FindShift1D, _obsI, modI, dx
      ExtractSubProfile, _obsI, modI, dx, obsI
      ExtractSubProfile, _obsSigma, modI, dx, obsSigma
      obsX.add, obsI
      obsSX.add, obsSI
     
      obsMax=max(obsI.flux)
      modMax=max(modI.flux)    
      ImaxObs_x[k]=obsMax
      ImaxMod[k]=modMax
      IthrObs_x[k]=obsMax*thr
      IthrMod_x[k]=modMax*thr
     
      u=where((obsI.flux gt (obsMax*thr)) or (modI.flux gt (modMax*thr)), ms)
      ItotalObs_x[k]=total(obsI.flux[u])*(obsI.x[1]-obsI.x[0])
      ItotalMod_x[k]=total(modI.flux[u])*(modI.x[1]-modI.x[0])
      CC_x[k]=c_correlate(obsI.flux[u], modI.flux[u], 0)
    
      chi_x[k]=mean(((modI.flux[u]-obsI.flux[u])/obsSigma.flux[u])^2)    ;\chi^2
      rho_x[k]=mean((modI.flux[u]/obsI.flux[u]-1)^2)    ;\rho^2
      eta_x[k]=mean(((modI.flux[u]-obsI.flux[u])/mean(obsI.flux[u]))^2)    ;\eta^2     
     endfor
    
     obj_destroy, modVmaps
    endif else begin    
     modX=obj_new('map')
     obsX=obj_new('map')
     obsSX=obj_new('map')
    
     for k=0, Nfreq-1 do begin
      modI=modImaps.getmap(k)
      _obsI=obsImaps.getmap(k)
      _obsSigma=obsSImaps.getmap(k)
    
      MakeLocalBeam, obsInfo, k, modI.dx, modI.dy, beam
      FixLocalBeam, beam, simbox.Nx, simbox.Ny
      modI.data=convol_fft(modI.data, beam)
      modX->setmap, k, modI
    
      if fixed_shifts then begin
       dx=xy_shift[0]
       dy=xy_shift[1]
      endif else FindShift, _obsI, modI, dx, dy
      ExtractSubmap, _obsI, modI, dx, dy, _obsI.id, obsI
      obsI=create_struct('shiftX', dx, $
                         'shiftY', dy, $
                         obsI)
      ExtractSubmap, _obsSigma, modI, dx, dy, _obsI.id+' sigma', obsSigma
      obsX->setmap, k, obsI
      obsSX->setmap, k, obsSigma
    
      obsMax=GetSmoothedMax(obsI, obsInfo.sx[k], obsInfo.sy[k])
      modMax=max(modI.data)    
      ImaxObs_x[k]=obsMax
      ImaxMod_x[k]=modMax
      IthrObs_x[k]=obsMax*thr
      IthrMod_x[k]=modMax*thr
    
      u=where((obsI.data gt (obsMax*thr)) or (modI.data gt (modMax*thr)), ms)
      ItotalObs_x[k]=total(obsI.data[u])*obsI.dx*obsI.dy
      ItotalMod_x[k]=total(modI.data[u])*modI.dx*modI.dy
      CC_x[k]=c_correlate(obsI.data[u], modI.data[u], 0)
          
      chi_x[k]=mean(((modI.data[u]-obsI.data[u])/obsSigma.data[u])^2)    ;\chi^2
      eta_x[k]=mean(((modI.data[u]-obsI.data[u])/mean(obsI.data[u]))^2)    ;\eta^2
      rho_x[k]=mean((modI.data[u]/obsI.data[u]-1)^2)    ;\rho^2
     endfor 
    endelse
    
    l=(Qx gt Qb) ? ib : ib-1
    
    Qgrid=[Qgrid[0 : l], Qx, Qgrid[l+1 : NQ-1]]
    modImages=[modImages[0 : l], modImaps, modImages[l+1 : NQ-1]]
    modImagesConv=[modImagesConv[0 : l], modX, modImagesConv[l+1 : NQ-1]]
    obsImages=[obsImages[0 : l], obsX, obsImages[l+1 : NQ-1]]
    obsImagesSigma=[obsImagesSigma[0 : l], obsSX, obsImagesSigma[l+1 : NQ-1]]
    
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
    
    rho2=dblarr(NQ+1, Nfreq)
    rho2[0 : l, *]=rho[0 : l, *]
    rho2[l+1, *]=rho_x
    rho2[l+2 : NQ, *]=rho[l+1 : NQ-1, *]
    rho=rho2
    
    eta2=dblarr(NQ+1, Nfreq)
    eta2[0 : l, *]=eta[0 : l, *]
    eta2[l+1, *]=eta_x
    eta2[l+2 : NQ, *]=eta[l+1 : NQ-1, *]
    eta=eta2
    
    ItotalObs2=dblarr(NQ+1, Nfreq)
    ItotalObs2[0 : l, *]=ItotalObs[0 : l, *]
    ItotalObs2[l+1, *]=ItotalObs_x
    ItotalObs2[l+2 : NQ, *]=ItotalObs[l+1 : NQ-1, *]
    ItotalObs=ItotalObs2
    
    ItotalMod2=dblarr(NQ+1, Nfreq)
    ItotalMod2[0 : l, *]=ItotalMod[0 : l, *]
    ItotalMod2[l+1, *]=ItotalMod_x
    ItotalMod2[l+2 : NQ, *]=ItotalMod[l+1 : NQ-1, *]
    ItotalMod=ItotalMod2

    ImaxObs2=dblarr(NQ+1, Nfreq)
    ImaxObs2[0 : l, *]=ImaxObs[0 : l, *]
    ImaxObs2[l+1, *]=ImaxObs_x
    ImaxObs2[l+2 : NQ, *]=ImaxObs[l+1 : NQ-1, *]
    ImaxObs=ImaxObs2
    
    ImaxMod2=dblarr(NQ+1, Nfreq)
    ImaxMod2[0 : l, *]=ImaxMod[0 : l, *]
    ImaxMod2[l+1, *]=ImaxMod_x
    ImaxMod2[l+2 : NQ, *]=ImaxMod[l+1 : NQ-1, *]
    ImaxMod=ImaxMod2

    IthrObs2=dblarr(NQ+1, Nfreq)
    IthrObs2[0 : l, *]=IthrObs[0 : l, *]
    IthrObs2[l+1, *]=IthrObs_x
    IthrObs2[l+2 : NQ, *]=IthrObs[l+1 : NQ-1, *]
    IthrObs=IthrObs2
    
    IthrMod2=dblarr(NQ+1, Nfreq)
    IthrMod2[0 : l, *]=IthrMod[0 : l, *]
    IthrMod2[l+1, *]=IthrMod_x
    IthrMod2[l+2 : NQ, *]=IthrMod[l+1 : NQ-1, *]
    IthrMod=IthrMod2
    
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
  rhoArr[j]=rho[ib, j]
  etaArr[j]=eta[ib, j]
  ItotalObsArr[j]=ItotalObs[ib, j]
  ItotalModArr[j]=ItotalMod[ib, j]
  ImaxObsArr[j]=ImaxObs[ib, j]
  ImaxModArr[j]=ImaxMod[ib, j]
  IthrObsArr[j]=IthrObs[ib, j]
  IthrModArr[j]=IthrMod[ib, j]    
  CCarr[j]=CC[ib, j]  
  bestQarr[j]=Qgrid[ib]
  modFlagArr[j, *]=flags[ib, *]
  
  bmap=modImages[ib]
  m=bmap.getmap(j)
  modImageArr->setmap, j, m
  
  if RATAN_on then begin
   bmap=modImagesConv[ib]
   m=bmap[j]
   modImageConvArr.add, m
   bmap=obsImages[ib]
   m=bmap[j]
   obsImageArr.add, m
   bmap=obsImagesSigma[ib]
   m=bmap[j]
   obsImageSigmaArr.add, m
  endif else begin
   bmap=modImagesConv[ib]
   m=bmap.getmap(j)
   modImageConvArr->setmap, j, m
   bmap=obsImages[ib]
   m=bmap.getmap(j)
   obsImageArr->setmap, j, m
   bmap=obsImagesSigma[ib]
   m=bmap.getmap(j)
   obsImageSigmaArr->setmap, j, m
  endelse     
 endif else begin
  bmap=modImages[0]
  m=bmap.getmap(j)
  m.data[*]=0
  modImageArr->setmap, j, m
  
  if RATAN_on then begin
   bmap=modImagesConv[0]
   m=bmap[j]
   m.flux[*]=0
   modImageConvArr.add, m
   _obsI=obsImaps[j]
   ExtractSubProfile, _obsI, m, 0, obsI
   obsImageArr.add, obsI
   _obsSI=obsSImaps[j]
   ExtractSubProfile, _obsSI, m, 0, obsSI
   obsImageSigmaArr.add, obsSI
  endif else begin
   modImageConvArr->setmap, j, m
   _obsI=obsImaps.getmap(j)
   ExtractSubmap, _obsI, m, 0, 0, _obsI.id, obsI
   obsImageArr->setmap, j, obsI
   _obsSI=obsSImaps.getmap(j)
   ExtractSubmap, _obsSI, m, 0, 0, _obsSI.id, obsSI
   obsImageSigmaArr->setmap, j, obsSI  
  endelse
 endelse
 
 allQ=Qgrid
 allMetrics=mtr
end

pro FindBestFitQ, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ ;input
                  a, b, Qstart, Qstep, iso, thr, metric, MultiFreq_on, fixed_shifts, xy_shift, $ ;input
                  freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $ ;output
                  ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ ;output
                  obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ ;output
                  modFlagArr, allQ, allMetrics, $ ;extra output
                  loud=loud
 forward_function MakeSimulationBox

 if MultiFreq_on then $
  FindBestFitQmf, libname, model, ebtel, simbox, obsImaps, obsSImaps, obsInfo, $ 
                  a, b, Qstart, Qstep, iso, thr, metric, fixed_shifts, xy_shift, $         
                  freqList, bestQarr, chiArr, rhoArr, etaArr, CCarr, $           
                  ItotalObsArr, ItotalModArr, ImaxObsArr, ImaxModArr, IthrObsArr, IthrModArr, $ 
                  obsImageArr, obsImageSigmaArr, modImageArr, modImageConvArr, $ 
                  modFlagArr, allQ, allMetrics, loud=loud $
 else begin
  Nfreq=obsInfo.Nfreq
  
  RATAN_on=obsInfo.id eq 'RATAN'
  
  bestQarr=dblarr(Nfreq)
  chiArr=dblarr(Nfreq)
  rhoArr=dblarr(Nfreq)
  etaArr=dblarr(Nfreq)
  ItotalObsArr=dblarr(Nfreq)
  ItotalModArr=dblarr(Nfreq)
  ImaxObsArr=dblarr(Nfreq)
  ImaxModArr=dblarr(Nfreq)
  IthrObsArr=dblarr(Nfreq)
  IthrModArr=dblarr(Nfreq)    
  CCarr=dblarr(Nfreq) 
  modImageArr=obj_new('map')
  modImageConvArr=RATAN_on ? list() : obj_new('map')
  obsImageArr=RATAN_on ? list() : obj_new('map')
  obsImageSigmaArr=RATAN_on ? list() : obj_new('map')
  modFlagArr=lonarr(Nfreq, 6)
  
  for i=0, Nfreq-1 do begin
   simbox_loc=MakeSimulationBox(simbox.xc, simbox.yc, simbox.dx, simbox.dy, simbox.Nx, simbox.Ny, obsInfo.freq[i], $
                                rot=RATAN_on ? obsInfo.rot[i] : 0d0) 
   if RATAN_on then begin
    obsImaps_loc=list()
    m=obsImaps[i]
    obsImaps_loc.add, m
    obsSImaps_loc=list()
    m=obsSImaps[i]
    obsSImaps_loc.add, m
    obsInfo_loc={id: 'RATAN', Nfreq: 1, freq: obsInfo.freq[i], rot: obsInfo.rot[i]}
   endif else begin                             
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
   endelse             
                
   FindBestFitQmf, libname, model, ebtel, simbox_loc, obsImaps_loc, obsSImaps_loc, obsInfo_loc, $ 
                   a, b, Qstart, Qstep, iso, thr, metric, fixed_shifts, xy_shift, $         
                   freqList_loc, bestQarr_loc, chiArr_loc, rhoArr_loc, etaArr_loc, CCarr_loc, $           
                   ItotalObsArr_loc, ItotalModArr_loc, ImaxObsArr_loc, ImaxModArr_loc, IthrObsArr_loc, IthrModArr_loc, $ 
                   obsImageArr_loc, obsImageSigmaArr_loc, modImageArr_loc, modImageConvArr_loc, $ 
                   modFlagArr_loc, allQ_loc, allMetrics_loc, loud=loud
                   
   bestQarr[i]=bestQarr_loc
   chiArr[i]=chiArr_loc
   rhoArr[i]=rhoArr_loc 
   etaArr[i]=etaArr_loc  
   ItotalObsArr[i]=ItotalObsArr_loc
   ItotalModArr[i]=ItotalModArr_loc
   ImaxObsArr[i]=ImaxObsArr_loc
   ImaxModArr[i]=ImaxModArr_loc
   IthrObsArr[i]=IthrObsArr_loc
   IthrModArr[i]=IthrModArr_loc   
   CCarr[i]=CCarr_loc 
   m=modImageArr_loc.getmap(0)
   modImageArr->setmap, i, m
   if RATAN_on then begin
    m=modImageConvArr_loc[0]
    modImageConvArr.add, m
    m=obsImageArr_loc[0]
    obsImageArr.add, m
    m=obsImageSigmaArr_loc[0]
    obsImageSigmaArr.add, m
   endif else begin
    m=modImageConvArr_loc.getmap(0)
    modImageConvArr->setmap, i, m
    m=obsImageArr_loc.getmap(0)
    obsImageArr->setmap, i, m
    m=obsImageSigmaArr_loc.getmap(0)
    obsImageSigmaArr->setmap, i, m
   endelse   
   modFlagArr[i, *]=modFlagArr[0, *]
  endfor
  
  freqList=obsInfo.freq
  allQ=0
  allMetrics=0
 endelse
end                  