pro LoadSRHlow, obsDir, obsId, Imap, Smap, Info, rm
 SRHid=['SRH0306_CH00_2.80GHz', $
        'SRH0306_CH01_3.10GHz', $
        'SRH0306_CH02_3.40GHz', $
        'SRH0306_CH03_3.90GHz', $
        'SRH0306_CH04_4.70GHz', $
        'SRH0306_CH05_5.60GHz']
        
 Nfreq=n_elements(SRHid)
 freq=dblarr(Nfreq)
 sx=dblarr(Nfreq)
 sy=dblarr(Nfreq)
 rho=dblarr(Nfreq)
 
 Imap=obj_new('map')
 Smap=obj_new('map')
 
 for i=0, Nfreq-1 do begin
  fname=file_search(obsDir+SRHid[i]+'_I_'+obsId+'*.fits')
  mreadfits, fname, index, data, /silent
  freq[i]=index.obs_d$freq
  sx[i]=index.beam_sx
  sy[i]=index.beam_sy
  rho[i]=index.beam_rho
  
  ProcessSRHimages, index, data, sigma, rm
  
  index2map, index, data, _obsI
  index2map, index, sigma, _obsSigma
  
  Imap->setmap, i, _obsI
  Smap->setmap, i, _obsSigma
 endfor
 
 Info={id: 'low', $
       Nfreq: Nfreq, freq: freq, sx: sx, sy: sy, rho: rho}
end