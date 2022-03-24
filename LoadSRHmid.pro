pro LoadSRHmid, obsFile, Imap, Smap, Info, rm
 restore, obsfile, /relaxed_structure_assignment
 
 Nfreq=n_elements(index_I)
 freq=dblarr(Nfreq)
 dx=dblarr(Nfreq)
 dy=dblarr(Nfreq)
 sx=dblarr(Nfreq)
 sy=dblarr(Nfreq)
 
 m=index_i[0]
 Ntimes=n_elements(m)
  
 Imap=obj_new('map')
 Smap=obj_new('map')
 
 for i=0, Nfreq-1 do begin
  index=index_I[i, Ntimes/2]
  data=median_I[i]
  freq[i]=index.obs_d$freq
  sx[i]=index.beam_sx
  sy[i]=index.beam_sy
  
  ProcessSRHimages, index, data, sigma, rm
  
  index2map, index, data, _obsI
  index2map, index, sigma, _obsSigma
  
  dx[i]=_obsI.dx
  dy[i]=_obsI.dy
  
  Imap->setmap, i, _obsI
  Smap->setmap, i, _obsSigma
 endfor
 
 Info={id: 'mid', $
       Nfreq: Nfreq, freq: freq, dx: dx, dy: dy, sx: sx, sy: sy, beam: median_psf}
end