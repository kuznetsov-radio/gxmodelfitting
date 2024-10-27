pro LoadObservations, fnames, Imap, Smap, Info
 Nfreq=n_elements(fnames)
 
 n_RATAN=0
 n_gen=0
 for i=0, Nfreq-1 do begin
  o=obj_new('IDL_Savefile', fnames[i])
  if InSav(o, 'instrument') then begin
   o->restore, 'instrument'
   if instrument eq 'RATAN' then n_RATAN+=1 else n_gen+=1
  endif else n_gen+=1
  obj_destroy, o
 endfor
 
 if n_gen eq Nfreq then begin
  freq1=dblarr(Nfreq)
  sx1=dblarr(Nfreq)
  sy1=dblarr(Nfreq)
  psf_dx1=dblarr(Nfreq)
  psf_dy1=dblarr(Nfreq)
  psf1=list()
  Imap1=obj_new('map')
  Smap1=obj_new('map')
 
  for i=0, Nfreq-1 do begin
   restore, fnames[i], /relaxed_structure_assignment
   obsI=ref->get(0, /map)
   obsSigma=ref->get(1, /map)
   beam=ref->get(2, /map)
  
   freq1[i]=obsI.freq
   sx1[i]=beam.a_beam
   sy1[i]=beam.b_beam
   psf_dx1[i]=beam.dx
   psf_dy1[i]=beam.dy
   psf1.add, beam.data
   Imap1->setmap, i, obsI
   Smap1->setmap, i, obsSigma
  endfor
 
  u=sort(freq1)
 
  freq=freq1[u]
  sx=sx1[u]
  sy=sy1[u]
  psf_dx=psf_dx1[u]
  psf_dy=psf_dy1[u]
 
  psf=list()
  Imap=obj_new('map')
  Smap=obj_new('map')
 
  for i=0, Nfreq-1 do begin
   psf.add, psf1[u[i]]
   Imap->setmap, i, Imap1.getmap(u[i])
   Smap->setmap, i, Smap1.getmap(u[i])
  endfor
 
  Info={id: 'gen2D', $
        Nfreq: Nfreq, freq: freq, sx: sx, sy: sy, $
        psf_dx: psf_dx, psf_dy: psf_dy, psf: psf}
 endif else if n_RATAN eq Nfreq then begin
  freq1=dblarr(Nfreq)
  rot1=dblarr(Nfreq)
  Imap1=list()
  
  for i=0, Nfreq-1 do begin
   restore, fnames[i], /relaxed_structure_assignment
   freq1[i]=ref.freq
   rot1[i]=ref.rot
   Imap1.add, ref
  endfor
  
  u=sort(freq1)
  freq=freq1[u]
  rot=rot1[u]
  
  Imap=list()
  Smap=list()
  for i=0, Nfreq-1 do begin
   Imap.add, Imap1[u[i]]
   Smap.add, Imap1[u[i]]
  endfor
  
  Info={id: 'RATAN', $
        Nfreq: Nfreq, freq: freq, rot: rot}
 endif else begin
  print, 'Data at all frequencies must be of the same type (either 1D or 2D)!'
  stop
 endelse      
end