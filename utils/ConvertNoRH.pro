pro ConvertNoRH
 obsDir='C:\MCloud\CoronalMW\AR-SRH\Data\NoRH\'
 flist=['ifa140202_022005', 'ifz140202_022005']
        
 Nfreq=n_elements(flist)
 
 for i=0, Nfreq-1 do begin
  ref=obj_new('map')
  
  fname=obsDir+flist[i]
  fits2map, fname, obsI
  norh_rd_img, fname, index, data
  
  obsI=create_struct('freq', double(index.norh.obs_freq), obsI)
  obsSigma=obsI
  
  beam=norh_beam(index)
  psf=make_map(beam, dx=index.norh.sec_per_pix, dy=index.norh.sec_per_pix, xc=0, yc=0)
  
  parms=BeamFitNoRH(fname, x1, y1, beam1, /quiet)
  psf=create_struct('a_beam', parms[0], $
                    'b_beam', parms[1], $
                    'phi_beam', parms[2], $
                    'freq', double(index.norh.obs_freq), psf)
 
  ref->setmap, 0, obsI
  ref->setmap, 1, obsSigma
  ref->setmap, 2, psf
  
  fname1=obsDir+flist[i]+'.sav'
  save, ref, filename=fname1, /compress
 endfor
end