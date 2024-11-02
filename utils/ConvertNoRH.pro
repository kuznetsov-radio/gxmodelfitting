pro ConvertNoRH, InDir, InFiles
;This program converts the Nobeyama Radioheliograph data (fits format, no file extension by default) 
;into .sav reference files accepted by the CHMP routines.
;Input parameters:
; InDir - the directory where the input files are located.
; InFiles - array of the input file names.
;Output: the program creates .sav files in the same InDir directory, with the same names as specified by InFiles,
;but with '.sav' extension added.

 Nfreq=n_elements(InFiles)
 
 for i=0, Nfreq-1 do begin
  ref=obj_new('map')
  
  fname=InDir+InFiles[i]
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
  
  fname1=InDir+InFiles[i]+'.sav'
  save, ref, filename=fname1, /compress
 endfor
end