pro MultiScanABlow
 imgDir='C:\MCloud\CoronalMW\AR-SRH\Data\Images\low\'
 obsDir='C:\MCloud\CoronalMW\AR-SRH\Data\SRH\low\'
 obsId='20210603' 
 model=LoadGXmodel('C:\MCloud\CoronalMW\AR-SRH\Data\Models\model20210603.gxm', $
                   230.0, 210.0, 256.0, 256.0, 128, 128)

 LoadSRHlow, obsDir, obsId, obsImaps, obsSImaps, obsInfo, 0
 
 a_arr=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
 b_arr=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
 
 for i=0, n_elements(a_arr)-1 do for j=0, n_elements(b_arr)-1 do begin
  a=a_arr[i]
  b=b_arr[j]
  print, 'Computing the best fit Q for a=', a, ', b=', b 
  
  Q0=(b gt 0) ? 1e-3 : 1e-3/16
 
  fname=imgDir+'bestQ_'+obsId+'_a'+string(a, format='(F5.2)')+'_b'+string(b, format='(F5.2)')+'.sav'
  
  if file_exist(fname) then begin
   o=obj_new('IDL_Savefile', fname)
   o->Restore, 'Qgrid'
   o->Restore, 'modIgrid'
   obj_destroy, o
   FindBestFitQ, model, obsImaps, obsInfo, a, b, Q0, 0.1, bestQarr, Qgrid, modIgrid, /reprocess
  endif else FindBestFitQ, model, obsImaps, obsInfo, a, b, Q0, 0.1, bestQarr, Qgrid, modIgrid

  save, obsId, a, b, bestQarr, Qgrid, modIgrid, filename=fname, /compress
 
  for k=0, n_elements(Qgrid)-1 do obj_destroy, modIgrid[k]
 endfor 
 
 obj_destroy, model
 obj_destroy, obsImaps
 obj_destroy, obsSImaps
end
