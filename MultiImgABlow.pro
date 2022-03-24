pro MultiImgABlow
 imgDir='C:\MCloud\CoronalMW\AR-SRH\Data\Images\low\'
 obsDir='C:\MCloud\CoronalMW\AR-SRH\Data\SRH\low\'
 obsId='20210603' 
 model=LoadGXmodel('C:\MCloud\CoronalMW\AR-SRH\Data\Models\model20210603.gxm', $
                   230.0, 210.0, 256.0, 256.0, 128, 128)

 LoadSRHlow, obsDir, obsId, obsImaps, obsSImaps, obsInfo, 3
 
 a_arr=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
 b_arr=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
 
 for i=0, n_elements(a_arr)-1 do for j=0, n_elements(b_arr)-1 do begin
  a=a_arr[i]
  b=b_arr[j]
  print, 'Computing images and metrics for a=', a, ', b=', b 
 
  fname_in=imgDir+'bestQ_'+obsId+'_a'+string(a, format='(F5.2)')+'_b'+string(b, format='(F5.2)')+'.sav'
  fname_out=imgDir+'img_'+obsId+'_a'+string(a, format='(F5.2)')+'_b'+string(b, format='(F5.2)')+'.sav'
  
  if ~file_exist(fname_out) then begin
   o=obj_new('IDL_Savefile', fname_in)
   o->Restore, 'bestQarr'
   obj_destroy, o 
  
   MakeCompareImages, model, obsImaps, obsSImaps, obsInfo, a, b, bestQarr, 0.1, $
                      modImagesI, modImagesV, rho_arr, chi_arr

   freq=obsInfo.freq
   save, obsId, a, b, bestQarr, freq, $
         modImagesI, modImagesV, rho_arr, chi_arr, $
         filename=fname_out, /compress
        
   obj_destroy, modImagesI
   obj_destroy, modImagesV
  endif else print, 'Images already exist.' 
 endfor 
 
 obj_destroy, model
 obj_destroy, obsImaps
 obj_destroy, obsSImaps
end