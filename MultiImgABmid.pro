pro MultiImgABmid
 imgDir='C:\MCloud\CoronalMW\AR-SRH\Data\Images\mid\'
 obsFile='C:\MCloud\CoronalMW\AR-SRH\Data\SRH\mid\srh_median.sav'
 obsId='20220109' 
 model=LoadGXmodel('C:\MCloud\CoronalMW\AR-SRH\Data\Models\model20220109.gxm', $
                   55.0, -450.0, 300.0, 300.0, 150, 150)

 LoadSRHmid, obsFile, obsImaps, obsSImaps, obsInfo, 0
 
 a_arr=[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]
 b_arr=[1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]
 
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

   fname_out=imgDir+'img_'+obsId+'_a'+string(a, format='(F5.2)')+'_b'+string(b, format='(F5.2)')+'.sav'
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
