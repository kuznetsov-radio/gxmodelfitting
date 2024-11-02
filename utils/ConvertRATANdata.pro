pro ConvertRATANdata, InFile
;This program converts the RATAN data (.dat) into .sav reference files accepted by the CHMP routines.
;Input parameter: 
; InFile - name of the RATAN data file (.dat).
;Output: the program creates .sav files in the current directory, one file per frequency channel.

 rtu_read_scans, InFile, scans_R, scans_L, xarc, freqs, pos_angle, index=index
 if index.shift_hmi ne 0 then xarc+=index.shift_hmi
 xarc-=index.xc
 scans=scans_L+scans_R
 
 instrument='RATAN'
 
 for i=0, n_elements(freqs)-1 do begin
  ref={flux: reform(scans[*, i]), $
       x: xarc, $
       freq: freqs[i]/1d9, $
       time: index.date_obs, $
       rot: pos_angle}
       
  save, instrument, ref, filename='RATAN_'+string(ref.freq, format='(F05.2)')+'GHz.sav'
 endfor
end