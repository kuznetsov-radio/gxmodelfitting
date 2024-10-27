pro ConvertRATANdata
 fname='C:\MCloud\CoronalMW\AR-SRH\Data\RATAN\RATAN_AR12723_20180930_090339_az0_SCANS__lin_appr.dat'
 
 ;----------------------------------------------------------
 
 rtu_read_scans, fname, scans_R, scans_L, xarc, freqs, pos_angle, index=index
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