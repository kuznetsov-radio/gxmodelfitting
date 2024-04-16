pro PlotSummaryFG
 restore, 'C:\MCloud\CoronalMW\AR-SRH\Data\Images\srh\mid\DDM_SS\Summary_eta_thr0.100_M_20230528_0355.sav'
 
 ;------------------------------------------------------------------
 
 a_arr=alist
 b_arr=blist
 freq=freqlist
 
 Na=n_elements(a_arr)
 Nb=n_elements(b_arr)
 Nfreq=n_elements(freq)
 
 dx=a_arr[1]-a_arr[0]
 xmin=min(a_arr)-dx
 xmax=max(a_arr)+dx
 dy=b_arr[1]-b_arr[0]
 ymin=min(b_arr)-dy
 ymax=max(b_arr)+dy
 
 xsize=8.0
 yw=7.0*(ymax-ymin+dy)/(xmax-xmin+dx)
 ysize=yw+1.0+0.87
 position=[0.17, 0.87/ysize, 0.97, (yw+0.87)/ysize]
 
 set_plot, 'ps'
 
 for k=0, Nfreq-1 do begin
  m=make_map(bestQ[*, *, k], $
             xc=(min(a_arr)+max(a_arr))/2, yc=(min(b_arr)+max(b_arr))/2, $
             dx=a_arr[1]-a_arr[0], dy=b_arr[1]-b_arr[0]) 
  device, file='FigBestQ_'+string(k, format='(I02)')+'_'+string(freq[k], format='(F05.2)')+'GHz.eps', $
          xsize=xsize, ysize=ysize, font_size=10, bits_per_pixel=8, /encapsulated, /color
  loadct, 13, /silent     
  plot_map, m, position=position, /isotropic, cbar=1, /log, $
            title='(!8Q!3!D0!N)!Dopt!N  '+obsId+'  '+string(freq[k], format='(F5.2)')+' GHz   ', $
            xtitle='!18a!3', ytitle='!18b!3'
  loadct, 0, /silent
  x=get_map_xp(m)
  y=get_map_yp(m)            
  u=where(~finite(m.data), g)
  if g gt 0 then plots, x[u], y[u], /data, psym=7, color=128
  device, /close      
  
  m=make_map(eta[*, *, k], $
             xc=(min(a_arr)+max(a_arr))/2, yc=(min(b_arr)+max(b_arr))/2, $
             dx=a_arr[1]-a_arr[0], dy=b_arr[1]-b_arr[0])
  device, file='FigEtaAB_'+string(k, format='(I02)')+'_'+string(freq[k], format='(F05.2)')+'GHz.eps', $
          xsize=xsize, ysize=ysize, font_size=10, bits_per_pixel=8, /encapsulated, /color
  loadct, 13, /silent
  plot_map, m, position=position, /isotropic, cbar=1, /log, $
            title='!7g!U!32!N  '+obsId+'  '+string(freq[k], format='(F5.2)')+' GHz   ', $
            xtitle='!18a!3', ytitle='!18b!3'
  loadct, 0, /silent
  x=get_map_xp(m)
  y=get_map_yp(m)          
  u=where(~finite(m.data), g)
  if g gt 0 then plots, x[u], y[u], /data, psym=7, color=128
  cm=min(m.data, g, /nan)
  plots, x[g], y[g], /data, psym=2, color=255  
  x=get_map_xp(m, /oned)
  y=get_map_yp(m, /oned)     
  contour, m.data, x, y, levels=min(m.data, /nan)*2, c_color=255, /overplot
  device, /close
 endfor 
end