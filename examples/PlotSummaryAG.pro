pro PlotSummaryAG
 restore, 'C:\MCloud\CoronalMW\AR-SRH\Data\Images\SRH\mid\DDM_SS\Summary_eta_thr0.100_M_20230528_0355_05.80GHz.sav'
 threshold_metric=2.0
 
 ;------------------------------------------------------------------
 
 a_arr=alist
 b_arr=blist
 freq=freqlist
 
 Na=n_elements(a_arr)
 Nb=n_elements(b_arr)
 Nfreq=n_elements(freq)
 
 set_plot, 'ps'
 
 for k=0, Nfreq-1 do begin
  dx=a_arr[1]-a_arr[0]
  dy=b_arr[1]-b_arr[0]
  m0=make_map(eta[*, *, k], xc=(min(a_arr)+max(a_arr))/2, yc=(min(b_arr)+max(b_arr))/2, dx=dx, dy=dy)
              
  m=m0
  u=where(m.data lt 0, g)
  if g gt 0 then m.data[u]=!values.d_NaN
  u=where(m.data gt (min(m.data, /nan)*threshold_metric), g)
  if g gt 0 then m.data[u]=!values.d_NaN
  
  x=get_map_xp(m)
  y=get_map_yp(m)
  u=where(finite(m.data))
  xmin=min(x[u])-dx
  xmax=max(x[u])+dx
  ymin=min(y[u])-dy
  ymax=max(y[u])+dy
  m=get_sub_map(m, xrange=[xmin, xmax], yrange=[ymin, ymax])
  m0=get_sub_map(m0, xrange=[xmin, xmax], yrange=[ymin, ymax])
  
  xsize=8.0
  yw=7.0*(ymax-ymin+dy)/(xmax-xmin+dx)
  ysize=yw+1.0+0.87
  position=[0.17, 0.87/ysize, 0.97, (yw+0.87)/ysize]
              
  device, file='FigEtaAB_'+string(k, format='(I02)')+'_'+string(freq[k], format='(F05.2)')+'GHz.eps', $
          xsize=xsize, ysize=ysize, font_size=10, bits_per_pixel=8, /encapsulated, /color
  loadct, 13, /silent
  plot_map, m, position=position, /isotropic, cbar=1, $
            title='!7g!U!32!N  '+obsId+'  '+string(freq[k], format='(F5.2)')+' GHz', $
            xtitle='!18a!3', ytitle='!18b!3'
            
  loadct, 0, /silent
  x=get_map_xp(m0)
  y=get_map_yp(m0)          
  u=where(~finite(m0.data), g)
  if g gt 0 then for i=0, g-1 do begin
   plots, x[u[i]]+dx*[-0.5, 0.5], y[u[i]]+dy*[-0.5, 0.5], /data, color=128 
   plots, x[u[i]]+dx*[-0.5, 0.5], y[u[i]]-dy*[-0.5, 0.5], /data, color=128
  endfor  
  
  u=where((m0.data lt 0) or (m0.data gt (min(m.data, /nan)*2)), g)
  if g gt 0 then plots, x[u], y[u], /data, psym=3, color=128
  
  cm=min(m.data, g, /nan)
  plots, x[g], y[g], /data, psym=2, color=255  
  
  device, /close
 endfor 
end