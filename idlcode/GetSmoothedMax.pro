pro Gauss2DrotC, u, q, f
 N=n_elements(u)/2
 x=u[0 : N-1]
 y=u[N : n_elements(u)-1]
 
 a=q[0]
 b=q[1]
 xc=q[2]
 yc=q[3]
 sx=q[4]
 sy=q[5]
 theta=q[6]
  
 xx= (x-xc)*cos(theta)+(y-yc)*sin(theta)
 yy=-(x-xc)*sin(theta)+(y-yc)*cos(theta)
 
 f=a*exp(-xx^2/2/sx^2-yy^2/2/sy^2)+b
end

pro Gauss2D, u, q, f
 N=n_elements(u)/2
 x=u[0 : N-1]
 y=u[N : n_elements(u)-1]
 
 a=q[0]
 b=q[1]
 xc=q[2]
 yc=q[3]
 s=q[4]
  
 xx=x-xc
 yy=y-yc
 
 f=a*exp(-xx^2/2/s^2-yy^2/2/s^2)+b
end

function GetSmoothedMax, m, sx, sy
 x=get_map_xp(m, /oned)
 y=get_map_yp(m, /oned)
 
 Imax=0d0
 for i=0, n_elements(x)-1 do for j=0, n_elements(y)-1 do if m.data[i, j] gt Imax then begin
  Imax=m.data[i, j]
  xc=x[i]
  yc=y[j]
 endif
 
 x=get_map_xp(m)
 y=get_map_yp(m)
 
 u=where(((x-xc)^2+(y-yc)^2) lt (sx*sy))
 
 if ((sx/sy) gt 1.01) || ((sy/sx) gt 1.01) then begin
  a=[Imax-min(m.data[u]), min(m.data[u]), xc, yc, sx, sy, 0]
  if n_elements(u) gt n_elements(a) then $
   q=curvefit([x[u], y[u]], m.data[u], w, a, function_name='Gauss2DrotC', /noderivative, /double, itmax=1000, status=st)
 endif else begin
  a=[Imax-min(m.data[u]), min(m.data[u]), xc, yc, sx]
  if n_elements(u) gt n_elements(a) then $
   q=curvefit([x[u], y[u]], m.data[u], w, a, function_name='Gauss2D', /noderivative, /double, itmax=1000, status=st)  
 endelse 
 
 return, (st eq 0) ? a[0]+a[1] : Imax
end