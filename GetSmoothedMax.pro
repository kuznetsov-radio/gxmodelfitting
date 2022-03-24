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
 
 a=[Imax-min(m.data[u]), min(m.data[u]), xc, yc, sx, sy, 0]
 q=curvefit([x[u], y[u]], m.data[u], w, a, function_name='Gauss2DrotC', /noderivative, /double, itmax=1000)
 
 return, a[0]+a[1]
end