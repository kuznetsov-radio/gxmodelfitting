pro DoubleGauss, x, a, f
 x1=a[0]
 A1=a[1]
 s1=a[2]
 x2=a[3]
 A2=a[4]
 s2=a[5]
 
 f=A1*exp(-(x-x1)^2/2/s1^2)+A2*exp(-(x-x2)^2/2/s2^2)
end

function map2sdev, map
 h=histogram(map.data, locations=x, binsize=10)
 u=where(x lt 6d4)
 h=h[u]
 x=x[u]
 
 a=[0d0, max(h), 1d3, 20d3, max(h), 1d3]
 q=curvefit(x, h, w, a, function_name='DoubleGauss', /noderivative, /double, itmax=1000)

 T_sky=a[0] 
 s_sky=a[2]
 T_quiet=a[3]
 s_quiet=a[5]

 sdev=map 
 sdev.data=sqrt((map.data*s_quiet/T_quiet)^2+s_sky^2)
 
 return, sdev
end