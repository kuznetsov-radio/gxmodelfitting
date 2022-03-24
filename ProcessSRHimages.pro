pro DoubleGauss, x, a, f
 x1=a[0]
 A1=a[1]
 s1=a[2]
 x2=a[3]
 A2=a[4]
 s2=a[5]
 
 f=A1*exp(-(x-x1)^2/2/s1^2)+A2*exp(-(x-x2)^2/2/s2^2)
end

pro ProcessSRHimages, index, data, sigma, m
;m: bit 0 (& 1): apply quiet Sun temperature by Landi
;   bit 1 (& 2): redetermine quiet Sun flux using histogram

 forward_function GetQS

 T_old=index.qsun_ref
 T_Landi=GetQS(index.obs_d$freq*1d9)

 if (m and 1) ne 0 then data*=(T_Landi/T_old)
 
 h=histogram(data, locations=x, binsize=10)
 u=where(x lt 6d4)
 h=h[u]
 x=x[u]
 
 a=[0d0, max(h), 1d3, ((m and 1) ne 0) ? T_Landi : T_old, max(h), 1d3]
 q=curvefit(x, h, w, a, function_name='DoubleGauss', /noderivative, /double, itmax=1000)

 T_sky=a[0] 
 s_sky=a[2]
 T_quiet=a[3]
 s_quiet=a[5]

 if (m and 2) ne 0 then begin
  k=(((m and 1) ne 0) ? T_Landi : T_old)/T_quiet
  data*=k
  T_sky*=k
  s_sky*=k
  T_quiet*=k
  s_quiet*=k
 endif
 
 sigma=sqrt((data*s_quiet/T_quiet)^2+s_sky^2)
end