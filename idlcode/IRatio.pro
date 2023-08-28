function IRatio, a, b
 qmin=a<b
 qmax=a>b
 return, (qmin le 0) ? 1d100 : qmax/qmin
end 
