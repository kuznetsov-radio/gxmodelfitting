pro ExtractSubProfile, obsscan, modscan, dx, newscan
 newscan={flux: interpol(obsscan.flux, obsscan.x, modscan.x+dx, /spline), $
          x: modscan.x, $
          freq: obsscan.freq, $
          time: obsscan.time, $
          shiftX: dx}
end