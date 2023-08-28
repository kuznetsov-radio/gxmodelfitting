pro FixLocalBeam, beam, Nx, Ny
 s=size(beam, /dimensions)
  
 if s[0] ge Nx then begin
  Nx1=(Nx/2-1)*2+1
  k=(s[0]-Nx1)/2
  beam=beam[k : k+Nx1-1, *] 
 endif
 
 if s[1] ge Ny then begin
  Ny1=(Ny/2-1)*2+1
  k=(s[1]-Ny1)/2
  beam=beam[*, k : k+Ny1-1]
 endif
end
