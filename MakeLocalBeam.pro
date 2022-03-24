pro MakeLocalBeam, obsInfo, j, dx, dy, beam
 if obsInfo.id eq 'low' then begin
  Nx=round(180.0/dx)
  Ny=round(180.0/dy)
  MakeSRHbeam, obsInfo.sx[j], obsInfo.sy[j], obsInfo.rho[j], Nx, Ny, dx, dy, x, y, beam 
 endif else if obsInfo.id eq 'mid' then begin
  beam0=obsInfo.beam
  beam0=beam0[j]
  bsize0=size(beam0, /dimensions)
  
  Nx=round(bsize0[0]*obsInfo.dx[j]/dx)
  Ny=round(bsize0[1]*obsInfo.dy[j]/dy)
  beam=congrid(beam0, Nx, Ny, /center, cubic=-0.5)
 endif 
 
 beam/=total(beam)
end