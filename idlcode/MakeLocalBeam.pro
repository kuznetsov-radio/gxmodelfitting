pro MakeLocalBeam, obsInfo, j, dx, dy, beam
 beam0=obsInfo.psf
 beam0=beam0[j]
 bsize0=size(beam0, /dimensions)
  
 Nx=round(bsize0[0]*obsInfo.psf_dx[j]/dx)
 Ny=round(bsize0[1]*obsInfo.psf_dy[j]/dy)
 beam=congrid(beam0, Nx, Ny, /center, cubic=-0.5)
 
 beam/=total(beam)
end