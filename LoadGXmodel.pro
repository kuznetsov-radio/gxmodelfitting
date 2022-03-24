function LoadGXmodel, filename, xc, yc, xfov, yfov, nx, ny
 model=gx_read(filename) 
 fovdata=model->SetFOV(xc=xc, yc=yc, xfov=xfov, yfov=yfov, nx=nx, ny=ny, /compute_grid)
 return, model
end