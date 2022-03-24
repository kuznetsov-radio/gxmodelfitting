pro ComputeGXmaps, model, freqlist, Q, a, b, Imap, Vmap, dem_on=dem_on
 ebtel_path=gx_findfile(keyword_set(dem_on) ? 'ebtel_runsGG_impulsive_ddm.sav' : 'ebtel.sav', folder='')
 renderer=gx_findfile('grff_dem_transferM.pro', folder='')
 
 q0_formula='q[0]'
 q_formula='q[0]*(B/q[1])^(q[3])/(L/q[2])^(q[4])'
 q_parms=[Q, 100.0, 1e9, a, b] ;q0, B0, L0, a, b

 dk=keyword_set(dem_on) ? 0 : 1 
 omap=gx_mwrender_ebtel(model, renderer, ebtel_path=ebtel_path, q_parms=q_parms, q_formula=q_formula, q0_formula=q0_formula, $
                        f_min=0.0, n_freq=n_elements(freqlist), freqlist=freqlist, DEM_key=dk, DDM_key=dk, gxcube=gxcube)
 
 Imap=obj_new('map')
 Vmap=obj_new('map')
 for i=0, n_elements(freqlist)-1 do begin
  Imap->setmap, i, omap->get(i*2, /map)
  Vmap->setmap, i, omap->get(i*2+1, /map)
 endfor
 obj_destroy, omap
end