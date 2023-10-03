pro ExtractSubmap, obsmap, modmap, dx, dy, id, newmap
 x1=get_map_xp(obsmap, /oned)
 y1=get_map_yp(obsmap, /oned)
 x2=get_map_xp(modmap)+dx
 y2=get_map_yp(modmap)+dy
 
 xc=(x2-x1[0])/(x1[1]-x1[0])
 yc=(y2-y1[0])/(y1[1]-y1[0])
 
 a=interpolate(obsmap.data, xc, yc, missing=0)
 
 newmap=modmap
 newmap.data=a
 newmap.id=id
end 