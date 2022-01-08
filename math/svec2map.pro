
function svec2map, v0, p0, scl, log = log
    on_error, 0
    
;    if n_params() ne 2 then message, 'need vec and origin data ...'
    
    ; x: azimuthal, y: radial.
    nrec = n_elements(v0)/2
    v1 = sqrt(v0[0]^2+v0[1]^2)
    
    ; p1: (lon,lat).
    p1 = reform(p0)
    
    ; xy1: (x,y) for p1.
    xy1 = map_proj_forward(p1[0],p1[1])     ; (lon,lat) -> (x,y).
    t = atan(-xy1[0],-xy1[1])
    cost = cos(t)
    sint = sin(t)
    
    v = sqrt(v0[0]^2+v0[1]^2)
    if keyword_set(log) then v = alog(v)
    a = atan(v0[1],v0[0])
    sina = sin(a)
    cosa = cos(a)

    ; xy2: (x2,y2) for p2.
    x2 = xy1[0]+(sina*sint+cosa*cost)*scl*v
    y2 = xy1[1]+(sina*cost-cosa*sint)*scl*v
    
    ; p2: (lon,lat).
    p2 = map_proj_inverse(x2,y2)               ; (x,y) -> (lon,lat).
    
    return, p2
end