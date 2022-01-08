
pro rbsp_uvw2gse, tname, uvws, newname = newname, inverse = inverse, $
    probe = probe, no_spice_load = no_spice_load

    compile_opt idl2
    on_error, 0
    
    get_data, tname, data = dat
    vec0 = dat.y            ; old vec.
    vec1 = double(vec0)     ; new vec.
    vx0 = vec0[*,0] & vy0 = vec0[*,1] & vz0 = vec0[*,2]
    
    ; get x_mgse, which is also w_sc in gse.
    if n_elements(uvws) eq 0 then begin         ; uvw is set, so use it.
        ; load spice kernel and get uvw.
        ; interpolation because spice runs slow and uvw varies slow.
        if n_elements(probe) eq 0 then prob = 'a' else prob = probe[0]
        if ~keyword_set(no_spice_load) then begin
            rbsp_efw_init
            rbsp_load_spice_kernels
        endif
        t0 = dat.x
        nrec = n_elements(t0)
        uvws = dblarr(nrec,9)
        for i = 0, nrec-1 do begin
            tstr = time_string(t0[i], tformat = 'YYYY-MM-DDThh:mm:ss.ffffff')
            cspice_str2et, tstr, et
            cspice_pxform, 'RBSP'+prob+'_SCIENCE', 'GSE', et, pxform
            uvws[i,*] = pxform[*]   ; the colums are (u)(v)(w) in gse.
        endfor
;        uvws = dblarr(nrec,9)
;        for i = 0, 8 do uvws[*,i] = interpol(uvw0[*,i], tt0, t0)
        
;        ; ca*cb ca*sb*sc-sa*cc  ca*sb*cc+sa*sc
;        ; sa*cb sa*sb*sc+ca*cc  sa*sb*cc-ca*sc
;        ; -sb   cb*sc           cb*cc
;        ; a-theta, b-phi, c-psi
;        
;        ; prepare angle.
;        sb = -uvw[2,0]
;        cb = sqrt(1-sb^2)
;        sc = uvw[2,1]/cb
;        cc = uvw[2,2]/cb
;        ca = uvw[0,0]/cb
;        sa = uvw[1,0]/cb
    endif
    
    ; rotation. break the matrix to do vectorized calc, fast.
    if ~keyword_set(inverse) then begin     ; M2, gse->uvw.
        vx1 = uvws[*,0]*vx0 + uvws[*,1]*vy0 + uvws[*,2]*vz0
        vy1 = uvws[*,3]*vx0 + uvws[*,4]*vy0 + uvws[*,5]*vz0
        vz1 = uvws[*,6]*vx0 + uvws[*,7]*vy0 + uvws[*,8]*vz0
    endif else begin                        ; M1, uvw->gse.
        vx1 = uvws[*,0]*vx0 + uvws[*,3]*vy0 + uvws[*,6]*vz0
        vy1 = uvws[*,1]*vx0 + uvws[*,4]*vy0 + uvws[*,7]*vz0
        vz1 = uvws[*,2]*vx0 + uvws[*,5]*vy0 + uvws[*,8]*vz0
    endelse
    
    if keyword_set(newname) then name = newname else $
    name = keyword_set(inverse)? tname+'_uvw': tname+'_gse'

    store_data, name, data = {x:dat.x, y:[[vx1],[vy1],[vz1]]}
    
end

rbsp_efw_init
tr = time_double(['2013-03-14/07:00','2013-03-14/07:30'])
timespan, tr[0], tr[1]-tr[0], /second
rbsp_load_spice_kernels
probe = 'b'
rbsp_load_spice_state, probe = probe, coord = 'gse', /no_spice_load
rbsp_uvw2gse, 'rbspb_state_pos_gse', newname = 'rbspb_state_pos_uvw', $
    probe = probe, /inverse
rbsp_uvw2gse, 'rbspb_state_pos_uvw', newname = 'rbspb_state_pos_gse2', $
    /no_spice_load, probe = probe

get_data, 'rbspb_state_pos_gse', t0, gse1
get_data, 'rbspb_state_pos_gse2', t0, gse2
store_data, 'rbspb_state_pos_dgse', t0, gse2-gse1
    
vars = ['rbspb_state_pos_*']
options, vars, 'colors', [6,4,2]
tplot, vars

print, max(gse1-gse2)

end