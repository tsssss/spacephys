;+
; Test 1D model along the RBSP-A field line.
;-

; z in Re.
function model_ne, z

    n0 = 6e4    ; cc.
    n1 = 1.34e7 ; cc.
    z0 = 318.   ; km.
    h = 383.    ; km.

    re = constant('re')
    zz = z*re
    t = (z0-zz)/h
    n = n0*exp(t)+n1*zz^(-1.55)
    return, n

end

function model_te, z

    re = constant('re')
    zz = z*re
    t1 = 0.0135 ; eV.
    h0 = 2000d  ; km.
    h0 = 3000d  ; km.   ; 200 eV at s/c.
    h0 = 3500d  ; km.   ; 50 eV at s/c.
    t0 = 1d     ; eV.
    t = zz/h0
    tl = t1*exp(t)+t0

    tps = 2e3   ; eV.
    dz = 0.3    ; Re.
    zps = 3.75  ; Re.
    w = 0.5*(1-tanh((z-zps)/dz))
    te = tl*w+tps*(1-w)
    te = tl

    return, te

end


function model_ro, z, target=target

    re = constant('re')
    zz = z*re
    z1 = 2370d  ; km.
    h1 = 1800d  ; km.
    tt = (zz-z1)/h1
    ro = 0.5*(1-tanh(tt))
    ro = (1+target)*0.5-tanh(tt)*(1-target)*0.5
    return, ro
    
end


    event_info = _2015_0218_02_load_data()
    test = 1


    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = time_double('2015-02-18/'+['02:05','02:15'])
    ndim = 3
    rad = constant('rad')
    deg = constant('deg')

    ; Get the s/c location.
    the_time = mean(time_range)
    r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
    rx = r_gsm[0]
    ry = r_gsm[1]
    rz = r_gsm[2]
    dir = -1    ; map to N-hem.


;---Get the locations of the field line.
    tilt = geopack_recalc(the_time)

    geopack_trace, rx,ry,rz, dir, 0, tilt=tilt, fx,fy,fz, igrf=0, fline=fline_gsm, $
        t89=0, t96=0, t01=0
    index = where(snorm(fline_gsm) gt 1)
    fline_gsm = fline_gsm[reverse(index),*]
    fline = cotran(fline_gsm, the_time, 'gsm2sm')
    nrec = n_elements(fline[*,0])
    z = snorm(fline)-1
    
    ; Distance along field line.
    ss = fltarr(nrec)
    ss[0] = z[0]
    for ii=1,nrec-1 do begin
        ds = snorm(fline[ii,*]-fline[ii-1,*])
        ss[ii] += ss[ii-1]+ds
    endfor

    ; Density and temperature.
    n_e = model_ne(z)
    t_e = model_te(z)
    
    the_n_e = event_info['num_dens']
    n_e = n_e+(the_n_e-n_e[-1])

    ; B.
    b_gsm = fltarr(nrec,ndim)
    for ii=0,nrec-1 do begin
        geopack_dip, fline_gsm[ii,0],fline_gsm[ii,1],fline_gsm[ii,2], bx,by,bz
        b_gsm[ii,*] = [bx,by,bz]
    endfor
    bmag = snorm(b_gsm)
    the_bmag = event_info['bmag']
    the_ratio = the_bmag/bmag[-1]
    norm_ratio = 1-tanh(5-z)
    norm_ratio = 1-norm_ratio/norm_ratio[-1]*(1-the_ratio)
    bmag = bmag*norm_ratio
    
    
    ; O ratio.
    r_o = model_ro(z, target=event_info['o_num_dens_ratio'])
    
    ; v_A.
    mu0 = 4*!dpi*1e-7
    m_p = 1.67e-27      ; kg.
    mass_ratio = 16*r_o+1*(1-r_o)
    v_a = bmag/sqrt(mu0*n_e*m_p*mass_ratio)*1e-9/1e3*1e-3
    
    ; v_eth.
    m_e = 0.91e-30      ; kg.
    v_eth = sqrt(t_e*1.6e-19/m_e)*1e-3


    sgopen, 0, xsize=4, ysize=6
    nvar = 6
    margins = [12,5,2,2]
    poss = sgcalcpos(nvar, margins=margins)
    
    !p.charsize = 1.5
    
    tpos = poss[*,0]
    plot, z+1, bmag, $
        ylog=1, ytitle='|B| (nT)', $
        xstyle=1, xtickformat='(A1)', $
        position=tpos, noerase=1

    tpos = poss[*,1]
    plot, z+1, n_e, $
        ylog=1, ytitle='n!De!N (cm!U-3!N)', $
        xstyle=1, xtickformat='(A1)', $
        position=tpos, noerase=1

    tpos = poss[*,2]
    plot, z+1, t_e, $
        ylog=1, ytitle='T!De!N (eV)', $
        xstyle=1, xtickformat='(A1)', $
        position=tpos, noerase=1
        
    tpos = poss[*,3]
    plot, z+1, v_eth, $
        ylog=1, ytitle='v!DTe!N (km/s)', $
        xstyle=1, xtickformat='(A1)', $
        position=tpos, noerase=1
    
    tpos = poss[*,4]
    plot, z+1, v_a, $
        ylog=1, ytitle='v!DA!N (km/s)', $
        xstyle=1, xtickformat='(A1)', $
        position=tpos, noerase=1
        
    tpos = poss[*,5]
    plot, z+1, v_a/v_eth, $
        ylog=1, ytitle='v!DA!N/v!Dth,e!N', $
        xstyle=1, xtitle='R (Re)', $
        position=tpos, noerase=1
    oplot, minmax(z)+1, [1,1], linestyle=1
    
end
