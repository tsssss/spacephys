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

function model_te, z, h0=h0

    re = constant('re')
    zz = z*re
    t1 = 0.0135 ; eV.

;    h0 = 2000d  ; km.
;    h0 = 3000d  ; km.   ; 200 eV at s/c.
;    h0 = 3500d  ; km.   ; 50 eV at s/c.

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




    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    test = 1


    probe = event_info['probe']
    prefix = event_info['prefix']
    time_range = event_info['time_range']
    ndim = 3
    rad = constant('rad')
    deg = constant('deg')

    ; Get the s/c location.
    the_time = mean(time_range)
    r_gsm = get_var_data(prefix+'r_gsm', at=the_time)
    xp = r_gsm[0]
    yp = r_gsm[1]
    zp = r_gsm[2]
    trace_dir = -1    ; map to N-hem.
    par = 2
    refine = 1
    h0 = 100d
    r0 = h0/constant('re')+1

;---Get the locations of the field line.
    model_setting = event_info['model_setting']
    model = model_setting['external_model']
    igrf = model_setting['igrf']
    ps = geopack_recalc(the_time)
    
    tmp = geopack_resolve_model(model)
    t89 = tmp.t89
    t96 = tmp.t96
    t01 = tmp.t01
    ts04 = tmp.ts04
    storm = tmp.storm

    geopack_trace, xp,yp,zp, trace_dir, par, tilt=ps, $
        xf,yf,zf, r0=r0, refine=refine, ionosphere=1, $
        fline=fline_gsm, $
        t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf
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
    t_p = model_te(z, h0=2500d)*5
    

    plasma_param = event_info['plasma_param']
    the_n_e = plasma_param['num_dens']
    n_e = n_e+(the_n_e-n_e[-1])

    ; B.
    b_gsm = fltarr(nrec,ndim)
    for ii=0,nrec-1 do begin
        geopack_dip, fline_gsm[ii,0],fline_gsm[ii,1],fline_gsm[ii,2], bx,by,bz
        b_gsm[ii,*] = [bx,by,bz]
    endfor
    bmag = snorm(b_gsm)
    the_bmag = plasma_param['bmag']
    the_ratio = the_bmag/bmag[-1]
    norm_ratio = 1-tanh(5-z)
    norm_ratio = 1-norm_ratio/norm_ratio[-1]*(1-the_ratio)
    bmag = bmag*norm_ratio


    ; O ratio.
    r_o = model_ro(z, target=plasma_param['o_num_dens_ratio'])

    ; v_A.
    mu0 = 4*!dpi*1e-7
    m_p = 1.67e-27      ; kg.
    mass_ratio = 16*r_o+1*(1-r_o)
    v_a = bmag/sqrt(mu0*n_e*m_p*mass_ratio)*1e-9/1e3*1e-3

    ; v_pth.
    m_p = 1.67e-27      ; kg.
;    m_p = plasma_param['avg_ion_mass']*m_p
    m_e2 = m_p*16
    v_pth = sqrt(t_p*1.6e-19/m_p)*1e-3
    v_oth = sqrt(t_p*1.6e-19/m_e2)*1e-3


    plot_file = join_path([srootdir(),'fig_2013_0501_1d_model.pdf'])
    if keyword_set(test) then plot_file = 0
    if keyword_set(test) then magn = 2 else magn = 1
    sgopen, plot_file, xsize=5, ysize=5, xchsz=xchsz, ychsz=ychsz, magn=magn
    nvar = 6
    lmarg = 1
    margins = [10,4,2,1]
    poss = sgcalcpos(nvar, margins=margins)


    color2 = sgcolor('gray')
    linestyle2 = 2
    xtitle = '|R| (Re)'
    !p.charsize = 1
    xticklen = -0.03
    yticklen = -0.01

    tpos = poss[*,0]
    ytitle = '(nT)'
    ytickv_log = make_bins([2,5],1)
    ytickv = 10^ytickv_log
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = '10!U'+string(ytickv_log,format='(I0)')
    index = where(ytickv_log eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(ytickv_log eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    index = where(ytickv_log eq -1, count)
    if count ne 0 then ytickn[index] = '0.1'
    index = where(ytickv_log eq 2, count)
    if count ne 0 then ytickn[index] = '100'
    index = where(ytickv_log eq -2, count)
    if count ne 0 then ytickn[index] = '0.01'
    ytickn[1:*:2] = ' '
    plot, z+1, bmag, $
        ystyle=1, ylog=1, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, ytitle=ytitle, yticklen=yticklen, $
        xstyle=1, xtickformat='(A1)', xticklen=xticklen, $
        position=tpos, noerase=1
    tx = xchsz*lmarg
    ty = tpos[3]-ychsz*0.7
    msg = 'a) |B|'
    xyouts, tx,ty,msg, normal=1

;    tx = tpos[0]+xchsz*1
;    ty = tpos[3]+ychsz*0.5
;    xyouts, tx,ty,normal=1, 'Plasma sheet'
;    tx = tpos[0]+xchsz*15
;    xyouts, tx,ty,normal=1, 'Plasmasphere', color=color2



    tpos = poss[*,1]
    ytitle = '(eV)'
    data = t_p
    yrange_log = alog10(minmax(data))
    yrange_log = [floor(yrange_log[0]),ceil(yrange_log[1])]
    ytickv_log = make_bins(yrange_log,1)
    ytickv = 10^ytickv_log
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = '10!U'+string(ytickv_log,format='(I0)')
    index = where(ytickv_log eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(ytickv_log eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    index = where(ytickv_log eq -1, count)
    if count ne 0 then ytickn[index] = '0.1'
    index = where(ytickv_log eq 2, count)
    if count ne 0 then ytickn[index] = '100'
    index = where(ytickv_log eq -2, count)
    if count ne 0 then ytickn[index] = '0.01'
    ytickn[1:*:2] = ' '
    plot, z+1, t_p, $
        ystyle=1, ylog=1, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, ytitle=ytitle, yticklen=yticklen, $
        xstyle=1, xtickformat='(A1)', xticklen=xticklen, $
        position=tpos, noerase=1
;    oplot, z+1, t_e2, color=color2, linestyle=linestyle2
    tx = xchsz*lmarg
    ty = tpos[3]-ychsz*0.7
    msg = 'b) T!DH+&O+'
    xyouts, tx,ty,msg, normal=1


    tpos = poss[*,2]
    ytitle = '(km/s)'
    v_oth = v_pth/sqrt(16)
    yrange_log = alog10(minmax([v_pth,v_oth]))
    yrange_log = [floor(yrange_log[0]),ceil(yrange_log[1])]
    ytickv_log = make_bins(yrange_log,1)
    ytickv = 10^ytickv_log
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = '10!U'+string(ytickv_log,format='(I0)')
    index = where(ytickv_log eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(ytickv_log eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    index = where(ytickv_log eq -1, count)
    if count ne 0 then ytickn[index] = '0.1'
    index = where(ytickv_log eq 2, count)
    if count ne 0 then ytickn[index] = '100'
    index = where(ytickv_log eq -2, count)
    if count ne 0 then ytickn[index] = '0.01'
    ytickn[1:*:2] = ' '
    plot, z+1, v_pth, $
        ystyle=1, ylog=1, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, ytitle=ytitle, yticklen=yticklen, $
        xstyle=1, xtickformat='(A1)', xticklen=xticklen, $
        position=tpos, noerase=1
    oplot, z+1, v_oth, color=color2, linestyle=linestyle2
    tx = xchsz*lmarg
    ty = tpos[3]-ychsz*0.7
    msg = 'c) v!DT'
    xyouts, tx,ty,msg, normal=1
    

    tpos = poss[*,3]
    ytitle = '(Hz)'
    q = 1.6e-19
    f_gp = q*bmag/m_p*1e-9
    f_go = f_g/16
    yrange_log = alog10(minmax([f_gp,f_go]))
    yrange_log = [floor(yrange_log[0]),ceil(yrange_log[1])]
    ytickv_log = make_bins(yrange_log,1)
    ytickv = 10^ytickv_log
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = '10!U'+string(ytickv_log,format='(I0)')
    index = where(ytickv_log eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(ytickv_log eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    index = where(ytickv_log eq -1, count)
    if count ne 0 then ytickn[index] = '0.1'
    index = where(ytickv_log eq 2, count)
    if count ne 0 then ytickn[index] = '100'
    index = where(ytickv_log eq -2, count)
    if count ne 0 then ytickn[index] = '0.01'
    ytickn[1:*:2] = ' '
    plot, z+1, f_gp, $
        ystyle=1, ylog=1, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, ytitle=ytitle, yticklen=yticklen, $
        xstyle=1, xtickformat='(A1)', xticklen=xticklen, $
        position=tpos, noerase=1
    oplot, z+1, f_go, color=color2, linestyle=linestyle2
    tx = xchsz*lmarg
    ty = tpos[3]-ychsz*0.7
    msg = 'd) f!Dg'
    xyouts, tx,ty,msg, normal=1


    tpos = poss[*,4]
    ytitle = '(km)'
    r_gp = v_pth/f_gp
    r_go = v_oth/f_go
    yrange_log = alog10(minmax([r_gp,r_go]))
    yrange_log = [floor(yrange_log[0]),ceil(yrange_log[1])]
    ytickv_log = make_bins(yrange_log,1)
    ytickv = 10^ytickv_log
    ytickv = 10^ytickv_log
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = '10!U'+string(ytickv_log,format='(I0)')
    index = where(ytickv_log eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(ytickv_log eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    index = where(ytickv_log eq -1, count)
    if count ne 0 then ytickn[index] = '0.1'
    index = where(ytickv_log eq 2, count)
    if count ne 0 then ytickn[index] = '100'
    index = where(ytickv_log eq -2, count)
    if count ne 0 then ytickn[index] = '0.01'
    ytickn[0:*:2] = ' '
    plot, z+1, r_gp, $
        ystyle=1, ylog=1, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, ytitle=ytitle, yticklen=yticklen, $
        xstyle=1, xtickformat='(A1)', xticklen=xticklen, $
        position=tpos, noerase=1
    oplot, z+1, r_go, color=color2, linestyle=linestyle2
    tx = xchsz*lmarg
    ty = tpos[3]-ychsz*0.7
    msg = 'e) R!Dg'
    xyouts, tx,ty,msg, normal=1
    
    
    
    tpos = poss[*,5]
    lambda_perp0 = 5e3
    lambda_perp = lambda_perp0/bmag*min(bmag)
    r = z+1
    ;lambda_perp = lambda_perp0/r*max(r)
    ytitle = '(#)'
    yrange_log = alog10(minmax([r_gp/lambda_perp,r_go/lambda_perp]))
    yrange_log = [floor(yrange_log[0]),ceil(yrange_log[1])]
    ytickv_log = make_bins(yrange_log,1)
    ytickv = 10^ytickv_log
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytickn = '10!U'+string(ytickv_log,format='(I0)')
    index = where(ytickv_log eq 0, count)
    if count ne 0 then ytickn[index] = '1'
    index = where(ytickv_log eq 1, count)
    if count ne 0 then ytickn[index] = '10'
    index = where(ytickv_log eq -1, count)
    if count ne 0 then ytickn[index] = '0.1'
    index = where(ytickv_log eq 2, count)
    if count ne 0 then ytickn[index] = '100'
    index = where(ytickv_log eq -2, count)
    if count ne 0 then ytickn[index] = '0.01'
    ytickn[0:*:2] = ' '
    
    plot, z+1, r_gp/lambda_perp, $
        ystyle=1, ylog=1, ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, ytitle=ytitle, yticklen=yticklen, $
        xstyle=1, xtickformat='', xticklen=xticklen, xtitle=xtitle, $
        position=tpos, noerase=1
    oplot, z+1, r_go/lambda_perp, color=color2, linestyle=linestyle2
    tx = xchsz*lmarg
    ty = tpos[3]-ychsz*0.7
    msg = 'f) R!Dg!N/'+tex2str('lambda')+'!D'+tex2str('perp')
    xyouts, tx,ty,msg, normal=1
    
    tpos = poss[*,-1]
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'H+'
    xyouts, tx,ty,normal=1, msg
    plots, tx+xchsz*3+xchsz*1*[0,1],ty+ychsz*0.3+[0,0], normal=1
    ty = tpos[3]-ychsz*2
    msg = 'O+'
    xyouts, tx,ty,normal=1, msg
    plots, tx+xchsz*3+xchsz*1*[0,1],ty+ychsz*0.3+[0,0], normal=1, linestyle=linestyle2, color=color2

    if keyword_set(test) then stop
    sgclose

end