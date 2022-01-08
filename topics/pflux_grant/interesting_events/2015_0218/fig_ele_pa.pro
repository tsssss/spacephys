;+
; Electron pitch angle distribution and related plasma variables.
;-


    event_info = _2015_0218_02_load_data()


;---Settings.
    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']
test = 0

    pa_times = time_double('2015-02-18/'+[$
        '02:07:20','02:08:28','02:09:37','02:10:45'])


    perp = '!9'+string(94b)+'!X'
    labfac = ['||',perp+',West',perp+',Out']
    rgb = constant('rgb')
    xyz = ['x','y','z']
    label_size = 0.8

    ticklen = -0.02
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 4
    tplot_options, 'labflag', -1
    tplot_options, 'ticklen', ticklen
    tplot_options, 'zticklen', -0.1
    tplot_options, 'thick', 1
    tplot_options, 'charsize', 1
    tplot_options, 'xcharsize', 1
    tplot_options, 'ycharsize', 1


;---Calc fig size.
    npa = n_elements(pa_times)
    pa_margins = [3,1,1,3]
    pa_margins = [6,2,8,2]
    poss = panel_pos(0, pansize=[1,1], margins=pa_margins, nxpan=npa, fig_size=fig_size, xpad=1)
    plot_vars = prefix+['e_en_spec','combo_density','combo_temp','combo_pressure','beta']
    dy = 0.7
    plot_ypans = [1,dy,dy,1,1]
    nvar = n_elements(plot_vars)
    fig_labels = (letters(nvar+1))[1:*]+'. '+['Ele','n','T','P',tex2str('beta')]
    the_poss = panel_pos(0, pansize=fig_size, margins=[0,0,0,0], ypans=[1,nvar*0.6], fig_size=fig_size)

    plot_file = join_path([srootdir(),'fig_ele_pa.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_size[0], ysize=fig_size[1], xchsz=xchsz, ychsz=ychsz



;---Plot PA 2d, adopted from plot_hope_l3_pitch2d.pro
    pa_poss = sgcalcpos(1,npa, region=the_poss[*,0], margins=pa_margins, xpad=1)
    log = 1
    unit = 'energy'
    hopel3 = sread_rbsp_hope_l3(time_range, probes=probe)
    type = 'electron'
    type0 = 'FEDU'
    zrange = [4.,10]
    npxl = 5
    rad = constant('rad')
    deg = constant('deg')
    mass0 = 1d/1836
    mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.
    ticklen = -0.02

    epidx = (type eq 'electron')? 'EPOCH_ELE': 'EPOCH_ION'
    epidx = where(tag_names(hopel3) eq epidx)
    enidx = (type eq 'electron')? 'HOPE_ENERGY_ELE': 'HOPE_ENERGY_ION'
    enidx = where(tag_names(hopel3) eq enidx)

    foreach time, pa_times, time_id do begin
        tmp = min(hopel3.(epidx)-stoepoch(time,'unix'),rec,/absolute)
        ut0 = sfmepoch(hopel3.(epidx)[rec],'unix')  ; time in the middle of a spin.
        i0 = where(tag_names(hopel3) eq type0)
        dtidx = (type eq 'electron')? 'EPOCH_ELE_DELTA': 'EPOCH_ION_DELTA'
        dtidx = where(tag_names(hopel3) eq dtidx)
        type_dt = hopel3.(dtidx)[rec]*1e-3

        ; level 3.
        datl3 = reform((hopel3.(i0))[rec,*,*])            ; in [nen,npa].
        enl3s = reform(hopel3.(enidx)[rec,*])
        pal3s = reform(hopel3.pitch_angle)
        npal3 = n_elements(pal3s)
        idx = where(datl3 eq -1e31, cnt)
        if cnt ne 0 then datl3[idx] = !values.d_nan
        idx = where(datl3 eq 0, cnt)
        if cnt ne 0 then datl3[idx] = !values.d_nan

    ;    ; remove duplicated energy bins.
        idx = uniq(enl3s,sort(enl3s))
        enl3s = enl3s[idx]
        nenl3 = n_elements(enl3s)
        datl3 = datl3[idx,*]

        ; the data for polar contour.
        tdat = transpose([[datl3],[datl3]])     ; in [2*npa,nen].
        tang = [pal3s,360-pal3s]
        case unit of
            'energy': begin
                tdis = enl3s
                xtitl = 'E (eV)'
                end
            'velocity': begin
                tdis = sqrt(2*enl3s/mass0)*1e-3
                xtitl = 'V (km/s)'
                end
        endcase
        if keyword_set(log) then begin
            tdis = alog10(tdis)
            xtitl = 'Log!D10!N '+xtitl
        endif

        tang = tang # ((bytarr(nenl3)+1)+smkarthm(0,0.001,nenl3,'n'))
        tdis = tdis ## (bytarr(2*npal3)+1)

        tang = tang*rad

        ; remove nan.
        idx = where(finite(tdat,/nan))
        tdat[idx] = 0

        idx = where(tdat ne 0)
        min0 = min(tdat[idx],/nan)
        max0 = max(tdat[idx],/nan)
        nztick = 10
        if n_elements(zrange) eq 0 then zrange = [floor(alog10(min0)),ceil(alog10(max0))-2]

        tpos = pa_poss[*,time_id]
        sgindexcolor, 43, file = 'ct2'
        no_cb = time_id ne npa-1
        sgdistr2d, tdat, tang, tdis, position=tpos, zrange=zrange, $
            ncolor=10, $
            no_cb=no_cb, no_data_point=1, no_axis=1

        xrange = [-1,1]*max(tdis,/nan)
        yrange = xrange
        xtitle = 'Para '+xtitl
        ytitle = 'Perp '+xtitl
        ticklen = -0.03
        xticklen = ticklen
        yticklen = ticklen
        if time_id eq 0 then begin
            ytickformat=''
        endif else begin
            ytickformat='(A1)'
            ytitle = ' '
        endelse
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
            ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
            noerase=1, nodata=1, position=tpos, color=sgcolor('black')
        tx = tpos[0]+xchsz*0.25
        ty = tpos[3]-ychsz*0.7
        xyouts, tx,ty,/normal, 'a-'+string(time_id+1,format='(I0)'), color=sgcolor('black')

        tx = tpos[0]
        ty = tpos[3]+ychsz*0.3
        ;msg = 'RBSP-'+strupcase(probe)+' Electron!C'+$
        msg = $
            time_string(ut0-type_dt,tformat='hh:mm:ss')+' - '+$
            time_string(ut0+type_dt,tformat='hh:mm:ss')+' UT'
        xyouts, tx,ty,/normal, msg, charsize=label_size, color=sgcolor('black')
    endforeach





    ; En spec.
    tvar = prefix+'e_en_spec'
    yrange = [15,4e4]
    yminor = 10
    log_ytickv = make_bins(alog10(yrange), 1, /inner)
    ytickv = 10^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)')
    options, tvar, 'ytitle', 'Energy!C(eV)'
    options, tvar, 'yrange', yrange
    options, tvar, 'yticks', yticks
    options, tvar, 'ytickv', ytickv
    options, tvar, 'yminor', yminor
    options, tvar, 'ytickname', ytickn
    options, tvar, 'zcharsize', label_size
    zrange = [1e4,1e10]
    log_zrange = alog10(zrange)
    ztickv = 10^make_bins(log_zrange, 1)
    zticks = n_elements(ztickv)
    ztickn = '10!U'+string(alog10(ztickv),format='(I0)')
    ztickn[0:*:2] = ' '
    options, tvar, 'zrange', zrange
    options, tvar, 'zticks', zticks
    options, tvar, 'ztickv', ztickv
    options, tvar, 'ztickname', ztickn


    ; Pitch angle spec.
    options, prefix+'e_pa_spec', 'yrange', [0,180]
    options, prefix+'e_pa_spec', 'yticks', 2
    options, prefix+'e_pa_spec', 'yminor', 3
    options, prefix+'e_pa_spec', 'ylog', 0
    options, prefix+'e_pa_spec', 'ystyle', 1


    ; Density.
    efw_density = get_var_data(prefix+'efw_density', times=times, limits=limits)
    copy_data, prefix+'emfisis_density', prefix+'tmp_density'
    interp_time, prefix+'tmp_density', times
    emfisis_density = get_var_data(prefix+'tmp_density')
    copy_data, prefix+'ele_n', prefix+'tmp_density'
    interp_time, prefix+'tmp_density', times
    hope_density = get_var_data(prefix+'tmp_density')
    tvar = prefix+'combo_density'
    store_data, tvar, times, [[efw_density],[emfisis_density],[hope_density]], limits=limits
    options, tvar, 'labels', ['EFW','Emfisis','HOPE-e']
    options, tvar, 'colors', sgcolor(['black','red','silver'])
    options, tvar, 'yrange', [1,20]
    options, tvar, 'ylog', 1

    ; temperature.
    rbsp_read_ele_temp, time_range, probe=probe
    e_temp = get_var_data(prefix+'e_temp', at=times)
    p_temp = get_var_data(prefix+'p_temp', at=times)
    o_temp = get_var_data(prefix+'o_temp', at=times)
    tvar = prefix+'combo_temp'
    store_data, tvar, times, [[e_temp],[p_temp],[o_temp]], limits={$
        labels: ['e','H+','O+'], $
        colors: constant('rgb'), $
        ytitle: '(eV)', $
        ytickv: [1e3,1e4], $
        ytickname: '10!U'+['3','4'], $
        yticks: 1, $
        yminor: 10, $
        yrange: [5e2,2e4], $
        ylog: 1}

    ; pressure.
    ; Calc full density.
    dens = get_var_data(prefix+'efw_density', times=times, limits=lim, in=time_range)
    ntime = n_elements(times)
    store_data, prefix+'e_full_density', times, dens, limits=lim
    p_dens = get_var_data(prefix+'p_density', at=times, limits=p_lim)
    o_dens = get_var_data(prefix+'o_density', at=times, limits=o_lim)
    store_data, prefix+'p_full_density', times, dens*p_dens/(p_dens+o_dens), limits=p_lim
    store_data, prefix+'o_full_density', times, dens*o_dens/(p_dens+o_dens), limits=o_lim

    foreach species, ['e','p','o'] do begin
        rbsp_read_hope_moments, time_range, probe=probe, species=species
        dens = get_var_data(prefix+species+'_density', at=times)
        tavg = get_var_data(prefix+species+'_t_avg', at=times)
        p_thermal = dens*tavg*1e6*1.6e-19*1e9 ; nPa.

        full_dens = get_var_data(prefix+species+'_full_density')
        tavg = (species eq 'e')? 200: 50
        dp = (full_dens-dens)*tavg*1e6*1.6e-19*1e9 ; nPa.
        p_thermal += dp

        store_data, prefix+species+'_p_thermal', times, p_thermal, limits={$
            ytitle: '(nPa)', labels:species+' P thermal'}
    endforeach

    b_gsm = get_var_data(prefix+'b_gsm', at=times)
    bmag = snorm(b_gsm)
    p_mag = (bmag*1e-9)^2*0.5/(4*!dpi*1e-7)*1e9 ; nPa.
    store_data, prefix+'p_mag', times, p_mag, limits={$
        ytitle: '(nPa)', labels:'P mag'}


    p_total = fltarr(ntime)
    foreach var, prefix+['p_mag',['e','p','o']+'_p_thermal'] do begin
;        interp_time, var, to=var0
        p_total += get_var_data(var)
    endforeach
    store_data, prefix+'p_total', times, p_total, limits={ytitle:'(nPa)', labels:'P total'}

    var = prefix+'combo_pressure'
    the_vars = prefix+['p_mag',['e','p','o']+'_p_thermal','p_total']
    stplot_merge, the_vars, newname=var, $
        labels=['P!DB','P '+['e-','H+','O+'], 'P Total'], $
        colors=sgcolor(['red','blue','green','purple','black'])
    options, var, 'ytitle', '(nPa)'
    options, var, 'yrange', [-1,5]+0.5
    options, var, 'ytickv', [0,4]
    options, var, 'yticks', 1
    options, var, 'yminor', 4
    options, var, 'ystyle', 1


    ; Beta.
    p_mag = get_var_data(prefix+'p_mag', at=times)
    pressure = get_var_data(prefix+'combo_pressure', at=times)
    ndim = n_elements(pressure[0,*])
    beta = fltarr(ntime,ndim)
    for ii=0,ndim-1 do begin
        beta[*,ii] = pressure[*,ii]/p_mag
        if ii eq ndim-1 then beta[*,ii] -= 1
    endfor
    var = prefix+'beta'
    str_beta = tex2str('beta')
    store_data, var, times, beta[*,1:*], limits={ytitle:'(#)', constant:1}
    options, var, 'labels', str_beta+' '+['e-','H+','O+','Total']
    options, var, 'colors', sgcolor(['blue','green','purple','black'])
    options, var, 'yrange', [0,2]-0.5
    options, var, 'ytickv', [0,2]
    options, var, 'yticks', 1
    options, var, 'yminor', 4


    ; MHD velocity.
    p_v_fac = get_var_data(prefix+'p_v_fac', at=times)
    o_v_fac = get_var_data(prefix+'o_v_fac', at=times)
    p_dens = get_var_data(prefix+'p_density', at=times)
    o_dens = get_var_data(prefix+'o_density', at=times)
    p_ratio = p_dens/(p_dens+o_dens)
    o_ratio = 1-p_ratio
    p_mass = 1d
    o_mass = 16d
    avg_mass = p_ratio*p_mass+o_ratio*o_mass
    v_fac = fltarr(ntime,3)
    for ii=0,2 do begin
        v_fac[*,ii] = (p_ratio*p_mass*p_v_fac[*,ii]+o_ratio*o_mass*o_v_fac[*,ii])/avg_mass
    endfor
    store_data, prefix+'vi', times, v_fac, limits={$
        ytitle: '(km/s)', $
        ystyle: 1, $
        yrange: [-30,30], $
        yticks: 2, $
        yminor: 3, $
        labels: 'v!D'+labfac, $
        colors: rgb }



;---Plot other data.
    margins = [10,4,10,3]
    poss = sgcalcpos(nvar, region=the_poss[*,1], margins=margins, ypans=plot_ypans)
    tplot, plot_vars, trange=time_range, position=poss, noerase=1

    tpos = poss[*,0]
    xrange = time_range
    yrange = [0,1]
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos
    foreach time, pa_times, time_id do begin
        tmp = convert_coord(time,yrange[1], data=1, to_normal=1)
        ty = tmp[1]
        tx = tmp[0]
        plots, tx+[0,0], ty+[0,0.3]*ychsz, normal=1, color=sgcolor('black')
        msg = 'a-'+string(time_id+1,format='(I0)')
        xyouts, tx,ty+ychsz*0.5,normal=1, alignment=0.5, msg, color=sgcolor('black')
    endforeach

    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = xchsz*3
        ty = tpos[3]-ychsz*0.35
        xyouts,tx,ty,normal=1, fig_labels[ii], color=sgcolor('black')
    endfor

    if keyword_set(test) then stop
    sgclose
end
