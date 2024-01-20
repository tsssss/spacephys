;+
; Part of fig_2015_0416_0800_overview_v03. Note: v02 is later than v03.
; DMSP inverted-V and electron eflux.
;-


function fig_2015_0416_0800_overview_v03_dmsp_panel, dmsp_poss, event_info=event_info

    label_size = 0.8
    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]

    dmsp_info = event_info.dmsp.dmspf19
    prefix = dmsp_info['prefix']
    probe = dmsp_info['probe']
    dmsp_plot_time_range = time_double(['2015-04-16/08:01','2015-04-16/08:06'])
    dmsp_orbit_time_range = time_double(['2015-04-16/08:01','2015-04-16/08:06'])
    invertedv_times = time_double('2015-04-16/'+['08:02:45','08:03:22','08:03:37'])
    invertedv_color = sgcolor('salmon')
    invertedv_text = 'Inverted-V'
    dmsp_color = dmsp_info.sc_color
    ssusi_id = 'energy'
    ssusi_wavelength = strupcase(ssusi_id)

;---Load data.
    dmsp_mlt_image_var = dmsp_read_mlt_image(dmsp_plot_time_range+[-1800,0], probe=probe, id=ssusi_id)
    mlat_vars = dmsp_read_mlat_vars(dmsp_plot_time_range, probe=probe, errmsg=errmsg)
    ele_spec_var = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    db_xyz_var = dmsp_read_bfield_madrigal(dmsp_plot_time_range, probe=probe)
    r_var = dmsp_read_orbit(dmsp_plot_time_range, probe=probe)

    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]


    ; ssusi eflux along sc track.
    min_mlat = 50.
    mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
    mlons = get_var_data(mlon_var, at=the_times)
    ; use aacgm mlon.
    mlts = aacgm_mlon2mlt(mlons, the_times)
    sc_tt = (mlts*15-90)*constant('rad')
    sc_rr = (90-abs(mlats))/(90-min_mlat)
    sc_xx = sc_rr*cos(sc_tt)
    sc_yy = sc_rr*sin(sc_tt)
    ntime = n_elements(the_times)
    ssusi_eflux = fltarr(ntime)
    
    ; mlt image.
    get_data, dmsp_mlt_image_var, times, mlt_images, limits=lim
    tmp = min(times-mean(dmsp_plot_time_range), abs=1, time_id)
    mlt_image = reform(mlt_images[time_id,*,*])
    pixel_mlat = lim.pixel_mlat
    pixel_mlt = lim.pixel_mlt
    
;    ; use ssusi time per pixel.
;    files = dmsp_load_ssusi(dmsp_plot_time_range+[-120d,0]*60, probe=probe, errmsg=errmsg)
;    file = files[0]
;    time_var = 'UT_S'
;    date = dmsp_plot_time_range[0]
;    date = date-(date mod 86400d)
;    pixel_time = netcdf_read_var(time_var, filename=files[0])*3600+date
;    pixel_time = pixel_time[*]
;    mlt_image = mlt_image[*]
;    for ii=0,ntime-1 do begin
;        tmp = min(pixel_time-the_times[ii], abs=1, index)
;        ssusi_eflux[ii] = mlt_image[index]
;    endfor
;    stop


    ; use pixel location.
    min_mlat = 50
    pixel_tt = (pixel_mlt*15-90)*constant('rad')
    pixel_rr = (90-pixel_mlat)/(90-min_mlat)
    pixel_xx = pixel_rr*cos(pixel_tt)
    pixel_yy = pixel_rr*sin(pixel_tt)

    mlt_image = mlt_image[*]
    for ii=0,ntime-1 do begin
        dis = sqrt((pixel_xx-sc_xx[ii])^2+(pixel_yy-sc_yy[ii])^2)
        tmp = min(dis[*], abs=1, index)
        ssusi_eflux[ii] = mlt_image[index]
    endfor
    ssusi_eflux_var = prefix+'ssusi_eflux'
    index = where(ssusi_eflux[1:-1]-ssusi_eflux[0:-2] ne 0)
    store_data, ssusi_eflux_var, the_times[index+1], ssusi_eflux[index+1]
    dmsp_eflux_yrange = [1e-1,1e2]
    add_setting, ssusi_eflux_var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', dmsp_eflux_yrange, $
        'display_type', 'scalar', $
        'short_name', 'Aurora eflux', $
        'unit', 'mW/m!U2!N' )


    ; map eflux and convert unit.
    get_data, prefix+'e_eflux', times, eflux, limits=lim
    var = prefix+'e_eflux_map'
    cmap = 1.4  ; this is for 800 km.
    theta_loss_cones = [45,60]
    ndim = n_elements(theta_loss_cones)
    ntime = n_elements(times)
    eflux_map = fltarr(ntime,ndim+1)
    foreach theta_loss_cone, theta_loss_cones, lc_id do begin
        sr = !dpi*sin(theta_loss_cone*constant('rad'))^2
        eflux_map[*,lc_id+1] = eflux*cmap*sr*1.6e-19*1e4*1e3  ; convert from eV/cm^2-s-sr to mW/m^2.
    endforeach
    store_data, var, times, eflux_map
    yrange = [2e-1,9e2]
    constant = 10d^[0,1,2]
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', yrange, $
        'constant', constant, $
        'display_type', 'stack', $
        'labels', ['Aurora',tex2str('theta')+'!DLC!N='+string(theta_loss_cones,format='(I0)')+'!U'+tex2str('circ')], $
        'colors', sgcolor(['red','green','blue']), $
        'ytitle', '(mW/m!U2!N)' )
    

    ; convert dB from xyz to fac.
    db_xyz = get_var_data(db_xyz_var, times=times)
    db_fac = db_xyz
    fmlt = get_var_data(mlt_var, at=times)
    fmlat = get_var_data(mlat_var, at=times)
    theta = (fmlt*15-90)*constant('rad')
    n_hat = -[[cos(theta)],[sin(theta)]]
    w_hat = [[cos(theta-0.5*!dpi)],[sin(theta-0.5*!dpi)]]
    db_fac[*,0] = n_hat[*,0]*db_xyz[*,0]+n_hat[*,1]*db_xyz[*,1]
    db_fac[*,1] = w_hat[*,0]*db_xyz[*,0]+w_hat[*,1]*db_xyz[*,1]
    db_var = prefix+'db_fac'
    store_data, db_var, times, db_fac
    add_setting, db_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', '', $
        'coord_labels', fac_labels )
    var = db_var
    yrange = [-1,1]*300
    ytickv = [-1,0,1]*200
    yticks = n_elements(ytickv)-1
    yminor = 2
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', yticks
    options, var, 'yminor', yminor


;---Plot
    prefix = dmsp_info['prefix']
    dmsp_vars = prefix+['e_en_spec','e_eflux_map']
    dmsp_labels = ['a) e-','b) KEflux']
    ndmsp_var = n_elements(dmsp_vars)
    
    
    var = prefix+'e_en_spec'
    get_data, var, times, data, vals
    index = where(data eq 0 or finite(data,nan=1), count)
    if count ne 0 then begin
        data[index] = 0.01
        store_data, var, times, data, vals
    endif
    
    yrange = [3e2,3e4]
    log_ytickv = [3,4]
    ytickn = '10!U'+string(log_ytickv,format='(I0)')+'!N'
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    yminor = 9

    zrange = [1e5,1e9]
    ztickv = 10d^[5,6,7,8,9]
    zticks = n_elements(ztickv)-1
    zminor = 9
    ztickn = '10!U'+string([5,6,7,8],format='(I0)')
    zticklen = uniform_ticklen/xchsz*1/fig_size[0]

    options, var, 'constant', ytickv
    options, var, 'yrange', yrange
    options, var, 'ytickv', ytickv
    options, var, 'yticks', n_elements(ytickv)-1
    options, var, 'yminor', 10
    options, var, 'ytickname', ytickn
    options, var, 'ytitle', 'Energy!C(eV)'
    options, var, 'zcharsize', label_size
    options, var, 'zrange', zrange
    options, var, 'ztickname', ztickn
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zminor', zminor
    options, var, 'zticklen', zticklen
    options, var, 'color_table', 64
    
    the_poss = dmsp_poss
    for pid=0,ndmsp_var-1 do begin
        tpos = the_poss[*,pid]
        xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, dmsp_vars[pid], 'xticklen', xticklen
        options, dmsp_vars[pid], 'yticklen', yticklen
    endfor
    
    tplot_options, 'version', 2
    tplot, dmsp_vars, trange=dmsp_plot_time_range, position=the_poss, noerase=1, single_line_uttick=1
    timebar, invertedv_times, linestyle=1, color=invertedv_color
    
    ; Add inverted V labels.
    pid = where(dmsp_vars eq prefix+'e_en_spec')
    tpos = the_poss[*,pid]
    plot, dmsp_plot_time_range, [0,1], $
        xstyle=5, ystyle=5, position=tpos, nodata=1, noerase=1
;    foreach invertedv_time, invertedv_times do begin
;        tmp = convert_coord(invertedv_time,1, data=1, to_normal=1)
;        tx = tmp[0]
;        plots, tx,tpos[3],normal=1, psym=8, symsize=0.5, color=invertedv_color
;    endforeach
;    invertedv_time = mean(invertedv_times)
;    tmp = convert_coord(invertedv_time,1, data=1, to_normal=1)
;    tx = tmp[0]
;    ty = tpos[3]+ychsz*0.4
;    msg = invertedv_text
;    xyouts, tx,ty,normal=1, msg, charsize=sc_label_size, alignment=0.5, color=invertedv_color

    is_juno = 0
    label_color = sgcolor('red')
    txs = invertedv_times
    tx1 = mean(minmax(txs))
    tmp = convert_coord(tx1,yrange[0], data=1, to_normal=1)
    tx0 = tmp[0]
    ty0 = (is_juno eq 1)? tpos[3]-ychsz*1.5: tpos[1]+ychsz*1.5
    ty1 = (is_juno eq 1)? ty0-ychsz*2: ty0+ychsz*2
    foreach tx,txs do begin
        tmp = convert_coord(tx,yrange[0], data=1, to_normal=1)
        tx1 = tmp[0]
        arrow, tx0,ty0, tx1,ty1, normal=1, solid=1, color=label_color, hsize=hsize, thick=thick
    endforeach
    ty2 = (is_juno eq 1)? ty0+ychsz*0.2: ty0-ychsz*1
    msg = 'Inverted-V electron'
    xyouts, tx0,ty2, msg, normal=1, alignment=0.5, color=label_color;, charsize=label_size
    
    ; Add panel label.
    for pid=0,ndmsp_var-1 do begin
        tpos = the_poss[*,pid]
        tx = tpos[0]-xchsz*8
        ty = tpos[3]-ychsz*0.5
        msg = dmsp_labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor
    
    
    ; Add SSUSI eflux.
    e_eflux_var = prefix+'e_eflux_map'
    invertedv_time_range = time_double('2015-04-16/'+['08:02:20','08:04:00'])
    pid = where(dmsp_vars eq e_eflux_var, count)
    if count ne 0 then begin
        xrange = dmsp_plot_time_range
        yrange = get_setting(e_eflux_var,'yrange')
        ylog = get_setting(e_eflux_var,'ylog', exist)
        if ~exist then ylog=0
        tpos = the_poss[*,pid]
        plot, xrange, yrange, ylog=ylog, xlog=0, $
            xstyle=5, ystyle=5, position=tpos, nodata=1, noerase=1
        ssusi_eflux = get_var_data(prefix+'ssusi_eflux', times=times, in=invertedv_time_range)
        oplot, times, ssusi_eflux, psym=6, symsize=0.25, color=sgcolor('red')
;        foreach time, invertedv_times do oplot, time+[0,0], yrange, color=invertedv_color, linestyle=1
    endif


end



function fig_2015_0416_0800_extended_data_dmsp_v01, event_info=event_info, test=test

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)

;---Plot file.
    test_plot_panels = 0
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_extended_data_dmsp_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0

;---Figure out panel size.
    fig_size = [4,2.5]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    
    margins = [9,2.5,8,1.5]
    dmsp_poss = sgcalcpos(2, margins=margins)


;---Test panels.
    if keyword_set(test_plot_panels) then begin
        panel_list = list()
        panel_list.add, cartoon_pos
        panel_list.add, config_pos
        panel_list.add, right_pos
        ;panel_list.add, left_pos
        ;panel_list.add, middle_pos
        for ii=0,nasi_panel-1 do panel_list.add, reform(asi_poss[*,ii])
        for ii=0,nleft_panel-1 do panel_list.add, reform(left_poss[*,ii+1])

        foreach tpos, panel_list do begin
            xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, $
                xtickformat='(A1)', ytickformat='(A1)', $
                xticklen=xticklen, yticklen=yticklen, $
                nodata=1, noerase=1, position=tpos
        endforeach
    endif
    

;---Panels.
    tpos = fig_2015_0416_0800_overview_v03_dmsp_panel(dmsp_poss, event_info=event_info)

;---Done.
    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end

test = 0
print, fig_2015_0416_0800_extended_data_dmsp_v01(event_info=event_info, test=test)
end