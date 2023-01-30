;+
; Show the arase data.
;-

function fig_2013_0501_arase_v01, plot_file, event_info=event_info

test = 1

;---Load data.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
    plasma_param = event_info['plasma_param']
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,'fig_2013_0501_arase_v01.pdf'])

;---Settings.
    prefix = 'arase_'
    fac_labels = event_info['fac_labels']
    time_range = time_double(['2017-05-20/13:00','2017-05-20/14:00'])+[-1,1]*3600
    bar_times = make_bins(time_range,600, inner=1)
    label_size = 0.8


;---Load data.
    time_step = 1.
    common_times = make_bins(time_range, time_step)

    b_gsm_var = arase_read_bfield(time_range, coord='gsm', resolution='64hz')
    interp_time, b_gsm_var, common_times


    r_gsm_var = arase_read_orbit(time_range)
    bmod_var = geopack_read_bfield(time_range, r_var=r_gsm_var,t89_par=2, igrf=1)
    calc_db, b_gsm_var, bmod_var, smooth=1800d
    b0_gsm_var = prefix+'b0_gsm'
    b1_gsm_var = prefix+'db_gsm'

    vinfo = geopack_trace_to_ionosphere(r_gsm_var, models='t89', igrf=1)
    get_data, b_gsm_var, times, b_gsm
    bf_gsm_var = prefix+'bf_gsm_t89'
    bf_gsm = get_var_data(bf_gsm_var, at=times)
    cmap = snorm(bf_gsm)/snorm(b_gsm)
    cmap_var = prefix+'cmap'
    store_data, cmap_var, times, cmap

    ; Convert B0 from gsm to dsi.
    b0_j2000_var = prefix+'b0_j2000'
    b0_dsi_var = prefix+'b0_dsi'
    get_data, b0_gsm_var, times, b0_gsm
    timespan, time_range[0], total(time_range*[-1,1]), seconds=1
    spd_cotrans, b0_gsm_var, b0_j2000_var, in_coord='gsm', out_coord='j2000'
    erg_cotrans, b0_j2000_var, b0_dsi_var, in_coord='j2000', out_coord='dsi'

    ; Calculate EdotB.
    e_dsi_var = arase_read_efield(time_range, coord='dsi', no_edotb=1)
    get_data, e_dsi_var, times, e_dsi
    b_dsi = get_var_data(b0_dsi_var, at=times)

    edotb_dsi_var = prefix+'edotb_dsi'
    e_dsi[*,2] = -(e_dsi[*,0]*b_dsi[*,0]+e_dsi[*,1]*b_dsi[*,1])/b_dsi[*,2]
    store_data, edotb_dsi_var, times, e_dsi
    edotb_gsm_var = prefix+'edotb_gsm'
    edotb_j2000_var = prefix+'edotb_j2000'
    erg_cotrans, edotb_dsi_var, edotb_j2000_var, in_coord='dsi', out_coord='j2000'
    spd_cotrans, edotb_j2000_var, edotb_gsm_var, in_coord='j2000', out_coord='gsm'


    add_setting, b1_gsm_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'unit', 'nT', $
        'coord', 'GSM' )
    add_setting, edotb_gsm_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'E', $
        'unit', 'mV/m', $
        'coord', 'GSM' )


    ; Define FAC.
    define_fac, b0_gsm_var, r_gsm_var
    db_fac_var = prefix+'db_fac'
    to_fac, b1_gsm_var, output=db_fac_var
    de_fac_var = prefix+'de_fac'
    to_fac, edotb_gsm_var, output=de_fac_var

    interp_time, de_fac_var, common_times
    interp_time, db_fac_var, common_times


    ; Calculate Poynting flux.
    pf_fac_var = prefix+'pf_fac'
    stplot_calc_pflux_mor, de_fac_var, db_fac_var, pf_fac_var


    add_setting, de_fac_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dE', $
        'unit', 'mV/m', $
        'yrange', [-1,1]*14, $
        'ytickv', [-1,0,1]*10, $
        'yticks', 2, $
        'yminor', 5, $
        'colors', constant('rgb'), $
        'constant', 0, $
        'coord', 'FAC', $
        'coord_labels', fac_labels )
    add_setting, db_fac_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'dB', $
        'yrange', [-1,1]*14, $
        'ytickv', [-1,0,1]*10, $
        'yticks', 2, $
        'yminor', 5, $
        'unit', 'nT', $
        'colors', constant('rgb'), $
        'constant', 0, $
        'coord', 'FAC', $
        'coord_labels', fac_labels )
    add_setting, pf_fac_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'S', $
        'unit', 'mW/m!U2!N', $
        'yrange', [-3,9], $
        'ytickv', [0,5], $
        'yticks', 1, $
        'yminor', 5, $
        'colors', constant('rgb'), $
        'constant', 0, $
        'coord', 'FAC', $
        'coord_labels', fac_labels ) 

    pf_fac = get_var_data(pf_fac_var, at=common_times, limits=lim)
    pf_map_var = prefix+'pf_map_fac'
    pf_map = pf_fac
    cmap = get_var_data(cmap_var, at=common_times)
    ndim = 3
    for ii=0,ndim-1 do pf_map[*,ii] *= cmap
    store_data, pf_map_var, common_times, pf_map, limits=lim


    ; Calculate the E/B ratio.
    ebr_var = stplot_calc_ebratio(time_range, e_var=de_fac_var, b_var=db_fac_var)


    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, xsize=6, ysize=4, xchsz=xchsz, ychsz=ychsz

    spec_var = prefix+'pf_fac_mor_spec_1'
    spec_var1 = prefix+'pf_fac_mor_spec_1_plot'
    vars = [prefix+['de_fac','db_fac','pf_map_fac'],spec_var1]
    nvar = n_elements(vars)
    yrange = [0.5,100]  ; mHz.

    get_data, spec_var, times, pfspec, ps, limits=lim
    store_data, spec_var1, times, pfspec*1e3, 1d3/ps, limits=lim
    add_setting, spec_var1, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'S!D||!N', $
        'color_table', 66, $
        'zrange', [-1,1]*3, $
        'ztickv', [-2,0,2], $
        'zticks', 2, $
        'zminor', 2, $
        'zticklen', -0.3, $
        'ytitle', 'Freq!C(mHz)', $
        'yrange', yrange, $
        'unit', tex2str('mu')+'W/m!U2!N' )
        
    poss = sgcalcpos(nvar, margins=[10,4,25,1])
    fig_labels = letters(nvar+1)+')'

    tplot, vars, trange=time_range, position=poss
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        tx = tpos[0]-xchsz*8
        ty = tpos[3]-ychsz*0.7
        xyouts, tx,ty,normal=1, fig_labels[ii]
    endfor

    tpos = poss[*,2]
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.9
    xyouts, tx,ty,normal=1, 'Norm to 100 km altitude!CFiltered in 0.5mHz-0.2Hz', charsize=label_size


    tpos = poss[*,nvar-1]
    tpos[0] = tpos[2]+xchsz*10
    tpos[2] = 1-xchsz*2.5
    get_data, ebr_var, ps, ebr
    xrange = [100,1e4]
    xtitle = 'E/B ratio (km/s)'
    plot, ebr, 1d3/ps, $
        ystyle=1, ylog=1, yrange=yrange, $
        xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, $
        position=tpos, noerase=1, $
        xticklen=-0.02, yticklen=-0.01

    bmag = mean(snorm(get_var_data(prefix+'b0_gsm',in=time_range)),nan=1)
    n0 = 1
    m0 = 2.5
    va = 22.0*bmag/sqrt(n0*m0)
    plots, va+[0,0], yrange, linestyle=1
    tmp = convert_coord(va, yrange[0], data=1, to_normal=1)
    tx = tmp[0]
    ty = tmp[1]-ychsz*label_size*0.9
    msg = 'v!DA!N'
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=sgcolor('red')
    xyouts, tpos[0], tpos[3]+ychsz*0.5, fig_labels[-1], normal=1
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    xyouts, tx,ty-ychsz*0,normal=1, '|B| = '+string(bmag,format='(I0)')+' nT', charsize=label_size
    xyouts, tx,ty-ychsz*1,normal=1, 'n ~ '+string(n0,format='(I0)')+' cm!U-3!N', charsize=label_size
    xyouts, tx,ty-ychsz*2,normal=1, 'm ~ '+string(m0,format='(I0)')+' m!DH!N', charsize=label_size


    if keyword_set(test) then stop
    sgclose



end

print, fig_2013_0501_arase_v01(event_info=event_info)
end