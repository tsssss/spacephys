;+
; Test Tracers ephemeris from SPICE kernel and from L1B.
;-


function test_tracers_ephemeris_spice_and_l1b, time_range, probe=probe, test=test
    compile_opt idl2

    prefix = 'ts'+probe+'_'
    default_coord = 'geo'   ; To be determined. Yes.
    r_var_default = prefix+'r_'+default_coord

;---From SPICE.
    tracers_load_spice_kernel, time_range, probe=probe

    r_var_spice = r_var_default+'_spice'
    time_step = 1.
    times = make_bins(time_range, time_step)
    ut0 = time_string(times[0],tformat='YYYY-MM-DDThh:mm:ss')
    cspice_str2et, ut0, et0
    ets = et0+times-ut0
    uts = time_string(times,tformat='YYYY-MM-DDThh:mm:ss.ffffff')
    cspice_str2et, uts, ets
    target = 'TS'+probe
    observer = 'EARTH'
    frame = 'GEO'
    abcoor = 'NONE'
    cspice_spkezr, target, ets, frame, abcoor, observer, state, local_time
    r_default = transpose(state[0:2,*])
;    v_geo = transpose(state[3:5,*])
    re = constant('re')
    settings = dictionary('coord',default_coord)
    r_var_spice = var_store(r_var_spice, r_default/re, times, id='position', settings=settings)
;    v_var = var_store(prefix+'v_geo', v_geo, times, id='velocity', settings=settings)

;---From L1B.
    r_var_l1b = r_var_default+'_l1b'
    files = tracers_load_aci(time_range, probe=probe, errmsg=errmsg, id='l1b%ipd_x292')

    var_list = list()
    in_vars = 'location'
    time_var = 'Epoch'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2000' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg

    locs = var_get_data(in_vars, times=times)
    re = constant('re')
    rad = constant('rad')
    diss = locs[*,0]/re
    colats = locs[*,1]*rad
    lons = locs[*,2]*rad
    ntime = n_elements(times)
    ndim = 3
    r_default = fltarr(ntime,ndim)
    r_default[*,2] = diss*cos(colats)
    r_xys = diss*sin(colats)
    r_default[*,0] = r_xys*cos(lons)
    r_default[*,1] = r_xys*sin(lons)
    settings = dictionary('coord',default_coord)
    r_var_l1b = var_store(r_var_l1b, r_default, times, id='position', settings=settings)

  
;---Plot.
    plot_dir = srootdir()
    time_str = time_string(time_range[0],tformat='YYYY_MMDD_hh')
    version_str = 'v01'
    plot_file = join_path([plot_dir,'test_tracers_ephemeris_spice_and_l1b'+time_str+'_'+version_str+'.pdf'])
    if keyword_set(test) then plot_file = 0
    fig_size = [12,8]
    sgopen, plot_file, size=fig_size

    r_vars = [r_var_spice, r_var_l1b]
    options, r_vars, yrange=[-1,1]*1.2, constant=[0]

    ndim = 3
    r_vec_l1b = var_get_data(r_var_l1b, times=times, settings=settings)
    r_vec_spice = var_get_data(r_var_spice, at=times)
    dr_vec = r_vec_l1b - r_vec_spice
    components = constant('xyz')
    colors = settings['colors']
    for ii=0,ndim-1 do begin
        var = prefix+'dr_'+default_coord+'_'+components[ii]
        var = var_store(var, dr_vec[*,ii], times)
        add_setting, var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', settings['unit'], $
            'colors', colors[ii], $
            'psym', 3, $
            'short_name', strupcase(default_coord)+' '+'dR'+components[ii] )
    endfor
    dr_vars = prefix+'dr_'+default_coord+'_'+components
    options, dr_vars, yrange=[-1,1]*0.1, constant=[0]

    plot_vars = [r_vars,dr_vars]
    tplot, plot_vars, trange=time_range
    if keyword_set(test) then stop
    sgclose


;---Time shift.
    plot, deriv(times, r_vec_spice[*,0]), dr_vec[*,0], nodata=1
    fit_xxs = []
    fit_yys = []
    for ii=0,ndim-1 do begin
        txs = deriv(times, r_vec_spice[*,ii])
        tys = dr_vec[*,ii]
        fit_xxs = [fit_xxs, txs]
        fit_yys = [fit_yys, tys]
        plots, txs, tys, psym=3, color=colors[ii]
    endfor
    res = linfit(fit_xxs, fit_yys)
    dt = -res[1]
    plots, fit_xxs, fit_xxs*dt, color=sgcolor('purple')

    plot_file = join_path([plot_dir,'test_tracers_ephemeris_spice_and_l1b_time_shift_'+time_str+'_'+version_str+'.pdf'])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    msg = 'Time shift on SPICE: '+strtrim(string(dt, format='(F8.1)'),2)+' s'

    r_vec_spice = var_get_data(r_var_spice, times=times, settings=settings)
    times = times+dt
    r_var_spice_shifted = var_store(r_var_spice+'_shifted', r_vec_spice, times, id='position', settings=settings)

    r_vec_l1b = var_get_data(r_var_l1b, times=times, settings=settings)
    r_vec_spice = var_get_data(r_var_spice_shifted, at=times)
    dr_vec = r_vec_l1b - r_vec_spice
    components = constant('xyz')
    colors = settings['colors']
    for ii=0,ndim-1 do begin
        var = prefix+'dr_'+default_coord+'_'+components[ii]
        var = var_store(var, dr_vec[*,ii], times)
        add_setting, var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', settings['unit'], $
            'colors', colors[ii], $
            'psym', 3, $
            'short_name', strupcase(default_coord)+' '+'dR'+components[ii] )
    endfor
    dr_vars = prefix+'dr_'+default_coord+'_'+components
    options, dr_vars, yrange=[-1,1]*0.1, constant=[0]

    r_vars = [r_var_spice_shifted, r_var_l1b]
    plot_vars = [r_vars,dr_vars]
    tplot, plot_vars, trange=time_range, get_plot_position=plot_poss
    tpos = combine_pos(plot_poss)
    tx = tpos[0]+xchsz*0.0
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty, msg, normal=1
    if keyword_set(test) then stop
    sgclose
    

    return, plot_file
end


; file = '/Volumes/data/tracers/SOC/spice/metakernels/predict_metakernel_current.tm'
; cspice_furnsh, file
; 
; ; Check currently loaded kernels.
; cspice_ktotal, 'ALL', nkernel
; kernels = list()
; for ii=0,nkernel-1 do begin
;     cspice_kdata, ii, 'ALL', file, file_type, source, handle, found
;     kernels.add, file
; end
; stop
; 
; tplot_time=time_double('2025-08-05/00:00:00')
; cspice_str2et, '2025-08-05T00:00:00', EphemTime
; tplot_to_spice_time=ephemtime-tplot_time
; t_1sec = tplot_time+[0,1]
; ets = t_1sec+tplot_to_spice_time
; SC = 'TS2'
; cspice_spkezr,SC,ets,'GEO','NONE','EARTH',ts_geo,lttime
; 
; pos_geo={x:t_1sec,y:transpose(ts_geo[0:2,*])}
; vel_geo={x:t_1sec,y:transpose(ts_geo[3:5,*])}
; 
; tplot_time=time_double('2025-08-05/00:00:00')
; cspice_str2et, '2025-08-05T00:00:00', EphemTime
; tplot_to_spice_time=ephemtime-tplot_time
; cspice_pxform, 'GEO', 'TS2_TSS',ets, RotMatGEOtoTSS
; 
; 
; 
; stop


compile_opt idl2
probe = '2'
test = 1
time_range = time_double(['2025-12-11/17:00','2025-12-11/21:00'])
dates = ['2025-08-30','2025-12-01','2026-01-01']
foreach date, dates do begin
    time_range = time_double(date+'/'+['00:00','02:00'])
    print, test_tracers_spin_axis_direction(time_range, probe=probe, test=test)
    stop
    print, test_tracers_ephemeris_spice_and_l1b(time_range, probe=probe, test=test)
endforeach



end