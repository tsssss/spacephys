;+
; Figure is the same as v03, just to update the codes and added comments.
;-
function micro_injection_fig_pad_for_noon_v01, input_event_id, probe=probe, $
    plot_dir=plot_dir, test=test, get_name=get_name, update=update

    errmsg = ''
    retval = !null
    version = 'v01'
    project= micro_injection_load_project()
    project_id = project.id

    if n_elements(input_event_id) eq 2 then begin
        time_range = time_double(input_event_id)
        event_id = time_string(time_range[0],tformat='YYYY_MMDD_hh')
    endif else begin
        event_id = input_event_id
    endelse
    event = project_get_event(project, id=event_id)
    time_range = event.time_range
    time_range = time_double(['2016-03-05/17:30','2016-03-05/19:40'])
    if n_elements(event) eq 0 then message, 'Inconsistency ...'

    if n_elements(plot_dir) eq 0 then plot_dir = event.plot_dir
    base = project_id+'_fig_pad_'+event_id+'_mms'+probe+'_'+version+'.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(get_name) then return, plot_file
    print, plot_file
    if keyword_set(update) then file_delete, plot_file, allow_nonexist=1
    if keyword_set(test) then begin
        plot_file = 0
    endif else begin
        if file_test(plot_file) eq 1 then begin
            print, plot_file+' exists, skip ...'
            return, plot_file
        endif
    endelse


;---Load data.
    prefix = 'mms'+probe+'_'
    default_coord = 'gsm'
    mission_probe = 'mms'+probe

    ; Particle related vars.
    ele_pad_var = lets_read_this(func='mms_read_pad_ele_all', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load e thermal or kev data ...'
        return, retval
    endif

    ion_pad_var = lets_read_this(func='mms_read_pad_ion_all', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion thermal data ...'
        return, retval
    endif
    ion_kev_pad_var = lets_read_this(func='mms_read_pad_ion_kev', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion kev data ...'
        return, retval
    endif
    ion_pad_var = lets_read_this(func='mms_read_pad_ion_thermal_fpi', $
        time_range, probe=mission_probe, errmsg=errmsg)
    ion_pad_var = rename_var(ion_pad_var, output=prefix+'p_pad_thermal')


;---Settings.
    unit_type = 'nflux'
    ele_pa_en_low = 50
    ele_pa_en = 800
    ion_pa_en_low = 80
    ion_pa_en = 3000
    ion_pa_high = 2.5e4
    ele_pa_high = ion_pa_high

    
;---Derived data.
    suffix = '_'+['thermal','kev','all']
    pad_vars = prefix+['e_pad'+suffix,'p_pad'+suffix]
    unit0 = '/cm!U2!N-s-sr-keV'
    foreach pad_var, pad_vars do begin
        var1 = pad_var+'_flux'

        ; Adjust unit.
        data = get_var_data(pad_var, times=times, limits=lim)
        vals = lim.en_centers
        foreach val, vals, vid do begin
            if unit_type eq 'nflux' then begin
                ; do nothing.
                unit = '#'+unit0
            endif else if unit_type eq 'eflux' then begin
                data[*,*,vid] *= val*1e-3
                unit = 'keV'+unit0
            endif else if unit_type eq 'xflux' then begin
                data[*,*,vid] *= sqrt(val*1e-3)
                unit = 'keV!U0.5!N'+unit0
            endif
        endforeach

        ; adjust flux for ion.
        if pad_var eq prefix+'p_pad_all' then begin
            index = where(vals ge 2.5e4)
            data[*,*,index] *= 5
        endif

        ; remove some useless energy bins
        index = where(vals le 1.5e5)
        data = data[*,*,index]
        vals = vals[index]

        ; save data.
        store_data, var1, times, data, limits=lim
        options, var1, unit=unit, en_centers=vals
    endforeach


    ; en spec.
    en_spec_vars = list()
    pa_spec_vars = list()
    foreach pad_var, pad_vars do begin
        en_spec_vars.add, pad_get_en_spec(pad_var=pad_var)
        var = pad_get_pa_spec(pad_var=pad_var)
        pa_spec_vars.add, var
        options, var, 'energy_range', minmax(get_var_setting(pad_var,'en_centers'))
    endforeach
    en_spec_vars = en_spec_vars.toarray()
    pa_spec_vars = pa_spec_vars.toarray()

    ; pa spec.
    ; electron low and mid.
    pad_var = prefix+'e_pad_thermal'
    ele_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ele_pa_en_low,ele_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range
    ; electron high.
    energy_range = [ele_pa_en,ele_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    
    ; ion low and mid.
    pad_var = prefix+'ion_pad_thermal'
    ion_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ion_pa_en_low,ion_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range
    ; ion high.
    energy_range = [ion_pa_en,ion_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    

;---Settings for en and pa spec.
    log_ytickv = [2,3,4]
    yticks = n_elements(log_ytickv)-1
    zrange = [1e3,1e8]
    log_ztickv = [3,4,5,6,7,8]
    zticks = n_elements(log_ztickv)-1
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    ztickn[0:*:2] = ' '
    vars = prefix+'e_en_spec_thermal'
    options, vars, zrange=zrange, constant=[ele_pa_en_low,ele_pa_en], $
        ytickv=10d^log_ytickv, yticks=yticks, ytickname='10!U'+string(log_ytickv,format='(I0)'), $
        ztickv=10d^log_ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9
        
    vars = prefix+'e_en_spec_kev'
    log_ytickv = [4,5]
    yfactor = 5
    yrange = yfactor*10d^log_ytickv
    ytickn = string(yfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ytickv,format='(I0)')
    ytickv = 1e5
    yticks = n_elements(ytickv)-1
    ytickn = '10!U5'
    log_zrange = [1,5]-3
    zfactor = 5
    zrange = zfactor*10d^log_zrange
    log_ztickv = make_bins(log_zrange,1, inner=1)
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ztickv,format='(I0)')
    index = where(log_ztickv eq 1, count)
    if count ne 0 then ztickn[index] = '10'
    ztickn[1:*:2] = ' '
    options, vars, zrange=zrange, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, ytickname=ytickn, $
        ztickv=ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9

    vars = prefix+'ion_pa_spec_thermal_low'
    zfactor = 1.3
    log_zrange = [4,6]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, vars, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    
    vars = prefix+'ion_pa_spec_thermal_high'
    zfactor = 8
    log_zrange = [4,5]-1
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, vars, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    vars = prefix+'p_en_spec_thermal'
    log_zrange = [4,6]
    zfactor = 1
    zrange = zfactor*10d^log_zrange
    log_yrange = [1,4]
    yfactor = 2
    yrange = yfactor*10d^log_yrange
    log_ytickv = make_bins(minmax(alog10(yrange)),1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    options, vars, constant=[ion_pa_en_low,ion_pa_en], $
        yrange=yrange, ytickv=ytickv, yticks=yticks, ytickname=ytickn, yminor=9, $
        zrange=zrange

;---Settings for electron PA specs.
    var = prefix+'e_pa_spec_kev'
    zfactor = 3.5
    log_zrange = [2,3]-2
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    var = prefix+'e_pa_spec_thermal_high'
    zfactor = 2.8
    log_zrange = [5,6]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    var = prefix+'e_pa_spec_thermal_low'
    zfactor = 1
    log_zrange = [6,8]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    
;---Load more data
    ; fields.
    ; B field related vars.
    field_time_range = time_range+[-1,1]*30.*60
    b_gsm_var = lets_read_this(func='mms_read_bfield', $
        field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load B field data ...'
        return, retval
    endif
    e_gsm_var = lets_read_this(func='mms_read_efield', $
        field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load E field data ...'
        return, retval
    endif


    ; Orbit related vars.
    orbit_time_range = time_range+[-1,1]*30*60
    r_gsm_var = lets_read_this(func='mms_read_orbit', $
        orbit_time_range, probe=mission_probe, coord=default_coord)
    print, 'Loading '+r_gsm_var+' ...'
    options, r_gsm_var, mission_probe='mms'+probe, coord='gsm'
    mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
    foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'

    ; Model related vars.
    external_models = ['t89','t96','t01','t04s']
    internal_models = ['dipole','igrf']
    hemispheres = ['north','south']

    foreach external_model, external_models do begin
        foreach internal_model, internal_models do begin
            ; The B at sc position.
            suffix = '_'+internal_model+'_'+external_model
            bmod_var = prefix+'bmod_gsm'+suffix
            if tnames(bmod_var) ne '' then del_data, bmod_var
            bmod_var = lets_read_geopack_bfield(var_info=bmod_var, $
                orbit_var=r_gsm_var, time_var=orbit_time_var, $
                internal_model=internal_model, external_model=external_model, save_to=data_file, update=update)
            print, 'Loading '+bmod_var+' ...'
        endforeach
    endforeach
    
    
    ; B model.
    b0_window = 15.*60
    bmod_var = prefix+'bmod_gsm_igrf_t89'
    b_vars = lets_decompose_bfield(b0_window=b0_window, b_var=b_gsm_var, bmod_var=bmod_var)
    b0_gsm_var = b_vars['b0']
    b_elev_var = lets_calc_vec_elev(b_gsm_var, coord='sm')
    bmod_elev_var = lets_calc_vec_elev(bmod_var, coord='sm', var_info=prefix+'bmod_elev')
    db_elev_var = lets_subtract_vars(b_elev_var, bmod_elev_var, save_to=prefix+'db_elev')
    options, db_elev_var, constant=0, yrange=[-1,1]*90

    b0_var = b_vars['b0']
    b1_var = b_vars['b1']

    ion_vel_var = lets_read_this(func='mms_read_ion_vel', $
        time_range, probe=mission_probe, errmsg=errmsg)

    ; Convert to FAC.
    options, [r_gsm_var,b0_var], mission='mms'
    q_fac_var = lets_define_fac(r_var=r_gsm_var, b_var=b0_var, fac_coord=fac_coord)
    coord_msgs = [default_coord,fac_coord]
    fac_vars = list()
    foreach var, [b1_var,e_gsm_var,ion_vel_var] do begin
        options, var, mission='mms'
        fac_vars.add, lets_cotran(coord_msgs, input=var, q_var=q_fac_var)
    endforeach
    
    get_data, prefix+'b1_mms_fac', times, b1_fac, limits=lim
    b0_mag = snorm(get_var_data(b0_var, at=times))
    b1_fac[*,0] += b0_mag
    store_data, prefix+'b_mms_fac', times, b1_fac, limits=lim
    

    ; Vexb.    
    b0_mag = snorm(get_var_data(b0_var, times=times))
    e_fac = get_var_data(prefix+'e_mms_fac', at=times)
    ntime = n_elements(times)
    ndim = 3
    b_fac = fltarr(ntime,ndim)
    b_fac[*,0] = b0_mag
    vexb_fac = vec_cross(e_fac,b_fac)
    coef = 1e3/b0_mag^2
    for ii=0,ndim-1 do vexb_fac[*,ii] *= coef
    var = prefix+'vexb_mms_fac'
    store_data, var, times, vexb_fac
    add_setting, var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'V', $
        'unit', 'km/s' )
    
    ; Settings.
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    vars = prefix+['b_mms_fac','e_mms_fac','vexb_mms_fac']
    options, vars, 'labels', fac_labels
    
    var = prefix+'b_mms_fac'
    options, var, yrange=[-15,15], ytickv=[-10,0,10], yminor=5, yticks=2, constant=[-10,0,10]
    var = prefix+'e_mms_fac'
    options, var, yrange=[-1,1]*20, ytickv=[-1,0,1]*15, yminor=5, yticks=2, constant=[0]
    var = prefix+'vexb_mms_fac'
    options, var, yrange=[-1,1]*350, ytickv=[-1,0,1]*300, yminor=3, yticks=2, constant=[-1,0,1]*200
    


;---Settings for plots.
    pa_suffix = ''
    plot_vars = prefix+[$
        'e_en_spec_kev','e_en_spec_thermal', $
        'e_pa_spec_kev','e_pa_spec_thermal_high','e_pa_spec_thermal_low', $
        'p_en_spec_thermal', $
        'ion_pa_spec_thermal_high','ion_pa_spec_thermal_low','b_mms_fac','e_mms_fac','vexb_mms_fac']
    panel_labels = [$
        'e- EN high','e- EN',$
        'e- high'+pa_suffix,'e- mid'+pa_suffix,'e- low'+pa_suffix,$
        'Ion EN', $
        'Ion mid'+pa_suffix,'Ion low'+pa_suffix, $
        'B FAC','E FAC','V!DExB!N FAC']
    
    
    tplot_options, 'tickinterval', 1200
    pansize = [6,0.75]
    fig_info = sgplot(plot_vars, panel_labels=panel_labels, $
        xrange=time_range, filename=plot_file, pansize=pansize)
    
    pa_vars = prefix+[['e','ion']+'_pa_spec_kev','e_pa_spec_thermal_'+['low','high'],'ion_pa_spec_thermal_'+['low','high']]
    panel_info = fig_info['panel_info']
    xchsz = fig_info.xchsz
    ychsz = fig_info.ychsz
    foreach var, pa_vars do begin
        if not panel_info.haskey(var) then continue
        the_info = panel_info[var]
        tpos = the_info['position']
        energy_range = get_var_setting(var,'energy_range')*1e-3
        msg = strjoin(strtrim(string(energy_range,format='(F5.1)'),2),'-')+' keV'
        index = strpos(var, 'kev')
        if index[0] ne -1 then msg = strjoin(strtrim(string(energy_range,format='(I0)'),2),'-')+' keV'
        tx = tpos[2]-xchsz*4.5
        ty = tpos[3]-ychsz*1
        polyfill, tpos[2]-[0.5,8.5,8.5,0.5,0.5]*xchsz, tpos[3]-[0.2,0.2,1.2,1.2,0.2]*ychsz, normal=1, color=sgcolor('white')
        xyouts, tx,ty,msg, normal=1, alignment=0.5, color=sgcolor('black')
    endforeach
    
    
    
    pad_times = [$
        '2015-09-01/18:12:40',$
        '2015-09-01/18:18:00',$
        '2015-09-01/18:59:20',$
        '2015-09-01/19:06:00',$
        ;'2015-09-01/19:17:40',$
        '2015-09-01/19:22:00',$
        '2015-09-01/19:27:00',$
        '2015-09-01/19:58:00',$
        '2015-09-01/20:02:40',$
        '2015-09-01/20:30:20',$
        '2015-09-01/20:38:20' ]
    event['pad_times'] = time_double(pad_times)
    if event.haskey('pad_times') then begin
        bar_times = event['pad_times']
        panel_info = fig_info.panel_info
        tpos = (panel_info[plot_vars[0]])['position']
        tpos[1] = ((panel_info[plot_vars[-1]])['position'])[1]
        yrange = [0,1]
        set_axis, position=tpos, xrange=time_range, yrange=yrange
        foreach bar_time, bar_times do begin
            txs = time_double(bar_time)+[0,0]
            plots, txs, yrange, color=sgcolor('red'), linestyle=3, data=1
        endforeach
    endif
    
    if keyword_set(test) then stop
    sgclose


    return, plot_file

end

probe = '4'
event_id = '2016_0305_17'
print, micro_injection_fig_pad_for_noon_v01(event_id, probe=probe, test=1, update=1)
end