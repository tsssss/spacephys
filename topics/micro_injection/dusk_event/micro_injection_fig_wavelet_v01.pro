;+
; Do wavelet transform for various quantities.
;-

function micro_injection_fig_wavelet_v01, input_event_id, test=test, get_name=get_name, update=update, errmsg=errmsg

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
    if n_elements(event) eq 0 then message, 'Inconsistency ...'

    if n_elements(plot_dir) eq 0 then plot_dir = event.plot_dir
    base = project_id+'_fig_wavelet_'+event_id+'_'+version+'.pdf'
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
    mission = 'mms'
    default_coord = 'gsm'
    external_model = 't89'
    internal_model = 'igrf'
    fac_coord = mission+'_fac'
    fac_labels = ['b','w','o']
    b0_window = 1200d
    probes = string(findgen(4)+1,format='(I0)')
    nprobe = n_elements(probes)
    colors = sgcolor(['red','green','blue','purple'])
    labels = strupcase('mms'+probes)
    comps = constant('xyz')
    ncomp = n_elements(comps)
    
    foreach probe, probes do begin
        prefix = mission+probe+'_'
        source = ['mms',probe]
        r_var = lets_read('orbit', time_range, source=source, coord=default_coord)
        b_var = lets_read('bfield', time_range, source=source, coord=default_coord)
        e_var = lets_read('efield', time_range, source=source, coord=default_coord)
        u_var = lets_read('ion_vel', time_range, source=source, coord=default_coord)

        bmod_var = lets_read_geopack_bfield(orbit_var=r_var, external_model=external_model, internal_model=internal_model)
        b_vars = lets_decompose_bfield(b_var=b_var, b0_window=b0_window, bmod_var=bmod_var)
        b0_var = b_vars['b0']
        b1_var = b_vars['b1']

        ; Convert to FAC.
        q_fac_var = lets_define_fac(r_var=r_var, b_var=b0_var, fac_coord=fac_coord)
        coord_msgs = [default_coord,fac_coord]
        fac_vars = list()
        foreach var, [b1_var,e_var,u_var] do begin
            fac_vars.add, lets_cotran(coord_msgs, input=var, q_var=q_fac_var)
        endforeach

        ; Calc wavelet.
        wanted_comp = 'w'
        wanted_index = where_pro(fac_labels, 'eq', wanted_comp)
        scale_info = {s0:2d, s1:4000, dj:1d/8, ns:0d }
        ct = 40
        zrange = [1e-2,1e3]
        foreach var, fac_vars do begin
            vec = get_var_data(var, times=times, settings=settings)
            dat = vec[*,wanted_index]
            field_var = var+'_tmp'
            store_data, field_var, times, dat
            spec_var = stplot_mor_new(field_var, scale_info=scale_info)
            get_data, spec_var, times, data, freqs
            store_data, spec_var, times, data, freqs*1e3
            unit = settings.unit
            short_name = settings.short_name
            add_setting, spec_var, smart=1, dictionary($
                'requested_time_range', time_range, $
                'no_interp', 1, $
                'display_type', 'spec', $
                'unit', unit, $
                'ytitle', 'Freq (mHz)', $
                'yrange', minmax(freqs*1e3), $
                'ylog', 1, $
                'zlog', 1, $
                'short_name', short_name )
            ztitle = short_name+'('+unit+')!U2!N'
            options, spec_var, ztitle=ztitle, color_table=ct, zrange=zrange
            tmp = strpos(var, 'u_mms_fac')
            if tmp[0] ne -1 then options, spec_var, zrange=zrange*10
            tmp = rename_var(spec_var,output=var+'_mor')
        endforeach
    endforeach


;---Make plot.

    prefix = 'mms4_'
    plot_vars = prefix+['b1','e','u']+'_mms_fac_mor'
    labels = ['dB','E','Ion Vel']+'!D'+tex2str('perp')+',west!N'
    
    vars = prefix+['b1','e','u']+'_mms_fac'
    options, vars, 'labels', ['||',tex2str('perp')+','+['west','out']]

    plot_vars = prefix+['b1','e','u']+'_mms_fac_mor'
    plot_vars = prefix+[$
        'b1_mms_fac'+['','_mor'], $
        'e_mms_fac'+['','_mor'], $
        'u_mms_fac'+['','_mor'] ]
    labels = ['dB','E','Ion Vel']+'!D'+tex2str('perp')+',west!N'
    comp = ''
    labels = [$
        'dB'+comp+['','!C    Morlet'], $
        'E'+comp+['','!C    Morlet'], $
        'Ion Vel'+comp+['','!C    Morlet'] ]        

    plot_tr = time_range

    tplot_options, 'tickinterval', 1200
    tmp = sgplot(plot_vars, filename=plot_file, xrange=time_range, panel_labels=labels)
    if keyword_set(test) then stop
    sgclose

    return, plot_file


end

input_event_id = '2015_0901_11'
test = 1
input_event_id = '2015_0901_10'
print, micro_injection_fig_wavelet_v01(input_event_id, test=test, update=1)
end