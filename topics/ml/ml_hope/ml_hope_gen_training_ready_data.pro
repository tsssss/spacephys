;+
; Generate training data for Python.
;-

pro ml_hope_gen_training_ready_data, input_time_range, $
    energy=energy, $
    history_shifts=history_shifts, $
    filename=out_file

    files = ml_hope_load_training_data(input_time_range)
    all_vars = cdf_vars(files[0])

;---Collect info for loading data.
    time_var = 'unix_time'
    if n_elements(probes) eq 0 then probes = ['a','b']
    if n_elements(all_species) eq 0 then all_species = ['p','o']
    en_spec_vars = list()
    pos_vars = list()
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        en_spec_vars.add, prefix+all_species+'_en_spec', extract=1
        pos_vars.add, prefix+['lshell','mlt','mlat'], extract=1
    endforeach
    en_spec_vars = en_spec_vars.toarray()
    pos_vars = pos_vars.toarray()

    index = where(stregex(all_vars,'omni_') ne -1)
    omni_vars = all_vars[index]

;---Load data.
    data_time_range = time_double(input_time_range)
    var_list = list()
    var_list.add, dictionary('in_vars', [en_spec_vars,pos_vars])
    read_vars, data_time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return

    if n_elements(history_shifts) eq 0 then history_shifts = [0]
    history_time_range = data_time_range+[-max(abs(history_shifts)),0]
    var_list = list()
    var_list.add, dictionary('in_vars', omni_vars)
    read_vars, history_time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return


;---Some additional processing.
    ; For flux vars.
    flux_vars = list()
    foreach en_spec_var, en_spec_vars do begin
        probe = strmid(en_spec_var,4,1)
        species = strmid(en_spec_var,6,1)
        prefix = 'rbsp'+probe+'_'
        energy_bins = ml_rbsp_hope_read_energy_bins(probe=probe, species=species)
        get_data, en_spec_var, times, en_spec
        energy_var = species+'_energy'
        energy_bins = cdf_read_var(energy_var, filename=files[0])
        energy_strs = string(round(energy_bins),format='(I0)')
        min_energy = 20
        foreach energy, energy_bins, energy_id do begin
            if energy le min_energy then continue
            energy_str = energy_strs[energy_id]
            out_var = prefix+species+'_flux_'+energy_str+'eV'
            store_data, out_var, times, en_spec[*,energy_id]
            add_setting, out_var, smart=1, dictionary($
                'display_type', 'scalar', $
                'short_name', species+' '+energy_str+'eV', $
                'unit', '#/cm!U2!N-s-sr-keV', $
                'energy', energy, $
                'ylog', 1 )
            flux_vars.add, out_var
        endforeach
    endforeach
    flux_vars = flux_vars.toarray()


    ; For pos vars.
    pos_vars = list(pos_vars, extract=1)
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        mlt_var = prefix+'mlt'
        get_data, mlt_var, times, mlt
        mlt_rad = mlt*15*constant('rad')

        pos_vars.remove, pos_vars.where(mlt_var)
        foreach type, ['cos','sin'] do begin
            var = prefix+type+'mlt'
            store_data, var, times, call_function(type,mlt_rad)
            add_setting, var, smart=1, dictionary($
                'display_type', 'scalar', $
                'unit', '#', $
                'short_name', type+'(MLT)' )
            pos_vars.add, var
        endforeach
    endforeach
    pos_vars = pos_vars.toarray()
    


    ; For omni vars.
    foreach var, omni_vars do begin
        vatt = cdf_read_setting(var, filename=files[0])
        options, var, 'unit', vatt.units
        short_name = strupcase((strsplit(var,'_',extract=1))[1:*])
        options, var, 'short_name', short_name
    endforeach
    
    omni_vars = list(omni_vars, extract=1)
    get_data, 'omni_au', times, au
    get_data, 'omni_al', times, al
    ae = au-al
    ao = (au+al)*0.5
    var = 'omni_ae'
    store_data, var, times, ae
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'unit', 'nT', $
        'short_name', 'AE' )
    omni_vars.add, var

    var = 'omni_ao'
    store_data, var, times, ao
    add_setting, var, smart=1, dictionary($
        'display_type', 'scalar', $
        'unit', 'nT', $
        'short_name', 'AO' )
    omni_vars.add, var

    xyz = constant('xyz')
    units = ['nT','km/s']
    foreach type, ['b','v'], type_id do begin
        in_var = 'omni_sw_'+type+'_gse'
        get_data, in_var, times, vec
        unit = units[type_id]
        omni_vars.remove, omni_vars.where(in_var)
        foreach comp, xyz, comp_id do begin
            out_var = 'omni_sw_'+type+comp+'_gse'
            store_data, out_var, times, vec[*,comp_id]
            add_setting, out_var, smart=1, dictionary($
                'display_type', 'scalar', $
                'unit', unit, $
                'short_name', 'SW '+strupcase(type)+comp+' GSE' )
            omni_vars.add, out_var
        endforeach
    endforeach
    omni_vars = omni_vars.toarray()
    


;---Deal with normalization.
    vars = [pos_vars,omni_vars]
    foreach var, vars do begin
        get_data, var, times, data
        avg = mean(data, nan=1)
        std = stddev(data, nan=1)
        options, var, 'avg', avg
        options, var, 'std', std
        print, var, avg, std
    endforeach


;---Deal with history.
    ; history shift = 1 sec means to use data at t0-1 at t0.
    get_data, flux_vars[0], common_times
    history_vars = list()
    foreach var, omni_vars do begin
        get_data, var, times, data, limits=lim
        foreach history_shift, history_shifts do begin
            out_var = var+'_'+string(history_shift,format='(I0)')+'s'
            store_data, out_var, common_times, interpol(data,times+history_shift, common_times)
            add_setting, out_var, smart=1, dictionary($
                'display_type', 'scalar', $
                'short_name', lim.short_name, $
                'unit', lim.unit, $
                'time_shift', time_shift, $
                'avg', lim.avg, $
                'std', lim.std )
            history_vars.add, out_var
        endforeach
    endforeach
    history_vars = history_vars.toarray()


;---Prepare to merge probes.
    prefix = 'rbsp'+probes[0]+'_'
    target_vars = sort_uniq(strmid(flux_vars,strlen(prefix)))
    probe_vars = sort_uniq(strmid(pos_vars,strlen(prefix)))
    
    uniq_vars = list()
    uniq_vars.add, history_vars, extract=1
    
    ; time of year.
    get_data, pos_vars[0], times
    time0 = time_double(time_string(times[0],tformat='YYYY'))
    toy_rad = (times-time0)/(365.25d*86400)*2*!dpi
    foreach type, ['cos','sin'] do begin
        var = type+'toy'
        data = call_function(type,toy_rad)
        store_data, var, times, data
        avg = mean(data, nan=1)
        std = stddev(data, nan=1)
        add_setting, var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', '#', $
            'avg', avg, $
            'std', std, $
            'short_name', type+'(TOY)' )
        uniq_vars.add, var
    endforeach
    uniq_vars = uniq_vars.toarray()
    
;---Merge data for probes.
    nprobe = n_elements(probes)
    get_data, uniq_vars[0], times
    ntime = n_elements(times)
    nrec = ntime*nprobe

    common_times = dblarr(nprobe,ntime)
    foreach probe, probes, probe_id do begin
        common_times[probe_id,*] = times
    endforeach
    common_times = reform(common_times, nrec)

    foreach probe_var, probe_vars do begin
        data = fltarr(nprobe,ntime)
        foreach probe, probes, probe_id do begin
            prefix = 'rbsp'+probe+'_'
            var = prefix+probe_var
            data[probe_id,*] = get_var_data(var, limits=lim)
        endforeach
        data = reform(data, nrec)
        store_data, probe_var, common_times, data, limits=lim
    endforeach

    foreach target_var, target_vars do begin
        data = fltarr(nprobe,ntime)
        foreach probe, probes, probe_id do begin
            prefix = 'rbsp'+probe+'_'
            var = prefix+target_var
            data[probe_id,*] = get_var_data(var, limits=lim)
        endforeach
        data = reform(data, nrec)
        store_data, target_var, common_times, data, limits=lim
    endforeach

    foreach uniq_var, uniq_vars do begin
        data = fltarr(nprobe,ntime)
        foreach probe, probes, probe_id do begin
            data[probe_id,*] = get_var_data(uniq_var, limits=lim)
        endforeach
        data = reform(data, nrec)
        store_data, uniq_var, common_times, data, limits=lim
    endforeach
    
    time_var = 'unix_time'
    save_vars = [uniq_vars,probe_vars,target_vars]
    if file_test(out_file) eq 1 then file_delete, out_file
    stplot2cdf, save_vars, time_var=time_var, filename=out_file
    

end

time_range = ['2012-11-01','2018-03-01']
;time_range = ['2012-11-01','2013-01-01']
history_shifts = smkarthm(0,6*24,6, 'dx')*3600
;history_shifts = smkarthm(0,6,6, 'dx')*3600
out_file = join_path([sdiskdir('data'),'sdata','ml_hope','training_ready_data',$
    'ml_hope_training_ready_data.cdf'])
ml_hope_gen_training_ready_data, time_range, history_shifts=history_shifts, filename=out_file
end