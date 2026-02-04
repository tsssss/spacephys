;+
; Read the CWT for keV electrons.
;
; energy=. Energys in eV.
;-

function micro_injection_read_cwt_kev_electron_gen_file_v01, time_range, $
    mission_probe=mission_probe, energy=energy, filename=data_file, $
    errmsg=errmsg, local_root=local_root


    compile_opt idl2
    on_error, 0
    errmsg = ''
    retval = !null

    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    energy_kev = energy*1e-3
    energy_str = string(energy_kev, format='(I0)')+'kev'

    routine = mission+'_read_kev_electron_cdaweb'
    if n_elements(time_range) eq 1 then begin
        secofday = constant('secofday')
        date = time_range[0]-(time_range[0] mod secofday)
        time_range = time_double(date)+[0,secofday]
    endif
    project_info = micro_injection_stat_load_project()
    scale_info = project_info['scale_info']
    data_time_range = time_range+[-1,1]*scale_info['s1']*4
    en_spec_var = call_function(routine, data_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval
    
    
    ; Get the flux for the wanted energy.
    fluxs = get_var_data(en_spec_var, times=times, energys, settings=settings)
    tmp = min(energys-energy_kev, abs=1, energy_index)
    if energy_index eq 0 or energy_index eq n_elements(energys) then begin
        errmsg = 'Energy '+string(energy_kev)+'keV not found in the data.'
        return, retval
    endif
    the_fluxs = fluxs[*,energy_index]
    flux_var = prefix+'kev_ele_'+energy_str
    store_data, flux_var, times, the_fluxs
    unit = settings['unit']
    options, flux_var, ylog=1, ytitle='e- flux!C('+unit+')', labels=energy_str, unit=unit
    
    log_flux_var = prefix+'log_kev_ele_'+energy_str
    var = var_store(log_flux_var, alog10(the_fluxs), times)
    options, log_flux_var, ylog=0, ytitle='Log e- flux!C('+unit+')', labels=energy_str, unit='Log10 '+unit
    
    
    time_step = project_info['common_time_step']
    common_times = make_bins(data_time_range, time_step)
    foreach var, [flux_var,log_flux_var] do begin
        flux = get_var_data(var, at=common_times)
        if var eq log_flux_var then begin
            index = where(finite(flux))
            flux = interp(flux[index],common_times[index],common_times)
        endif
        var = var_store(var, flux, common_times)
    endforeach

    ; Calculate cwt.
    time_var = 'unix_time'
    spec_var = stplot_mor_new(flux_var, scale_info=scale_info)
    options, spec_var, zlog=1, zrange=[1e0,1e8]
    log_spec_var = stplot_mor_new(log_flux_var, scale_info=scale_info)
    options, log_spec_var, zlog=1, zrange=[1e-4,1e4]
    
    
    ; Trim to the wanted time_range.
    vars = [flux_var,spec_var,log_flux_var,log_spec_var]
    foreach var, vars do begin
        data = get_var_data(var, vals, times=times, in=time_range, limits=lim)
        store_data, var, times, data, vals, limits=lim
    endforeach
    if file_test(data_file) eq 1 then file_delete, data_file
    stplot2cdf, vars, filename=data_file, time_var=time_var

    return, data_file

end

;+
; Save CWT (complex). v01 doesn't save it.
;-
function micro_injection_read_cwt_kev_electron_gen_file_v02, time_range, $
    mission_probe=mission_probe, energy=energy, filename=data_file, $
    errmsg=errmsg, local_root=local_root


    compile_opt idl2
    on_error, 0
    errmsg = ''
    retval = !null

    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    energy_kev = energy*1e-3
    energy_str = string(energy_kev, format='(I0)')+'kev'

    routine = mission+'_read_kev_electron_cdaweb'
    if n_elements(time_range) eq 1 then begin
        secofday = constant('secofday')
        date = time_range[0]-(time_range[0] mod secofday)
        time_range = time_double(date)+[0,secofday]
    endif
    project_info = micro_injection_stat_load_project()
    scale_info = project_info['scale_info']
    data_time_range = time_range+[-1,1]*scale_info['s1']*4
    en_spec_var = call_function(routine, data_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval


    ; Get the flux for the wanted energy.
    fluxs = get_var_data(en_spec_var, times=times, energys, settings=settings)
    tmp = min(energys-energy_kev, abs=1, energy_index)
    if energy_index eq 0 or energy_index eq n_elements(energys) then begin
        errmsg = 'Energy '+string(energy_kev)+'keV not found in the data.'
        return, retval
    endif
    the_fluxs = fluxs[*,energy_index]
    flux_var = prefix+'kev_ele_'+energy_str
    store_data, flux_var, times, the_fluxs
    unit = settings['unit']
    options, flux_var, ylog=1, ytitle='e- flux!C('+unit+')', labels=energy_str, unit=unit

    log_flux_var = prefix+'log_kev_ele_'+energy_str
    var = var_store(log_flux_var, alog10(the_fluxs), times)
    options, log_flux_var, ylog=0, ytitle='Log e- flux!C('+unit+')', labels=energy_str, unit='Log10 '+unit


    time_step = project_info['common_time_step']
    common_times = make_bins(data_time_range, time_step)
    foreach var, [flux_var,log_flux_var] do begin
        flux = get_var_data(var, at=common_times)
        if var eq log_flux_var then begin
            index = where(finite(flux))
            flux = interp(flux[index],common_times[index],common_times)
        endif
        var = var_store(var, flux, common_times)
    endforeach

    ; Calculate cwt.
    time_var = 'unix_time'
    spec_var = stplot_mor_new(flux_var, scale_info=scale_info)
    options, spec_var, zlog=1, zrange=[1e0,1e8]
    log_spec_var = stplot_mor_new(log_flux_var, scale_info=scale_info)
    options, log_spec_var, zlog=1, zrange=[1e-4,1e4]
    cwt_var = flux_var+'_cwt'
    log_cwt_var = log_flux_var+'_cwt'


    ; Trim to the wanted time_range.
    vars = [flux_var,spec_var,log_flux_var,log_spec_var,cwt_var,log_cwt_var]
    foreach var, vars do begin
        data = get_var_data(var, vals, times=times, in=time_range, limits=lim)
        store_data, var, times, data, vals, limits=lim
    endforeach
    if file_test(data_file) eq 1 then file_delete, data_file
    stplot2cdf, vars, filename=data_file, time_var=time_var

    return, data_file

end




function micro_injection_read_cwt_kev_electron, input_time_range, mission_probe=mission_probe, energy=energys, version=version, $
    errmsg=errmsg, update=update, get_name=get_name, log=log

    errmsg = ''
    retval = !null
    if n_elements(version) eq 0 then version = 'v02'

    project_info = micro_injection_stat_load_project()
    sample_energys = project_info['sample_energys']
    
    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    mission_probe_str = mission+'_'+probe

    if n_elements(energys) eq 0 then energys = sample_energys*1e3   ; eV.
    energy_strs = string(energys*1e-3, format='(I0)')+'kev'
    var_info = prefix+'cwt_kev_electron_'+energy_strs
    if keyword_set(log) then var_info = prefix+'cwt_log_kev_electron_'+energy_strs
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then is_success = delete_var_from_memory(var_info)
    time_range = time_double(input_time_range)
    if ~check_if_update_memory(var_info, time_range) then return, var_info


    foreach energy_str, energy_strs, eid do begin
    ;---Load files.
        base_name = mission_probe_str+'_cwt_kev_electron_'+energy_str+'_%Y_%m%d_'+version+'.cdf'
        local_root = project_info['local_root']
        local_path = [local_root,'cwt_kev_electron',mission_probe_str,'%Y',energy_str]
        request = dictionary($
            'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
            'cadence', 'day' )

        files = prepare_files(request=request, errmsg=errmsg, $
            file_times=file_times, time=time_range, nonexist_files=nonexist_files)
        if n_elements(nonexist_files) ne 0 then begin
            foreach file, request.nonexist_files do begin
                file_time = file.file_time
                local_file = file.local_file
                routine = 'micro_injection_read_cwt_kev_electron_gen_file_'+version
                file = call_function(routine, file_time, mission_probe=mission_probe, energy=energys[eid], filename=local_file)
            endforeach
            files = prepare_files(request=request, errmsg=errmsg, $
                file_times=file_times, time=time, nonexist_files=nonexist_files)
        endif


    ;---Read data.
        suffix = ['','_mor','_cwt']
        if keyword_set(log) then begin
            in_vars = prefix+'log_kev_ele_'+energy_str+suffix
            out_vars = prefix+'cwt_log_kev_electron_'+energy_str+suffix
        endif else begin
            in_vars = prefix+'kev_ele_'+energy_str+suffix
            out_vars = prefix+'cwt_kev_electron_'+energy_str+suffix
        endelse
        vars = var_read(in_vars, var_info=out_vars, file=files)
        options, vars[0], labels=energy_str
    endforeach


    return, var_info

end

missions = ['mms']

search_trs = micro_injection_load_search_time_range()
nsearch_tr = n_elements(search_trs[*,0])
dates = list()
for ii=0,nsearch_tr-1 do begin
    tr = reform(search_trs[ii,*])
    nday = total(tr*[-1,1])/constant('secofday')
    dates.add, findgen(nday)*constant('secofday')+tr[0], extract=1
endfor

foreach date, dates do begin
;if time_double(date) lt time_double('2016-01-21') then continue
    time_range = date+[0,constant('secofday')]
    foreach mission, missions do begin
        probes = call_function(mission+'_probes')
        foreach probe, probes do begin
            mission_probe = [mission,probe]
            mission_probe_str = mission+'_'+probe
            del_data, '*'
            var = micro_injection_read_cwt_kev_electron(time_range, mission_probe=mission_probe)
            print, 'Processed '+mission_probe_str+' for '+time_string(date)+' ...'
        endforeach
    endforeach
endforeach

stop


time_range = time_double(['2017-05-19','2017-05-20'])
time_range = ['2016-03-05/12:50','2016-03-06/02:40']
time_range = ['2016-11-01','2016-11-02']
mission_probe = ['mms','4']
;
;;routine = 'mms_read_kev_electron'
;;routine2 = routine+'_cdaweb'
;;probe = '4'
;;flux_var = call_function(routine, time_range, probe=probe, spec=1)
;;flux_var2 = call_function(routine2, time_range, probe=probe)
;;stop
;
var = micro_injection_read_cwt_kev_electron(time_range, mission_probe=mission_probe)
stop

time_range = time_double(['2017-05-19','2017-05-20'])
mission_probe = ['mms','4']
energy = 48e3   ; eV
data_file = join_path([homedir(),'test_micro_injection_cwt_kev_electron.cdf'])
file = micro_injection_read_cwt_kev_electron_gen_file_v02(time_range, mission_probe=mission_probe, energy=energy, filename=data_file)
end