;+
; Interface to read mlt, dtilt, bmod for certain probe and time.
;-
pro azim_dp_read_theta, time, probe=probe, datatype=datatype, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, _extra=ex

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','azim_dp'])
    if n_elements(version) eq 0 then version = 'v01'
    prefix = probe+'_'


;---Check if already have the data.
    theta_var = prefix+'theta'
    if ~check_if_update(theta_var, time) then return


;---Init settings.
    base_name = 'azim_dp_theta_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,'theta','%Y']
    ; MLT.
    request = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', prefix+'theta', $
                'time_var_name', 'bfield_ut', $
                'time_var_type', 'unix')))


;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)

    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            azim_dp_read_theta_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif


;---Check probe. If not checked, then run _gen_file to check.
    foreach file, files do begin
        foreach var_list, request.var_list do begin
            foreach in_var, var_list.in_vars do begin
                if cdf_has_var(in_var, filename=file) then continue
                the_time = cdf_read_var('dummy_ut', filename=file)
                azim_dp_read_theta_gen_file, the_time, probe=probe, filename=file
            endforeach
        endforeach
    endforeach


;---Read data from files and save to memory.
    foreach var_list, request.var_list do begin
        in_vars = var_list.in_vars
        foreach in_var, in_vars do begin
            xxs = []
            yys = []
            foreach file, files do begin
                cdf_load_var, in_var, filename=file
                get_data, in_var, txx, tyy, limits=lim
                xxs = [xxs, txx]
                yys = [yys, tyy]
            endforeach
            index = lazy_where(xxs, '[]', time, count=count)
            if count eq 0 then begin
                del_data, in_var
                errmsg = 'No data in given time ...'
            endif else begin
                xxs = xxs[index]
                yys = yys[index,*,*,*,*,*,*]
                store_data, in_var, xxs, yys, limits=lim
            endelse
        endforeach
    endforeach

end

time_range = time_double(['2019-08-01','2019-09-01'])
probes = [$
    'rbsp'+['a'], $
    'th'+['a','d','e'], $
    'g'+['14','15','16','17'], $
    'mms1']
    
    ; No data, need to be careful when interpolate.
    time_range = time_double(['2019-08-27','2019-08-28'])
    probes = ['mms1']

    time_range = time_double(['2019-09-06','2019-09-07'])
    probes = 'g'+['14','15','16','17']

foreach probe, probes do azim_dp_read_theta, time_range, probe=probe
end
