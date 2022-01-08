;+
; Read Polar hydra ion and ele energy flux.
; Need to generate the files using polar_gen_iowa_hydra_moments.pro first.
;-

pro polar_read_eflux, time, errmsg=errmsg, probe=probe

    prefix = 'po_'
    local_root = join_path([default_local_root(),'sdata','polar','hydra'])
    version = 'v01'

    valid_range = time_double(['1996-01-01','2009-01-01'])
    base_name = 'po_hyd_moments_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,'%Y']

    vars = ['ion','ele']+'_energy_flux'
    request = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', vars, $
                'out_vars', 'po_'+vars, $
                'time_var_name', 'ut_sec', $
                'time_var_type', 'unix')))

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)
    nfile = n_elements(files)
    if nfile ne 0 then begin
        empty_file_flags = bytarr(nfile)
        foreach file, files, ii do begin
            empty_file = 0
            foreach var, vars do if ~cdf_has_var(var, filename=file) then empty_file = 1
            empty_file_flags[ii] = empty_file
        endforeach
        index = where(empty_file_flags eq 0, count)
        if count eq 0 then files = !null else files = files[index]
    endif
    nfile = n_elements(files)
    if nfile eq 0 then begin
        errmsg = 'No data file'
        return
    endif


;---Read data from files and save to memory.
    read_files, time, files=files, request=request


    types = ['ion','ele']
    vars = 'po_'+types+'_energy_flux'
    fillval = !values.f_nan
    foreach var, vars, ii do begin
        get_data, var, times, data
        if min(data) eq max(data) then begin
            errmsg = 'No valid data'
            return
        endif
        index = where(data eq -1e31, count)
        if count ne 0 then begin
            data[index] = fillval
            store_data, var, times, data
        endif
        add_setting, var, /smart, dictionary($
            'display_type', 'scalar', $
            'unit', '(mW/m!U2!N)', $
            'short_name', 'In-situ '+tex2str('Gamma')+'!D'+types[ii] )
    endforeach
end
