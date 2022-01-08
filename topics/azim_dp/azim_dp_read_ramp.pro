;+
; Read ramp for given time.
;-

function azim_dp_read_ramp, time, probe=probe, mlt_range=mlt_range, errmsg=errmsg

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','azim_dp'])
    if n_elements(version) eq 0 then version = 'v01'
    prefix = probe+'_'


;---Init settings.
    base_name = 'azim_dp_ramp_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,'ramp','%Y']
    type_dispatch = dictionary()
    request = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', prefix+'ramp_width', $
                'time_var_name', prefix+'ramp_time', $
                'time_var_type', 'unix')))


;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)

    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            azim_dp_read_ramp_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif


;---Check probe. If not checked, then run _gen_file to check.
    foreach file, files do begin
        in_var = prefix+'ramp_width'
        if cdf_has_var(in_var, filename=file) then continue
        the_time = cdf_read_var('dummy_ut', filename=file)
        azim_dp_read_ramp_gen_file, the_time, probe=probe, filename=file
    endforeach


;---Read data from files and save to memory.
    vars = ['time','width','height','time_range','theta_range','r_sm','mlt']
    ramp_info = dictionary()
    foreach var, vars do ramp_info[var] = []
    foreach file, files do begin
        foreach var, vars do begin
            the_var = prefix+'ramp_'+var
            ramp_info[var] = [ramp_info[var], cdf_read_var(the_var, filename=file)]
        endforeach
    endforeach


;---Convert to list of dict.
    nramp = n_elements(ramp_info['time'])
    ramp_list = list()
    for ramp_id=0,nramp-1 do begin
        the_ramp = dictionary()
        foreach var, vars do begin
            val = (ramp_info[var])[ramp_id,*,*]
            val = val[*]
            if n_elements(val) eq 1 then val = val[0]
            the_ramp[var] = val
        endforeach
        ramp_list.add, the_ramp
    endfor



    
;---Remove dummy ramps.
    retval = list()
    widths = fltarr(nramp)
    foreach ramp, ramp_list, ramp_id do widths[ramp_id] = ramp.width
    index = where(finite(widths), nramp)
    if nramp eq 0 then begin
        return, retval
    endif else begin
        ramp_list = ramp_list[index]
    endelse
    
    ; Filter time.
    times = dblarr(nramp)
    foreach ramp, ramp_list, ramp_id do times[ramp_id] = ramp.time
    index = lazy_where(times, '[]', time, count=nramp)
    if nramp eq 0 then return, list()
    ramp_list = ramp_list[index]
    
    ; Filter MLT.
    if n_elements(mlt_range) eq 2 then begin
        mlts = fltarr(nramp)
        foreach ramp, ramp_list, ramp_id do begin
            mlts[ramp_id] = ramp.mlt
        endforeach
        index = lazy_where(mlts, '[]', mlt_range, count=nramp)
        if nramp eq 0 then return, retval
        ramp_list = ramp_list[index]
    endif
    
    foreach ramp, ramp_list do ramp['probe'] = probe
    
    return, ramp_list


end


secofday = constant('secofday')
time_range = time_double('2014-08-28')+[0,secofday]
time_range = time_double(['2019-08-04','2019-08-06'])
probes = ['tha','thd','the','rbspa','rbspb','g14','g15','g16','g17','mms1']
foreach probe, probes do begin
    ramp = azim_dp_read_ramp(time_range, probe=probe)
endforeach
end
