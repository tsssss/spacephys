;+
; Training data:
;   5-min cadence; yearly CDF; common times; contain data for all energies, all species; all needed features (spacecraft position, omni data)
; I remove NaN in all data unless it is impossible.
; I do this step in IDL because it's easier to check the data quality in IDL using tplot than in Python.
;-

function ml_hope_load_training_data, input_time_range, id=datatype, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, resolution=resolution

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 86400d*120
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','ml_hope'])
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(resolution) eq 0 then resolution = ml_time_step()
    resolution_str = string(round(resolution/60),format='(I0)')+'min'

    if size(input_time_range[0],type=1) eq 7 then begin
        time_range = time_double(input_time_range)
    endif else begin
        time_range = input_time_range
    endelse

;---Init settings.
    valid_range = ['2012','2020']
    base_name = 'ml_hope_training_data_'+resolution_str+'_%Y_'+version+'.cdf'
    local_path = [local_root,'training_data_'+resolution_str]
    request = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file(sync=1)])), $
        'valid_range', time_double(valid_range), $
        'sync_threshold', sync_threshold, $
        'cadence', 'year', $
        'extension', fgetext(base_name) )


;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time_range, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            ml_hope_load_training_data_gen_file, file_time, filename=local_file, resolution=resolution
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time_range, nonexist_files=nonexist_files)
    endif
    
    if n_elements(files) eq 0 then return, '' else return, files

end

years = make_bins([2012,2019],1)
year_strs = string(years,format='(I0)')
foreach year, years do begin
    time_range = string(year+[0,1],format='(I4)')
    file = ml_hope_load_training_data(time_range)
    print, file
endforeach
end