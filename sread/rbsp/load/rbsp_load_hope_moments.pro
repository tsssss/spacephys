;+
; Read HOPE moments.
;
; input_time_range. A time range in unix time or string. Set time to find files
;   automatically, or set files to read data in them directly.
; id=. A string sets the data type to read. Check supported ids by setting
;   print_datatype.
; print_datatype=. A boolean. Set to print all supported ids.
; probe=. A string set the probe to read data for.
; local_root=. A string to set the local root directory.
; local_files=. A string or an array of N full file names. Set to fine
;   tuning the files to read data from.
; file_times=. An array of N times. Set to fine tuning the times of the files.
;-

function rbsp_load_hope_moments, input_time_range, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root

    compile_opt idl2
    on_error, 0
    errmsg = ''


;---Check inputs.
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])
    if n_elements(remote_root) eq 0 then remote_root = 'https://cdaweb.gsfc.nasa.gov/pub/data/rbsp'
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(release) eq 0 then release = 'rel04'  ; updated 2019-06.
    if n_elements(datatype) eq 0 then datatype = 'l2%moments'

    if size(input_time_range[0],type=1) eq 7 then begin
        time_range = time_double(input_time_range)
    endif else begin
        time_range = input_time_range
    endelse

;---Init settings.
    type_dispatch = hash()
    ; Moment file, same as HOPE L2 file.
    case probe of
        'a': valid_range = ['2012-09-07','2019-10-14/24:00']
        'b': valid_range = ['2012-09-06','2019-07-16/24:00']
    endcase
    base_name = 'rbsp'+probe+'_hope_moments_'+release+'_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,'rbsp'+probe,'hope_moments_'+release,'%Y']
    type_dispatch['l2%moments'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)]) ), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name) )

    if keyword_set(print_datatype) then begin
        print, 'Suported data type: '
        ids = type_dispatch.keys()
        foreach id, ids do print, '  * '+id
        return, ids
    endif


;---Dispatch patterns.
    if n_elements(datatype) eq 0 then begin
        errmsg = handle_error('No input datatype ...')
        return, ''
    endif
    if not type_dispatch.haskey(datatype) then begin
        errmsg = handle_error('Do not support type '+datatype+' yet ...')
        return, ''
    endif
    request = type_dispatch[datatype]

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time_range, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            rbsp_load_hope_moments_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

    if n_elements(files) eq 0 then return, '' else return, files

end


time_range = ['2013-06-07','2013-06-09']
probe = 'a'
files = rbsp_load_hope_moments(time_range, probe=probe)
end
