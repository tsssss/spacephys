;+
; :Purpose: Get ROI time range that encloses the input time range.
;-
function tracers_get_roi_time_range, input_time_range, probe=probe, errmsg=errmsg, $
    local_root=local_root, remote_root=remote_root, update=update
    compile_opt idl2
    on_error, 0

    errmsg = ''
    retval = !null

;---Check inputs.
    if check_valid_probe(probe, tracers_get_probes()) eq 0 then return, retval
    if n_elements(local_root) eq 0 then local_root = tracers_get_local_root() 
    if n_elements(remote_root) eq 0 then remote_root = tracers_get_remote_root()

    probe_str = strupcase('ts'+probe)
    prefix = 'ts'+probe+'_'
    base = prefix+'roi-list.csv'
    roi_remote_file = join_path([remote_root,'SOC',probe_str,'events','roi_intervals',base])
    roi_local_file = join_path([local_root,'SOC',probe_str,'events','roi_intervals',base])
    update_file = 0
    if file_test(roi_local_file) eq 0 then update_file = 1
    if keyword_set(update) then update_file = 1
    if update_file then begin
        username = tracers_get_remote_root('username')
        password = tracers_get_remote_root('password')
        download_file, roi_local_file, roi_remote_file, $
            username=username, password=password
    endif

  
;---Read file and find time range.
    if file_test(roi_local_file) eq 0 then begin
        errmsg = 'Problem in downloading ROI file ...'
        return, retval
    endif
    lines = read_all_lines(roi_local_file)
    nheader = 6
    roi_lines = lines[nheader:*]
    nroi = n_elements(roi_lines)
    roi_trs = strarr(nroi,2)
    for ii=0,nroi-1 do begin
        infos = strsplit(roi_lines[ii],',',extract=1)
        roi_trs[ii,*] = infos[0:1]
    endfor
    roi_trs = time_double(roi_trs,tformat='YYYY-MM-DDThh:mm:ss')
    input_tr = time_double(input_time_range)
    if n_elements(input_time_range) eq 1 then input_tr = [0,0]+input_tr[0]
    index = where(roi_trs[*,0] le input_tr[0], count)
    if count eq 0 then begin
        errmsg = 'No ROI before requested time range ...'
        return, retval
    endif
    i0 = index[count-1]
    index = where(roi_trs[*,1] ge input_tr[1], count)
    if count eq 0 then begin
        errmsg = 'No ROI after requested time range ...'
        return, retval
    endif
    i1 = index[0]
    wanted_trs = roi_trs[i0:i1,*]
    return, wanted_trs

end

compile_opt idl2
tr = ['2025-11-22','2025-11-23']
tr = ['2025-10-11/00:51:09']
trs = tracers_get_roi_time_range(tr, probe='1')
if n_elements(trs) eq 0 then stop
print, time_string(trs)
end
