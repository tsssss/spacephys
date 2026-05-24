;+
; :Purpose: Initializes SPICE kernels for tracers.
; :Returns: Success or not.
; :Keywords:
;   probe: in, optional, load generic kernels by default.
;   local_root: in, optional, location to save the data.
;   remote_root: in, optional, url of the data.
;   update: in, optional, set to update kernels.
;-
function tracers_init_spice, probe=probe, $
    local_root=local_root, remote_root=remote_root, update=update
    compile_opt idl2
    on_error, 0

    errmsg = ''
    retval = 0

;---Check inputs.
    if n_elements(probe) eq 0 then probe = []
    if n_elements(local_root) eq 0 then local_root = tracers_get_local_root() 
    if n_elements(remote_root) eq 0 then remote_root = tracers_get_remote_root()

    spice_remote_path = join_path([remote_root,'SOC','spice'])
    spice_local_path = join_path([local_root,'SOC','spice'])

    username = tracers_get_remote_root('username')
    password = tracers_get_remote_root('password')

;---Download parent dir.
    spice_root_index = join_path([spice_local_path,'index.html'])
    if file_test(spice_root_index) eq 0 then begin
        download_file, spice_root_index, spice_remote_path, username=username, password=password
        if not file_test(spice_root_index) then begin
            errmsg = 'Please download manually from '+spice_remote_path+' ...'
            return, retval
        endif
    endif

    ; Read subdirectory names from index file.
    subdirs = download_parse_subdirs(spice_root_index)

    ; Collect all download dirs.
    download_dirs = list()
    if n_elements(probe) eq 0 then begin
        probe_dirs = 'ts'+tracers_get_probes()+'/'
        foreach dir, subdirs do begin
            index = where(probe_dirs eq dir, count)
            if count ne 0 then continue
            download_dirs.add, dir
        endforeach
    endif else begin
        valid_probes = tracers_get_probes()
        index = where(valid_probes eq probe, count)
        if count ne 0 then download_dirs.add, 'ts'+probe+'/'
    endelse
    

    ; Download probe and generic dirs recursively.
    foreach dir, download_dirs do begin
        remote_path = join_path([spice_remote_path,dir])
        local_path = join_path([spice_local_path,dir])
        tmp = download_recursively(local_path, remote_path, $
            username=username, password=password, update=update)
    endforeach

    return, 1

end


compile_opt idl2
foreach probe, [!null,tracers_get_probes()] do begin
    print, tracers_init_spice(probe=probe)
endforeach
end