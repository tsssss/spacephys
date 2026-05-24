;+
; :Returns: double in [3,3], rotation matrix from instrument fixed to TSCS.
; :Arguments:
;   instr_str: in, required, string. 
; :Keywords:
;   probe: in, required, string.
;   errmsg: out, optional, string.
;   verbose: in, optional, boolean.
;-

function tracers_get_m_instr2tscs, instr_str, probe=probe, errmsg=errmsg, verbose=verbose
    compile_opt idl2

    errmsg = ''
    retval = !null

    if n_elements(instr_str) eq 0 then begin
        errmsg = 'Instrument name is required.'
        return, retval
    endif

    instr = strlowcase(instr_str[0])
    valid_instr = tracers_get_valid_instr()
    index = where(valid_instr eq instr, count)
    if count eq 0 then begin
        errmsg = 'Unknown TRACERS instrument: '+instr
        return, retval
    endif

    prefix = 'ts'+probe+'_'

    ; Root dir.
    tracers_local_root = tracers_get_local_root()
    spice_local_root = join_path([tracers_local_root,'SOC','spice'])
    if not file_test(spice_local_root,directory=1) then begin
        errmsg = 'No SPICE directory ...'
        return, retval
    endif

    probe_dir = 'ts'+probe
    probe_local_root = join_path([spice_local_root,probe_dir])

    ; Collect kernels to load.
    files = list()

    ; Spacecraft frame.
    probe_tf_file = join_path([probe_local_root,'frame',prefix+'frame_v1.0.tf'])
    files.add, probe_tf_file
    tf_file = join_path([spice_local_root,'geophys','ts_geophys_frame_v1.0.tf'])
    files.add, tf_file
  
    ; Load kernels, skip if already loaded.
    kernels = spice_get_loaded_kernel()
    foreach file, files do begin
        index = where(kernels eq file, count)
        if count eq 0 then begin
            if keyword_set(verbose) then lprmsg, 'Loading kernel '+file+' ...'
            cspice_furnsh, file
        endif else begin
            if keyword_set(verbose) then lprmsg, 'Kernel '+file+' already loaded ...'
        end
    endforeach

    ; Get the rotation matrix from instruement fixed to TSCS.
;    ut0 = '2025-09-30T00:00:00'
;    cspice_str2et, ut0, et0
    et0 = 0d
    target_coord = 'tscs'
    orig_coord = instr
    target_frame = strupcase(prefix+target_coord)
    orig_frame = strupcase(prefix+orig_coord)
    cspice_pxform, orig_frame, target_frame, et0, m_instr2tscs

;    ; Force double, and fix sqrt(2)/2.
;    m_instr2tscs = double(m_instr2tscs)
;    sqrt2_over_2 = sqrt(0.5d0)
;    index = where(abs(abs(m_instr2tscs)-sqrt2_over_2) lt 1d-6, count)
;    if count gt 0 then begin
;        stop
;        signs = m_instr2tscs[index]/abs(m_instr2tscs[index])
;        m_instr2tscs[index] = sqrt2_over_2*signs
;    endif

    return, m_instr2tscs

end

compile_opt idl2
valid_instr = tracers_get_valid_instr()
probe = '2'
foreach instr_str, valid_instr do begin
    print, 'Getting rotation matrix for '+instr_str+' ...'
    m_instr2tscs = tracers_get_m_instr2tscs(instr_str, probe=probe)
    msg = 'Rotation matrix from '+strupcase(instr_str)+' to TSCS:'
    print, msg
    print, m_instr2tscs
endforeach

end
