;+
; :Purpose: Load SPICE kernels.
; :Arguments:
;   input_time_range: in, required. The time range.
; :Keywords:
;   probe: in, required. The probe.
;   orbit_type: in, optional. 'def' or 'predict'
;   errmsg: out, optional. Error message.
;-

pro tracers_load_spice_kernel, input_time_range, probe=probe, orbit_type=orbit_type, errmsg=errmsg
    compile_opt idl2
    errmsg = ''

    if n_elements(input_time_range) eq 0 then begin
        errmsg = 'No input time range ...'
        return
    end
    time_range = time_double(input_time_range)

    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        return
    end

    if not is_valid_probe(tracers_get_probes(),probe) then begin
        errmsg = 'Invalid probe ...'
        return
    endif

    ; Root dir.
    tracers_local_root = tracers_get_local_root()
    spice_local_root = join_path([tracers_local_root,'SOC','spice'])
    if not file_test(spice_local_root,directory=1) then begin
        errmsg = 'No SPICE directory ...'
        return
    endif
    naif_local_root = join_path([spice_local_root,'naif_kernels'])
    probe_dir = 'ts'+probe
    probe_local_root = join_path([spice_local_root,probe_dir])
    prefix = 'ts'+probe+'_'

    ; Collect kernels to load.
    files = list()

    ; Leap second.
    lsk_file = join_path([naif_local_root,'naif0012.tls'])
    lprmsg, 'Adding leap second kernel ...'
    files.add, lsk_file

    ; Earth's constant.
    pck_file = join_path([naif_local_root,'pck00011.tpc'])
    lprmsg, 'Adding planetary constants ...'
    files.add, pck_file

    ; SC clock.
    base = prefix+'sclk_v1.0.tf'
    sclk_file = join_path([probe_local_root,'sclk',base])
    lprmsg, 'Adding SC clock ...'
    files.add, sclk_file

    ; Planetary position and velocity.
    base = 'de430.bsp'
    bsp_file = join_path([naif_local_root,base])
    lprmsg, 'Adding planetary position and velocity ...'
    files.add, bsp_file

    ; Earth's positions and velocities.
    avail_files = file_search(join_path([naif_local_root,'*.bpc']))
    pattern = 'earth_000101_([0-9]+)_[0-9]+.bpc'
    tmp = stregex(avail_files, pattern, subexpr=1, extract=1)
    time_strs = reform(tmp[1,*])
    avail_files = reform(tmp[0,*])
    file_times = time_double(time_strs,tformat='hhmmdd')
    tmp = max(file_times, index)
    base = avail_files[index]
    bpc_file = join_path([naif_local_root,base])
    lprmsg, 'Adding Earth position and velocity ...'
    files.add, bpc_file
    

    ; Spacecraft frame.
    probe_tf_file = join_path([probe_local_root,'frame',prefix+'frame_v1.0.tf'])
    files.add, probe_tf_file
    tf_file = join_path([spice_local_root,'geophys','ts_geophys_frame_v1.0.tf'])
    files.add, tf_file

    ; Predict kernels for TSCS and TSB coordinates.
    coords = ['tscs','tsb']
    foreach coord, coords do begin
        kernel_path = join_path([probe_local_root,'predict',coord])
        pattern = prefix+coord+'_ck_fake-gnc_([0-9]{4}-[0-9]{2}-[0-9]{2})_v[0-9]{2}.bc'
        avail_files = file_search(join_path([kernel_path,'*.bc']))
        tmp = stregex(avail_files, pattern, subexpr=1, extract=1)
        time_strs = reform(tmp[1,*])
        avail_bases = reform(tmp[0,*])
        file_times = time_double(time_strs,tformat='YYYY-MM-DD')
        ; Sort in time.
        index = sort(file_times)
        file_times = file_times[index]
        avail_bases = avail_bases[index]
        index = where_pro(file_times, '>', time_range[1], count=count)
        if count eq 0 then begin
            errmsg = 'No position and attitude data found in given time_range ...'
            return
        endif
        i1 = index[0]-1
        index = where_pro(file_times, '<', time_range[0], count=count)
        if count eq 0 then begin
            errmsg = 'No position and attitude data found in given time_range ...'
            return
        endif
        i0 = index[count-1]
        wanted_bases = avail_bases[i0:i1]
        foreach base, wanted_bases do begin
            files.add, join_path([kernel_path,base])
        endforeach
    endforeach


    ; Spacecraft position and attitude.
    if n_elements(orbit_type) eq 0 then orbit_type = 'predict'
    bsp_root = join_path([probe_local_root,orbit_type,'spk'])
    avail_files = file_search(join_path([bsp_root,prefix+'peph_tle_*.bsp']))
    pattern = prefix+'peph_tle_([0-9]{4}-[0-9]{2}-[0-9]{2})_v[0-9]{2}.bsp'
    tmp = stregex(avail_files, pattern, subexpr=1, extract=1)
    time_strs = reform(tmp[1,*])
    avail_bases = reform(tmp[0,*])
    file_times = time_double(time_strs,tformat='YYYY-MM-DD')
    ; Sort in time.
    index = sort(file_times)
    file_times = file_times[index]
    avail_bases = avail_bases[index]
    index = where_pro(file_times, '>', time_range[1], count=count)
    if count eq 0 then begin
        errmsg = 'No position and attitude data found in given time_range ...'
        return
    endif
    i1 = index[0]-1
    index = where_pro(file_times, '<', time_range[0], count=count)
    if count eq 0 then begin
        errmsg = 'No position and attitude data found in given time_range ...'
        return
    endif
    i0 = index[count-1]
    wanted_bases = avail_bases[i0:i1]
    foreach base, wanted_bases do begin
        files.add, join_path([bsp_root,base])
    endforeach

  
    ; Load kernels, skip if already loaded.
    kernels = spice_get_loaded_kernel()
    foreach file, files do begin
        index = where(kernels eq file, count)
        if count eq 0 then begin
            lprmsg, 'Loading kernel '+file+' ...'
            cspice_furnsh, file
        endif else begin
            lprmsg, 'Kernel '+file+' already loaded ...'
        end
    endforeach
    


end


compile_opt idl2
probe = '1'
time_range = time_double(['2025-12-11','2025-12-12'])

prefix = 'ts'+probe+'_'
tracers_load_spice_kernel, time_range, probe=probe


time_step = 10.
times = make_bins(time_range, time_step)
ut0 = time_string(times[0],tformat='YYYY-MM-DDThh:mm:ss')
cspice_str2et, ut0, et0
ets = et0+times-ut0
uts = time_string(times,tformat='YYYY-MM-DDThh:mm:ss')
cspice_str2et, uts, ets
target = 'TS'+probe
observer = 'EARTH'
frame = 'GEO'
abcoor = 'NONE'
cspice_spkezr, target, ets, frame, abcoor, observer, state, local_time
r_geo = transpose(state[0:2,*])
v_geo = transpose(state[3:5,*])
re = constant('re')
settings = dictionary('coord','geo')
r_var = var_store(prefix+'r_geo', r_geo/re, times, id='position', settings=settings)
v_var = var_store(prefix+'v_geo', v_geo, times, id='velocity', settings=settings)
plot_vars = [r_var,v_var]
tplot, plot_vars, trange=time_range
stop

end