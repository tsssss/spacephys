;+
; 
;-
function test_tracers_cotrans, test=test
    compile_opt idl2

;---Settings.
    probe = '2'
    time_range = time_double(['2025-09-30/09:30','2025-09-30/10:00'])
    prefix = 'ts'+probe+'_'

;---Load kernels.
    tracers_load_spice_kernel, time_range, probe=probe
    kernels = spice_get_loaded_kernel()

    time_step = 1
    times = make_bins(time_range, time_step)
    ut0 = time_string(times[0],tformat='YYYY-MM-DDThh:mm:ss')
    cspice_str2et, ut0, et0
    ets = et0+times-ut0
    uts = time_string(times,tformat='YYYY-MM-DDThh:mm:ss.ffffff')
    cspice_str2et, uts, ets
    ntime = n_elements(times)
    ndim = 3
    default_coord = prefix+'tscs'
    target_frame = strupcase(default_coord)

;---Get position.
    target = 'TS'+probe
    observer = 'EARTH'
    abcoor = 'NONE'
    cspice_spkezr, target, ets, target_frame, abcoor, observer, state, local_time
    r_gse = transpose(state[0:2,*])
    re = constant('re')
    settings = dictionary('coord',default_coord, 'mission_probe', 'ts'+probe, 'requested_time_range', time_range)
    r_var_spice = var_store(prefix+'r_'+default_coord, r_gse/re, times, id='position', settings=settings)
  

    valid_instrs = tracers_get_valid_instr()
    foreach instr, valid_instrs do begin
        orig_frame = target+'_'+strupcase(instr)
        cspice_pxform, orig_frame, target_frame, ets[0], pxform
        xyz_orig = dblarr([ndim,ndim])
        xyz_orig[*,0] = [1,0,0d]
        xyz_orig[*,1] = [0,1,0d]
        xyz_orig[*,2] = [0,0,1d]
        m_instr2tscs = pxform[*,*,0]
        msg = 'Rotation matrix from '+orig_frame+' to '+target_frame+':'
        print, msg
        print, m_instr2tscs

        errmsg = ''
        func_name = 'ct_ts_'+instr+'2ts_tscs'
        vec_tscs = call_function(func_name, xyz_orig, probe=probe, errmsg=errmsg)
        if errmsg ne '' then begin
            print, func_name+' returned error: '+errmsg
            return, 0
        endif

        vec_tscs_expected = rotate_vector(xyz_orig, m_instr2tscs)
        if max(abs(vec_tscs-vec_tscs_expected)) gt 1e-6 then begin
            print, func_name+' returned unexpected TSCS vectors.'
            return, 0
        endif

        errmsg = ''
        vec_tscs = tracers_instr2ts_tscs(xyz_orig, instr_str=instr, probe=probe, errmsg=errmsg)
        if errmsg ne '' then begin
            print, 'tracers_instr2ts_tscs returned error: '+errmsg
            return, 0
        endif

        if max(abs(vec_tscs-vec_tscs_expected)) gt 1e-6 then begin
            print, 'tracers_instr2ts_tscs returned unexpected TSCS vectors for '+instr+'.'
            return, 0
        endif

        xyz_tscs = dblarr([ndim,ndim])
        xyz_tscs[*,0] = [1,0,0d]
        xyz_tscs[*,1] = [0,1,0d]
        xyz_tscs[*,2] = [0,0,1d]
        m_tscs2instr = transpose(m_instr2tscs)
        vec_instr_expected = rotate_vector(xyz_tscs, m_tscs2instr)

        errmsg = ''
        func_name = 'ct_ts_tscs2ts_'+instr
        vec_instr = call_function(func_name, xyz_tscs, probe=probe, errmsg=errmsg)
        if errmsg ne '' then begin
            print, func_name+' returned error: '+errmsg
            return, 0
        endif

        if max(abs(vec_instr-vec_instr_expected)) gt 1e-6 then begin
            print, func_name+' returned unexpected instrument vectors.'
            return, 0
        endif

        errmsg = ''
        vec_instr = tracers_ts_tscs2instr(xyz_tscs, instr_str=instr, probe=probe, errmsg=errmsg)
        if errmsg ne '' then begin
            print, 'tracers_ts_tscs2instr returned error: '+errmsg
            return, 0
        endif

        if max(abs(vec_instr-vec_instr_expected)) gt 1e-6 then begin
            print, 'tracers_ts_tscs2instr returned unexpected instrument vectors for '+instr+'.'
            return, 0
        endif
    endforeach

    return, 1
end


compile_opt idl2
test = 1

probe = '2'
time_range = time_double(['2025-09-30/09:30','2025-09-30/09:40'])
time_step = 60

mission = 'ts'
prefix = 'ts'+probe+'_'
times = make_bins(time_range, time_step)
ntime = n_elements(times)
ndim = 3
vec_mag = fltarr(ntime,ndim)
vec_mag[*,0] = smkarthm(-1,1,ntime,'n')
vec_mag[*,1] = smkarthm(-0.5,0.5,ntime,'n')
vec_mag[*,2] = 0.2
mag_var = prefix+'vec_mag'
settings = dictionary($
    'coord', 'ts_mag', $
    'mission_probe', 'ts'+probe, $
    'requested_time_range', time_range )
mag_var = var_store(mag_var, vec_mag, times, settings=settings)
vec_magic = cotran_pro(vec_mag, times, coord_msg=['ts_mag','ts_magic'], probe=probe, mission=mission, errmsg=errmsg)
settings['coord'] = 'ts_magic'
magic_var = prefix+'vec_magic'
magic_var = var_store(magic_var, vec_magic, times, settings=settings)
plot_vars = [mag_var, magic_var]
sgopen, 0, size=[6,6]
tplot, plot_vars, trange=time_range
; print, test_tracers_cotrans(test=test)
end
