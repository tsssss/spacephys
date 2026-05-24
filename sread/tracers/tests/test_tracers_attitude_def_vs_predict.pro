;+
; Test the difference between the predicted attitude vs definitive.
;-

; Test spin-axis pointing, c.f. rbsp_load_sa_pointing.pro.
function test_tracers_attitude_def_vs_predict, test=test
    compile_opt idl2

;---Settings.
    probe = '2'
    time_range = time_double(['2026-01-01','2026-01-01/00:20'])
    time_range = time_double(['2025-12-01','2025-12-01/00:20'])
    ; These are kernels manually checked for the above time range.
    def_kernel_file = '/Volumes/data/tracers/SOC/spice/ts2/def/tscs/hybrid_2026-02-15_v01/ts2_dtscs_2026-01-01_v0.1.0.bc'
    predict_kernel_file = '/Volumes/data/tracers/SOC/spice/ts2/predict/tscs/ts2_tscs_ck_fake-gnc_2025-12-15_v01.bc'

    probe = '2'
    time_range = time_double(['2025-12-22/15:25','2025-12-22/15:30'])
    prefix = 'ts'+probe+'_'

;---Load kernels.
    tracers_load_spice_kernel, time_range, probe=probe
    kernels = spice_get_loaded_kernel()

    time_step = 1
    times = make_bins(time_range, time_step)
    ut0 = time_string(times[0],tformat='YYYY-MM-DDThh:mm:ss')
    cspice_str2et, ut0, et0
    ets = et0+times-ut0
    uts = time_string(times,tformat='YYYY-MM-DDThh:mm:ss')
    cspice_str2et, uts, ets
    ntime = n_elements(times)
    ndim = 3
    default_coord = 'gse'
    target_frame = strupcase(default_coord)

;---Get position and B model.
    target = 'TS'+probe
    observer = 'EARTH'
    abcoor = 'NONE'
    cspice_spkezr, target, ets, target_frame, abcoor, observer, state, local_time
    r_gse = transpose(state[0:2,*])
    re = constant('re')
    settings = dictionary('coord',default_coord, 'mission_probe', 'ts'+probe, 'requested_time_range', time_range)
    r_var_spice = var_store(prefix+'r_'+default_coord, r_gse/re, times, id='position', settings=settings)

    bmod_var = lets_read_geopack_bfield(orbit_var=r_var_spice, internal_model='igrf', external_model='t89')

    

;---Get the rotation matrix from instruement fixed to GSE.
    foreach instr, ['mag','efi'] do begin
        orig_frame = target+'_'+strupcase(instr)
        cspice_pxform, orig_frame, target_frame, ets, pxform
        xyz_orig = dblarr([ndim,ndim])
        xyz_orig[*,0] = [1,0,0d]
        xyz_orig[*,1] = [0,1,0d]
        xyz_orig[*,2] = [0,0,1d]
        xyz_coord = dblarr(ntime,ndim,ndim)
        for ii=0,ntime-1 do begin
            for jj=0,ndim-1 do begin
                xyz_coord[ii,*,jj] = pxform[*,*,ii] ## xyz_orig[*,jj]
            endfor
        endfor

        components = constant('xyz')
        for ii=0,ndim-1 do begin
            var = prefix+instr+components[ii]+'_'+default_coord
            settings = dictionary('coord',default_coord)
            var = var_store(var, xyz_coord[*,*,ii], times, settings=settings)
            options, var, yrange=[-1,1]*1.2, constant=[0]
        endfor


    ;---Quaternion from TSCS to GSE.
        m_orig2gse = dblarr(ntime,ndim,ndim)
        for ii=0,ntime-1 do begin
            m_orig2gse[ii,*,*] = transpose(reform(pxform[*,*,ii]))
        endfor
        q_orig2gse = mtoq(m_orig2gse)
        orig_coord = 'ts_'+instr
        q_var = prefix+'q_'+orig_coord+'2'+default_coord
        settings = dictionary('coord', orig_coord+'2'+default_coord)
        q_var = var_store(q_var, q_orig2gse, times, settings=settings, id='quaternion')
    endforeach
    


;---Load B field and E field in the rotating frame (TSCS).
    e_var = tracers_read_efield(time_range, probe=probe)
    b_var = tracers_read_bfield(time_range, probe=probe)
    vars = [b_var, e_var]
    options, vars, mission_probe='ts'+probe, requested_time_range=time_range
    q_vars = prefix+'q_ts_'+['mag','efi']+'2'+default_coord
    new_vars = prefix+['b','e']+'_gsm'
    foreach var, vars, vid do begin
        vec_orig = var_get_data(var, times=times, settings=settings, in=time_range)
;        vec_orig[*,2] = 0
        q_var = q_vars[vid]
        q_orig2gse = var_get_data(q_var, times=uts)
        q_orig2gse = qslerp(q_orig2gse, uts, times)
        m_orig2gse = qtom(q_orig2gse)
        vec_gse = rotate_vector(vec_orig, m_orig2gse)
        vec_gsm = cotran_pro(vec_gse, times, coord_msg=['gse','gsm'])
        settings['coord'] = 'gsm'
        var = var_store(new_vars[vid], vec_gsm, times, settings=settings)
        vec_mag = snorm(vec_gsm)
        var = var_store(new_vars[vid]+'_norm', vec_mag, times)
        var = var_store(new_vars[vid]+'_z', vec_gsm[*,2], times)
    endforeach
    
;---dB.
;    b_gsm = var_get_data(prefix+'b_gsm', times=times, settings=settings)
;    bmod_gsm = var_get_data(prefix+'b_gsm_igrf_t89', at=times)
;    db_gsm = b_gsm - bmod_gsm
;    settings['short_name'] = 'dB'
;    db_var = var_store(prefix+'db_gsm', db_gsm, times, settings=settings)
    b_vars = lets_decompose_bfield(b_var=prefix+'b_gsm', bmod_var=bmod_var, b0_window=60d, update=1)
    b0_var = b_vars['b0']
    db_var = b_vars['b1']

    e_gsm = var_get_data(prefix+'e_gsm', times=times, settings=settings)
    dr = sdatarate(times)
    window = 120d
    width = window/dr
    for ii=0,ndim-1 do begin
        e_gsm[*,ii] -= smooth(e_gsm[*,ii], width, edge_mirror=1)
    endfor
    settings['short_name'] = 'dE'
    de_var = var_store(prefix+'de_gsm', e_gsm, times, settings=settings)
    options, de_var, constant=[0]


;---Plot.
    fig_size = [12,8]
    plot_dir = srootdir()
    version_str = 'v01'
    base = 'test_tracers_attitude_def_vs_predict_'+time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_'+version_str+'.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0
    sgopen, plot_file, size=fig_size

    plot_vars = prefix+['b_ts_mag','e_ts_efi','b_gsm','e_gsm','b_gsm_igrf_t89']
    plot_vars = [plot_vars,b0_var,db_var]
    plot_vars = prefix+['b_ts_mag','e_ts_efi','e_gsm_norm','de_gsm','b_gsm_norm','b_gsm']
;    plot_vars = prefix+['b_ts_mag','e_ts_efi','e_gsm_z','e_gsm_norm','e_gsm','b_gsm_z','b_gsm_norm','b_gsm']
    tplot, plot_vars, trange=time_range

    if keyword_set(test) then stop
    return, plot_file
end

compile_opt idl2
test = 1
print, test_tracers_attitude_def_vs_predict(test=test)
end