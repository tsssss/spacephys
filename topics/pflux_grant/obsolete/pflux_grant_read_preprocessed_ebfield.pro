;+
; Read preprocessed dE and dB.
; Save output data in rbspx_de_fac, rbspx_dedot0_fac, rbspx_db_fac.
; coord=. Default is 'fac', can be 'gsm'.
;-
pro pflux_grant_read_preprocessed_ebfield, time, probe=probe, id=datatype, coord=coord, $
    filename=data_file, errmsg=errmsg, local_root=local_root, project=project

    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])
    if n_elements(probe) ne 1 then message, 'Invalid probe ...'
    if n_elements(version) eq 0 then version = 'v01'
    index = where(['a','b'] eq probe, count)
    if count eq 0 then message, 'Invalid probe ...'
    rbspx = 'rbsp'+probe
    prefix = rbspx+'_'
    local_path = [local_root,rbspx,'ebfield','%Y']
    base_name = rbspx+'_ebfield_%Y_%m%d_'+version+'.cdf'
    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project.pflux_calc_setting
    valid_range = pflux_calc_setting[rbspx].time_range+[-1,1]*constant('secofday')

    if n_elements(datatype) eq 0 then datatype = 'ebfield'

    type_dispatch = hash()
    type_dispatch['b0'] = dictionary($
        'pattern', dictionary($
        'local_file', join_path([local_path,base_name]), $
        'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', valid_range, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', prefix+['b0_gsm'], $
                'time_var_name', 'ut_orbit', $
                'time_var_type', 'unix')))
    type_dispatch['ebfield'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', valid_range, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', prefix+['b1_gsm','e_uv','ew','bw_ratio'], $
                'time_var_name', 'ut_sec', $
                'time_var_type', 'unix'), $
            dictionary($
                'in_vars', prefix+['b0_gsm'], $
                'time_var_name', 'ut_orbit', $
                'time_var_type', 'unix')))
    type_dispatch['bw_ratio'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', valid_range, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', prefix+['bw_ratio'], $
                'time_var_name', 'ut_sec', $
                'time_var_type', 'unix')))

    request = type_dispatch[datatype]

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)


;---Read data from files and save to memory.
    read_files, time, files=files, request=request


;---Further settings.
    if datatype eq 'bw_ratio' then return

    xyz = constant('xyz')

    if datatype eq 'b0' then begin
        bvar = prefix+'b0_gsm'
        add_setting, bvar, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B', $
            coord: 'GSM', $
            coord_labels: xyz }
        return
    endif

    bvar = prefix+'b1_gsm'
    add_setting, bvar, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'dB', $
        coord: 'GSM', $
        coord_labels: xyz }

    e_uv_var = prefix+'e_uv'
    add_setting, e_uv_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E', $
        coord: 'UVW', $
        coord_labels: ['u','v'], $
        colors: sgcolor(['blue','red']) }


;---Remove bad data.
    bad_dates = ['2012-10-05','2015-03-16','2015-03-18','2015-03-19','2014-04-28','2014-04-29','2014-04-30']
    secofday = constant('secofday')
    fillval = !values.f_nan
    if probe eq 'a' then begin
        get_data, e_uv_var, times, e_uv
        foreach bad_date, bad_dates do begin
            index = lazy_where(times, '[]', bad_date+[0,secofday], count=count)
            if count ne 0 then e_uv[index,*] = fillval
        endforeach
    endif

;---Calculate E0.
    rbsp_read_quaternion, time, probe=probe
    get_data, e_uv_var, common_times, e_uv
    ncommon_time = n_elements(common_times)
    e_uvw = [[e_uv[*,0]],[e_uv[*,1]],[fltarr(ncommon_time)]]
    e_uvw_var = prefix+'e_uvw'
    store_data, e_uvw_var, common_times, e_uvw
    add_setting, e_uvw_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'dE', $
        coord: 'UVW', $
        coord_labels: ['u','v','w'] }
    q_var = prefix+'q_uvw2gsm'
    e_gsm_var = prefix+'e_gsm'
    rbsp_uvw2gsm, e_uvw_var, e_gsm_var, quaternion=q_var, probe=probe


;---Calculate Edot0.
    ew_var = prefix+'ew'
    e_uvw[*,2] = get_var_data(ew_var)
    bw_ratio_var = prefix+'bw_ratio'
    bw_ratio = get_var_data(bw_ratio_var, at=common_times)
    min_bw_ratio = 0.2  ; 10-15 deg.
    index = where(abs(bw_ratio) le min_bw_ratio, count)
    if count ne 0 then e_uvw[index,*] = !values.f_nan
    max_valid_e = 200.
    index = where(snorm(e_uvw) ge max_valid_e, count)
    if count ne 0 then e_uvw[index,*] = !values.f_nan
    edot0_uvw_var = prefix+'edot0_uvw'
    store_data, edot0_uvw_var, common_times, e_uvw
    add_setting, edot0_uvw_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'dE!S!Udot0!N!R', $
        coord: 'UVW', $
        coord_labels: ['u','v','w'] }
    edot0_gsm_var = prefix+'edot0_gsm'
    rbsp_uvw2gsm, edot0_uvw_var, edot0_gsm_var, quaternion=q_var, probe=probe

    if n_elements(coord) eq 0 then coord = 'fac'
    if coord eq 'gsm' then return

;---Convert to FAC.
    rbsp_read_orbit, time, probe=probe
    r_gsm_var = prefix+'r_gsm'
    common_time_step = total(common_times[0:1]*[-1,1])
    interp_time, r_gsm_var, common_times
    b0_gsm_var = prefix+'b0_gsm'
    interp_time, b0_gsm_var, common_times
    add_setting, b0_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'nT', $
        short_name: 'B0', $
        coord: 'GSM', $
        coord_labels: ['x','y','z'] }
    define_fac, b0_gsm_var, r_gsm_var
    db_gsm_var = prefix+'b1_gsm'
    in_vars = [edot0_gsm_var,e_gsm_var,db_gsm_var]
    out_vars = prefix+['dedot0_fac','de_fac','db_fac']
    foreach in_var, in_vars, ii do to_fac, in_var, to=out_vars[ii]

end

time_range = time_double(['2013-06-07','2013-06-08'])
time_range = time_double(['2013-05-01','2013-05-02'])
probe = 'b'
pflux_grant_read_preprocessed_ebfield, time_range, probe=probe
prefix = 'rbsp'+probe+'_'

de_var = prefix+'de_fac'
db_var = prefix+'db_fac'
pf_var = prefix+'pf_fac'
get_data, de_var, times, de
get_data, db_var, times, db
time_step = sdatarate(times)
width_small = 1/time_step
for ii=0,2 do begin
    de[*,ii] = smooth(de[*,ii], width_small, /nan, /edge_zero)
    db[*,ii] = smooth(db[*,ii], width_small, /nan, /edge_zero)
endfor
pf = spoynt(de,db)
store_data, pf_var, times, pf
add_setting, pf_var, /smart, {$
    display_type: 'vector', $
    unit: 'mW/m!U2!N', $
    short_name: 'S', $
    coord: 'FAC', $
    coord_labels: ['||','west','outward'] }

de_var = prefix+'dedot0_fac'
db_var = prefix+'db_fac'
pf_var = prefix+'pfdot0_fac'
get_data, de_var, times, de
get_data, db_var, times, db
time_step = sdatarate(times)
width_small = 1/time_step
for ii=0,2 do begin
    de[*,ii] = smooth(de[*,ii], width_small, /nan, /edge_zero)
    db[*,ii] = smooth(db[*,ii], width_small, /nan, /edge_zero)
endfor
pf = spoynt(de,db)
store_data, pf_var, times, pf
add_setting, pf_var, /smart, {$
    display_type: 'vector', $
    unit: 'mW/m!U2!N', $
    short_name: 'S', $
    coord: 'FAC', $
    coord_labels: ['||','west','outward'] }
end
