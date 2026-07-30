
function tracers_read_bfield_l1, input_time_range, probe=probe, $
    update=update, get_name=get_name, errmsg=errmsg, _extra=extra
    compile_opt idl2

    errmsg = ''
    retval = !null
    prefix = 'ts'+probe+'_'

    default_coord = 'ts_mag'
    var_info = prefix+'b_'+default_coord
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then del_data, var_info
    time_range = time_double(input_time_range)
    if ~check_if_update(var_info, time_range) then return, var_info

    ;files = tracers_load_mag(time_range, probe=probe, errmsg=errmsg, id='l1b%bdc-bor_x233')
    files = tracers_load_mag(time_range, probe=probe, errmsg=errmsg, id='l1b%bdc-roi_x232')
    if errmsg ne '' then return, retval

    var_list = list()
    ;in_vars = prefix+'l1b_bdc_bor'
    in_vars = prefix+'l1b_bdc_roi'
    time_var = 'Epoch'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2000' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    b_raws = var_get_data(in_vars, times=uts)
    dts = cdf_read_var('EpochOffset', filename=files[0])*1e-9   ; ns to sec.
    ndt = n_elements(dts)
    nut = n_elements(uts)
    ntime = nut*ndt
    ndim = 3
    times = dblarr(ndt,nut)
    b_vecs = dblarr(ndt,nut,ndim)
    for ii=0,nut-1 do begin
        times[*,ii] = uts[ii]+dts
        b_vecs[*,ii,*] = transpose(b_raws[ii,*,*])
    endfor
    times = reform(times, ntime)
    b_vecs = reform(b_vecs,[ntime,ndim])
    settings = dictionary('coord', default_coord)
    var_info = var_store(var_info, b_vecs, times)
    add_setting, var_info, smart=1, settings, id='bfield'

    return, var_info

end

function tracers_read_bfield, input_time_range, probe=probe, $
    update=update, get_name=get_name, errmsg=errmsg, _extra=extra
    compile_opt idl2

    errmsg = ''
    retval = !null
    prefix = 'ts'+probe+'_'

    default_coord = 'ts_fac'
    var_info = prefix+'b_'+default_coord
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then del_data, var_info
    time_range = time_double(input_time_range)
    if ~check_if_update(var_info, time_range) then return, var_info

    ; Load the files.
    files = tracers_load_mag(time_range, probe=probe, errmsg=errmsg, id='l2%bdc-16sps')
    if errmsg ne '' then return, retval

    var_list = list()
    in_vars = prefix+'l2_mag_16sps_fac_deltab'
    out_vars = var_info
    time_var = prefix+'l2_mag_16sps_epoch'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'out_vars', out_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2000' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    settings = dictionary($
        'coord', default_coord, $
        'requested_time_range', time_range )
    add_setting, var_info, smart=1, settings, id='bfield'

    return, var_info

end

compile_opt idl2
time_range = ['2025-11-22','2025-11-23']
time_range = ['2025-12-22','2025-12-23']
probe = '2'
var = tracers_read_bfield(time_range, probe=probe)
tplot, var, trange=time_range
end
