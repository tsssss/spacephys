
function tracers_read_efield, input_time_range, probe=probe, $
    update=update, get_name=get_name
    compile_opt idl2

    errmsg = ''
    retval = !null
    prefix = 'ts'+probe+'_'

    default_coord = 'ts_efi'
    var_info = prefix+'e_'+default_coord
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then del_data, var_info
    time_range = time_double(input_time_range)
    if ~check_if_update(var_info, time_range) then return, var_info

    files = tracers_load_efi(time_range, probe=probe, errmsg=errmsg, id='l1b%vdc-roi_x274')
    if errmsg ne '' then return, retval

    var_list = list()
    in_vars = prefix+'l1b_vdc_'+['xplus_roi','xminus_roi','yplus_roi','yminus_roi']
    time_var = 'Epoch'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2001' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    vxp = transpose(var_get_data(in_vars[0], times=uts))
    vxm = transpose(var_get_data(in_vars[1]))
    vyp = transpose(var_get_data(in_vars[2]))
    vym = transpose(var_get_data(in_vars[3]))
    dts = cdf_read_var('EpochOffset', filename=files[0])*1e-9   ; ns to sec.
    ndt = n_elements(dts)
    nut = n_elements(uts)
    ntime = nut*ndt
    ndim = 3
    times = dblarr(ndt,nut)
    for ii=0,nut-1 do begin
        times[*,ii] = uts[ii]+dts
    endfor
    times = reform(times, ntime)
    boom_lengths = [3,3,0d] ; m
    e_xyz = fltarr(ntime,ndim)
    e_xyz[*,0] = reform(vxp-vxm, ntime)/boom_lengths[0]
    e_xyz[*,1] = reform(vyp-vym, ntime)/boom_lengths[1]
    e_xyz *= 1e3    ; from V/m to mV/m
    settings = dictionary('coord', default_coord)
    var_info = var_store(var_info, e_xyz, times)
    add_setting, var_info, smart=1, settings, id='efield'
    
    return, var_info

end

compile_opt idl2
time_range = ['2025-11-22','2025-11-23']
time_range = ['2025-12-22/15:24','2025-12-22/15:45']
time_range = ['2025-12-22/15:24','2025-12-22/15:30']
probe = '2'
e_var = tracers_read_efield(time_range, probe=probe)
b_var = tracers_read_bfield(time_range, probe=probe)
spec_var = tracers_read_ion_en_spec(time_range, probe=probe)
r_var = tracers_read_orbit(time_range, probe=probe)
mlat_vars = lets_read_mlat_vars(r_var)
mlat_var = mlat_vars['mlat']
mlt_var = mlat_vars['mlt']
plot_vars = [spec_var,e_var,b_var,mlat_var,mlt_var]
tplot, plot_vars, trange=time_range
end
