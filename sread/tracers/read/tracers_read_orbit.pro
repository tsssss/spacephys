
function tracers_read_orbit, input_time_range, probe=probe, $
    update=update, get_name=get_name
    compile_opt idl2

    errmsg = ''
    retval = !null
    prefix = 'ts'+probe+'_'

    default_coord = 'geo'   ; To be determined.
    var_info = prefix+'r_'+default_coord
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then del_data, var_info
    time_range = time_double(input_time_range)
    if ~check_if_update(var_info, time_range) then return, var_info

    files = tracers_load_aci(time_range, probe=probe, errmsg=errmsg, id='l1b%ipd_x292')
    if errmsg ne '' then return, retval

    var_list = list()
    in_vars = 'location'
    time_var = 'Epoch'
    var_list.add, dictionary($
        'in_vars', in_vars, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2001' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    locs = var_get_data(in_vars, times=times)
    re = constant('re')
    rad = constant('rad')
    diss = locs[*,0]/re
    colats = locs[*,1]*rad
    lons = locs[*,2]*rad
    ntime = n_elements(times)
    ndim = 3
    r_geos = fltarr(ntime,ndim)
    r_geos[*,2] = diss*cos(colats)
    r_xys = diss*sin(colats)
    r_geos[*,0] = r_xys*cos(lons)
    r_geos[*,1] = r_xys*sin(lons)
    settings = dictionary('coord',default_coord)
    var_info = var_store(var_info, r_geos, times, id='position', settings=settings)

    return, var_info

end

compile_opt idl2
time_range = ['2025-12-22','2025-12-23']
time_range = ['2025-12-22/15:24','2025-12-22/15:45']
probe = '2'
spec_var = tracers_read_ion_en_spec(time_range, probe=probe)
r_var = tracers_read_orbit(time_range, probe=probe)
mlat_vars = lets_read_mlat_vars(r_var)
mlat_var = mlat_vars['mlat']
mlt_var = mlat_vars['mlt']
plot_vars = [spec_var,mlat_var,mlt_var]
tplot, plot_vars, trange=time_range
end