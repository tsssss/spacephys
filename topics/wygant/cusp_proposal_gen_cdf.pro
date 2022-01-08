;+
; Generate CDF containing ion and ele energy flux.
;-

data_dir = join_path([homedir(),'cusp_energy_flux'])
dates = time_double([$
    '2003-10-29', $
    '2001-11-06', $
    '2000-08-12', $
    '2004-11-08', $
    '2004-11-07', $
    '2003-08-18', $
    '2000-10-29', $
    '2000-07-16', $
    '2005-08-24', $
    '1998-08-26', $
    '2003-08-31', $
    '2001-04-12', $
    '2001-10-22', $
    '2001-11-06', $
    '2003-10-31', $
    '2003-11-04', $
    '2002-05-23', $
    '2002-03-24', $
    '2000-04-06'])
dates = sort_uniq(dates)


data_dir = join_path([homedir(),'cusp_energy_flux','%Y'])
full_time_range = time_double(['1997-01-01','2004-01-01'])
dates = make_bins(full_time_range,86400d)


data_file_pattern = 'po_hyd_energy_flux_%Y%m%d_v01.cdf'

foreach date, dates do begin
    lprmsg, 'Processing '+time_string(date)+'...'
    store_data, '*', /delete
    time_range = date+[0,constant('secofday')]
    polar_read_eflux, time_range, errmsg=errmsg

    if errmsg eq 'No data file' then continue

    data_file = join_path([data_dir,data_file_pattern])
    data_file = apply_time_to_pattern(data_file, date)
    ;if file_test(data_file) eq 1 then continue
    if file_test(data_file) eq 1 then file_delete, data_file
    lprmsg, 'Saving data to '+data_file+'...'

    get_data, 'po_ion_energy_flux', times
    utname = 'ut_hydra'
    settings = dictionary($
        'FIELDNAM', 'unix timestamp', $
        'VAR_TYPE', 'support_data', $
        'UNITS', 'sec')
    cdf_save_var, utname, value=times, filename=data_file
    cdf_save_setting, settings, var=utname, filename=data_file

    varname = 'po_ion_energy_flux'
    data = get_var_data('po_ion_energy_flux')
    settings = dictionary($
        'FIELDNAM', 'Ion energy flux', $
        'VAR_TYPE', 'data', $
        'UNITS', 'mW/m!U2!N', $
        'DEPEND_0', utname, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=data_file
    cdf_save_setting, settings, var=varname, filename=data_file

    varname = 'po_ele_energy_flux'
    data = get_var_data('po_ele_energy_flux')
    settings = dictionary($
        'FIELDNAM', 'Ele energy flux', $
        'VAR_TYPE', 'data', $
        'UNITS', 'mW/m!U2!N', $
        'DEPEND_0', utname, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=data_file
    cdf_save_setting, settings, var=varname, filename=data_file

;---Oribt data.
    polar_read_orbit, time_range

    utname = 'ut_orbit'
    get_data, 'po_r_gsm', times
    settings = dictionary($
        'FIELDNAM', 'unix timestamp', $
        'VAR_TYPE', 'support_data', $
        'UNITS', 'sec')
    cdf_save_var, utname, value=times, filename=data_file
    cdf_save_setting, settings, var=utname, filename=data_file

    varname = 'po_r_gsm'
    data = get_var_data(varname)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft position in GSM', $
        'VAR_TYPE', 'data', $
        'UNITS', 'Re', $
        'DEPEND_0', utname, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=data_file
    cdf_save_setting, settings, var=varname, filename=data_file

    varname = 'po_mlt'
    data = get_var_data(varname)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft MLT', $
        'VAR_TYPE', 'data', $
        'UNITS', 'hr', $
        'DEPEND_0', utname, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=data_file
    cdf_save_setting, settings, var=varname, filename=data_file

    varname = 'po_mlat'
    data = get_var_data(varname)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft MLat', $
        'VAR_TYPE', 'data', $
        'UNITS', 'deg', $
        'DEPEND_0', utname, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=data_file
    cdf_save_setting, settings, var=varname, filename=data_file

    varname = 'po_ilat'
    data = get_var_data(varname)
    settings = dictionary($
        'FIELDNAM', 'Spacecraft Invariant Lat', $
        'VAR_TYPE', 'data', $
        'UNITS', 'deg', $
        'DEPEND_0', utname, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=data_file
    cdf_save_setting, settings, var=varname, filename=data_file
endforeach

end
