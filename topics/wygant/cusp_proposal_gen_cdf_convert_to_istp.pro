;+
; Generate CDF containing ion and ele energy flux.
;-

data_dir = join_path([homedir(),'cusp_energy_flux'])
event_dates = time_double([$
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

data_dir = join_path([homedir(),'cusp_energy_flux','%Y'])
full_time_range = time_double(['1997-01-01','2004-01-01'])
dates = make_bins(full_time_range,86400d)
dates = [dates,event_dates]
dates = sort_uniq(dates)
data_file_pattern = 'po_hyd_energy_flux_%Y%m%d_v01.cdf'

foreach date, dates do begin
    lprmsg, 'Processing '+time_string(date)+'...'
    store_data, '*', /delete
    time_range = date+[0,constant('secofday')]
    polar_read_eflux, time_range, errmsg=errmsg

    if errmsg eq 'No data file' or errmsg eq 'No valid data' then continue

    data_file = join_path([data_dir,data_file_pattern])
    data_file = apply_time_to_pattern(data_file, date)
    ;if file_test(data_file) eq 1 then continue
    if file_test(data_file) eq 1 then file_delete, data_file
    lprmsg, 'Saving data to '+data_file+'...'
    path = file_dirname(data_file)
    if file_test(path) eq 0 then file_mkdir, path
    cdf_id = cdf_create(data_file, /col_major)

    global_settings = dictionary($
        'Project', 'ISTP>International Solar-Terrestrial Physics', $
        'Source_name', 'Polar>Polar Plasma Laboratory', $
        'Discipline', 'Space Physics>Magnetospheric Science', $
        'Data_type', 'moments-14sec>High Resolution data', $
        'Descriptor', 'HYDRA>Fast Plasma Analyzer', $
        'Data_version', '1', $
        'Logical_file_id', 'po_hyd_energy_flux_yyyymmdd_v01', $
        'PI_name', 'J. Scudder', $
        'PI_affiliation', 'U of Iowa', $
        'TEXT', 'Reference: HYDRA is a 3-Dimensional Electron and Ion Hot Plasma Instrument for the Polar Spacecraft of the GGS Mission, J. Scudder et al., Space Sci. Rev., 71,459-495, Feb. 1995. http://www-st.physics.uiowa.edu This data set contains survey electron and proton moments for the energy flux (parallel), at 13.8-second resolution as determined (0-20keV). Higher quality data products may be available from the P.I.', $
        'Instrument_type', 'Plasma and Solar Wind', $
        'Mission_group', 'Polar', $
        'Logical_source', 'po_hydra_energy_flux', $
        'Logical_source_description', 'Polar Fast Plasma Analyzer 13.8 second Resolution Moments', $
        'Time_resolution', '13.8 seconds' )
    cdf_save_setting, global_settings, filename=cdf_id

    get_data, 'po_ion_energy_flux', times
    utname = 'Epoch'
    epochs = stoepoch(times, 'unix')
    settings = dictionary($
        'FIELDNAM', 'Time since 0 AD', $
        'VAR_TYPE', 'support_data', $
        'CATDESC', 'Interval centered time', $
        'LABLAXIS', 'Epoch', $
        'UNITS', 'ms' )
        ;'VALIDMIN', 62987673600000.000, $
        ;'VALIDMAX', 63429523200000.000, $
        ;'FILLVAL', -1.00000e+31 )
    cdf_save_var, utname, value=epochs, filename=cdf_id, cdf_type='CDF_EPOCH'
    cdf_save_setting, settings, var=utname, filename=cdf_id

    fillval = float(-1e31)
    data_format = 'E11.3'
    varname = 'po_ion_energy_flux'
    data = float(get_var_data('po_ion_energy_flux'))
    settings = dictionary($
        'FIELDNAM', 'Ion energy flux', $
        'VAR_TYPE', 'data', $
        'CATDESC', 'Ion energy flux', $
        'LABLAXIS', 'Eflux ion', $
        'UNITS', 'mW/m!U2!N', $
        'DEPEND_0', utname, $
        'VALIDMIN', -10000.0, $
        'VALIDMAX', 10000.0, $
        'FILLVAL', fillval, $
        'FORMAT', data_format, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=cdf_id
    cdf_save_setting, settings, var=varname, filename=cdf_id

    varname = 'po_ele_energy_flux'
    data = float(get_var_data('po_ele_energy_flux'))
    settings = dictionary($
        'FIELDNAM', 'Ele energy flux', $
        'VAR_TYPE', 'data', $
        'CATDESC', 'Ele energy flux', $
        'LABLAXIS', 'Eflux ele', $
        'UNITS', 'mW/m!U2!N', $
        'DEPEND_0', utname, $
        'VALIDMIN', -10000.0, $
        'VALIDMAX', 10000.0, $
        'FILLVAL', fillval, $
        'FORMAT', data_format, $
        'DISPLAY_TYPE', 'time_series' )
    cdf_save_var, varname, value=data, filename=cdf_id
    cdf_save_setting, settings, var=varname, filename=cdf_id

;;---Oribt data.
;    polar_read_orbit, time_range
;
;    utname = 'epoch_orbit'
;    get_data, 'po_r_gsm', times
;    epochs = stoepoch(times, 'unix')
;    settings = dictionary($
;        'FIELDNAM', 'Time since 0 AD', $
;        'VAR_TYPE', 'support_data', $
;        'CATDESC', 'Interval centered time', $
;        'LABLAXIS', 'Epoch', $
;        'UNITS', 'ms', $
;        'VALIDMIN', 6.3367056e+13, $
;        'VALIDMAX', 6.3367142e+13, $
;        'FILLVAL', -1.00000e+31 )
;    cdf_save_var, utname, value=times, filename=cdf_id, cdf_type='CDF_EPOCH'
;    cdf_save_setting, settings, var=utname, filename=cdf_id
;
;    varname = 'po_r_gsm'
;    data = get_var_data(varname)
;    settings = dictionary($
;        'FIELDNAM', 'Spacecraft position in GSM', $
;        'VAR_TYPE', 'data', $
;        'UNITS', 'Re', $
;        'DEPEND_0', utname, $
;        'DISPLAY_TYPE', 'time_series' )
;    cdf_save_var, varname, value=data, filename=cdf_id
;    cdf_save_setting, settings, var=varname, filename=cdf_id
;
;    varname = 'po_mlt'
;    data = get_var_data(varname)
;    settings = dictionary($
;        'FIELDNAM', 'Spacecraft MLT', $
;        'VAR_TYPE', 'data', $
;        'UNITS', 'hr', $
;        'DEPEND_0', utname, $
;        'DISPLAY_TYPE', 'time_series' )
;    cdf_save_var, varname, value=data, filename=cdf_id
;    cdf_save_setting, settings, var=varname, filename=cdf_id
;
;    varname = 'po_mlat'
;    data = get_var_data(varname)
;    settings = dictionary($
;        'FIELDNAM', 'Spacecraft MLat', $
;        'VAR_TYPE', 'data', $
;        'UNITS', 'deg', $
;        'DEPEND_0', utname, $
;        'DISPLAY_TYPE', 'time_series' )
;    cdf_save_var, varname, value=data, filename=cdf_id
;    cdf_save_setting, settings, var=varname, filename=cdf_id
;
;    varname = 'po_ilat'
;    data = get_var_data(varname)
;    settings = dictionary($
;        'FIELDNAM', 'Spacecraft Invariant Lat', $
;        'VAR_TYPE', 'data', $
;        'UNITS', 'deg', $
;        'DEPEND_0', utname, $
;        'DISPLAY_TYPE', 'time_series' )
;    cdf_save_var, varname, value=data, filename=cdf_id
;    cdf_save_setting, settings, var=varname, filename=cdf_id

    cdf_close, cdf_id
endforeach

end
