;+
; Generate L2 burst1 data (E and B fields) v01 cdfs.
;-

pro rbsp_efw_phasef_gen_l2_burst1_v01, time_range, probe=probe, filename=file, log_file=log_file

    on_error, 0
    errmsg = ''

    msg = 'Processing '+file+' ...'
    lprmsg, msg, log_file

;---Check input.
    if n_elements(file) eq 0 then begin
        errmsg = 'cdf file is not set ...'
        lprmsg, errmsg, log_file
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        lprmsg, errmsg, log_file
        return
    endif
    if probe ne 'a' and probe ne 'b' then begin
        errmsg = 'Invalid probe: '+probe+' ...'
        lprmsg, errmsg, log_file
        return
    endif
    prefix = 'rbsp'+probe+'_'
    rbspx = 'rbsp'+probe


    vars = prefix+'efw_'+['eb1','mscb1']+'_mgse'
    del_data, vars
    rbsp_efw_read_l1_burst_efield, time_range, probe=probe, keep_spin_axis=1
    rbsp_efw_read_l1_burst_bfield, time_range, probe=probe
    nodata = 0
    foreach var, vars do begin
        if check_if_update(var) then nodata = 1
    endforeach
    if nodata then return

    ; Save data to file.
    vars = prefix+'efw_'+['eb1','mscb1']+'_mgse'
    version = 'v01'

    setting = dictionary($
        'HTTP_LINK', 'http://rbsp.space.umn.edu http://athena.jhuapl.edu', $
        'Logical_source', prefix+'efw-l2_burst1', $
        'LINK_TITLE', 'Daily Summary Plots and additional data Science Gateway (including Mission Rules of the Road for Data Usage)', $
        'PI_name', 'J. R. Wygant', $
        'Data_version', version, $
        'Mission_group', 'Van Allen Probes (RBSP)', $
        'Rules_of_use', 'Data for Scientific Use; please see EFW "Rules of the Road" (http://rbsp.space.umn.edu/data_policy.html) for details of the Van Allen Probes Mission and EFW Instrument "Rules of the Road", and the EFW-FAQ (http://rbsp.space.umn.edu/efw_faq.html) for details of and caveats on instrument performance and data quality.', $
        'File_naming_convention', 'source_datatype_descriptor_yyyyMMddthhmmss', $
        'Data_type', 'DC Efield and Bfield in MGSE Coord', $
        'MODS', '', $
        'Instrument_type', 'Electric Fields (space)', $
        'Source_name', 'RBSP-'+strupcase(probe)+'>Radiation Belt Storm Probe '+strupcase(probe), $
        'PI_affiliation', 'University of Minnesota', $
        'Time_resolution', 'UTC', $
        'Acknowledgement', "This work was supported by Van Allen Probes (RBSP) EFW funding provided by JHU/APL Contract No. 922613 under NASA's Prime Contract No. NNN06AA01C; EFW PI, J. R. Wygant, UMN.", $
        'Discipline', 'Space Physics>Magnetospheric Science', $
        'Generated_by', 'RBSP-EFW SOC', $
        'LINK_TEXT',  'EFW home page at Minnesota with Van Allen Probes', $
        'Generation_date', time_string(systime(1)), $
        'Logical_file_id', prefix+'_efw-l2_burst1_'+time_string(time_range[0],tformat='YYYYMMDDthhmmss')+'_'+version, $
        'Descriptor', 'EFW>Electric Field and Waves', $
        'TEXT', 'Contacts:  Tami.J.Kovalick@nasa.gov, Rita.C.Johnson@nasa.gov. Burst electric and magneti fields in the M-GSE coordinate system - see the EFW-FAQ (http://rbsp.space.umn.edu/efw_faq.html) for a description of the and the M-GSE coordinate system. The X-component of the E-field estimate, corresponding to the axial component in the spacecraft coordinate system, may contain DC offsets that are not physical. The nominal dynamic range of the E-field estimate is +/- 1 V/m in any component.', $
        'Project', 'RBSP>Radiation Belt Storm Probes' )
    cdf_save_setting, setting, filename=file


    e_var = prefix+'efw_eb1_mgse'
    get_data, e_var, times, e_mgse
    time_var = 'epoch_e'
    epochs = convert_time(times, from='unix', to='epoch16')
    cdf_save_var, time_var, value=epochs, filename=file, cdf_type='CDF_LONG_EPOCH'
    vatt = dictionary($
        'TIME_BASE', '0 AD', $
        'VAR_NOTES', 'Epoch tagged at the center of each interval, resolution varies from 512 to 16,384 Sampels/sec', $
        'MONOTON', 'INCREASE', $
        'FIELDNAM', 'Time, UTC', $
        'VALIDMIN', dcomplex(6.3397987e+10,0.0000000), $
        'DICT_KEY', 'time>Epoch', $
        'VAR_TYPE', 'support_data', $
        'SCALETYP', 'linear', $
        'FILLVAL', dcomplex(-1.0000000e+31,-1.0000000e+31), $
        'UNITS', 'ps (pico-second)', $
        'CATDESC', 'Time, UTC', $
        'VALIDMAX', dcomplex(6.6301114e+10,0.0000000) )
    cdf_save_setting, vatt, varname=time_var, filename=file

    cdf_save_var, e_var, value=float(e_mgse), filename=file
    setting = dictionary($
        'VAR_NOTES', 'Electric field in the MGSE coordinate system, DC offset removed', $
        'FIELDNAM', 'Efield(MGSE)', $
        'SCALETYP', 'linear', $
        'DEPEND_0', time_var, $
        'LABLAXIS', 'Efield(MGSE)', $
        'UNITS', 'mV/m', $
        'DISPLAY_TYPE', 'time_series', $
        'VALIDMAX', fltarr(3)+1e30, $
        'FORMAT', 'E13.6', $
        'LABL_PTR_1', 'efield_mgse_LABL', $
        'VALIDMIN', fltarr(3)-1e30, $
        'UNIT_PTR', 'efield_mgse_UNIT', $
        'VAR_TYPE', 'data', $
        'FILLVAL', -1.00000e+31, $
        'CATDESC', 'efield_in_corotation_frame_mgse' )
    cdf_save_setting, setting, varname=e_var, filename=file

    b_var = prefix+'efw_mscb1_mgse'
    b_mgse = get_var_data(b_var, at=times)
    cdf_save_var, b_var, value=float(b_mgse), filename=file

    setting = dictionary($
        'VAR_NOTES', 'Magnetic field in the MGSE coordinate system, DC offset removed', $
        'FIELDNAM', 'Bfield(nT)', $
        'SCALETYP', 'linear', $
        'DEPEND_0', time_var, $
        'LABLAXIS', 'Bfield(nT)', $
        'UNITS', 'nT', $
        'DISPLAY_TYPE', 'time_series', $
        'VALIDMAX', fltarr(3)+1e30, $
        'FORMAT', 'E13.6', $
        'LABL_PTR_1', 'bfield_mgse_LABL', $
        'VALIDMIN', fltarr(3)-1e30, $
        'UNIT_PTR', 'bfield_mgse_UNIT', $
        'VAR_TYPE', 'data', $
        'FILLVAL', -1.00000e+31, $
        'CATDESC', 'bfield_in_corotation_frame_mgse' )
    cdf_save_setting, setting, varname=b_var, filename=file


    xyz = constant('xyz')
    var = 'efield_mgse_LABL'
    cdf_save_var, var, value='E'+xyz+' MGSE', filename=file, save_as_one=1
    setting = dictionary($
        'FORMAT', 'A27', $
        'LABELAXIS', 'efield_mgse_LABL', $
        'FIELDNAM', 'efield_mgse_LABL', $
        'CATDESC', 'efield_mgse_LABL', $
        'VAR_TYPE', 'metadata' )
    cdf_save_setting, setting, varname=var, filename=file

    var = 'bfield_mgse_LABL'
    cdf_save_var, var, value='B'+xyz+' MGSE', filename=file, save_as_one=1
    setting = dictionary($
        'FORMAT', 'A27', $
        'LABELAXIS', 'bfield_mgse_LABL', $
        'FIELDNAM', 'bfield_mgse_LABL', $
        'CATDESC', 'bfield_mgse_LABL', $
        'VAR_TYPE', 'metadata' )
    cdf_save_setting, setting, varname=var, filename=file

    var = 'efield_mgse_UNIT'
    cdf_save_var, var, value=strarr(3)+'mV/m', filename=file, save_as_one=1
    setting = dictionary($
        'FORMAT', 'A4', $
        'LABELAXIS', 'efield_mgse_UNIT', $
        'FIELDNAM', 'efield_mgse_UNIT', $
        'CATDESC', 'efield_mgse_UNIT', $
        'VAR_TYPE', 'metadata' )
    cdf_save_setting, setting, varname=var, filename=file

    var = 'bfield_mgse_UNIT'
    cdf_save_var, var, value=strarr(3)+'nT', filename=file, save_as_one=1
    setting = dictionary($
        'FORMAT', 'A2', $
        'LABELAXIS', 'bfield_mgse_UNIT', $
        'FIELDNAM', 'bfield_mgse_UNIT', $
        'CATDESC', 'bfield_mgse_UNIT', $
        'VAR_TYPE', 'metadata' )
    cdf_save_setting, setting, varname=var, filename=file

end


date = '2013-06-07'
probe = 'a'

;secofday = constant('secofday')
;the_date = time_double(date)
;time_step = 15*60d
;the_times = make_bins(the_date+[0,secofday], time_step)
;ntime = n_elements(the_times)-1

time_range = time_double(['2013-06-07/02:15','2013-06-07/02:30'])
file = join_path([homedir(),'test_l2_burst1_v01.cdf'])
if file_test(file) eq 1 then file_delete, file
rbsp_efw_phasef_gen_l2_burst1_v01, time_range, probe=probe, filename=file
end
