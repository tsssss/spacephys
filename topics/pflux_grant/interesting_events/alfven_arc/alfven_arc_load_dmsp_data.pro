
function alfven_arc_load_dmsp_data, input_time_range, probe=probe, filename=data_file, _extra=ex

    time_range = time_double(input_time_range)
    if n_elements(probe) eq 0 then message, 'No input probe ...'

    if n_elements(data_file) eq 0 then begin
        base = time_string(time_range[0],tformat='YYYY_MMDD_hh')+'_dmsp'+probe+'_data_v01.cdf'
        data_file = join_path([googledir(),'works','pflux_grant','alfven_arc','data',base])
    endif

    if file_test(data_file) eq 0 then begin
        cdf_save_setting, 'probe', probe, filename=data_file
        cdf_save_setting, 'prefix', 'dmsp'+probe+'_', filename=data_file
        cdf_save_setting, 'time_range', time_range, filename=data_file
    endif
    event_info = cdf_read_setting(filename=data_file)
    event_info['data_file'] = data_file

    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']


;---SSUSI.
    mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, errmsg=errmsg, id='lbhs')

stop
    return, data_file

end


time_range = ['2015-03-12/08:45','2015-03-12/09:10']
probe = 'f19'
print, alfven_arc_load_dmsp_data(time_range, probe=probe)
end