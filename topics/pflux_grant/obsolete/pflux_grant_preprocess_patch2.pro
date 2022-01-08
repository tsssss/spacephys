;+
; Convert double data to float.
; Need to write data to a temporary file, cdf_del_var does not change file size.
;-

pro pflux_grant_preprocess_patch2, probe=probe, project=project

    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project.pflux_calc_setting

    if n_elements(probe) eq 0 then return
    if (where(probe eq ['a','b']))[0] eq -1 then return
    rbspx = 'rbsp'+probe
    mission_time_range = pflux_calc_setting[rbspx].time_range;>time_double('2013-10-29')
    secofday = constant('secofday')
    nday = total(mission_time_range*[-1,1])/secofday
    dates = mission_time_range[0]+findgen(nday)*secofday

;---Some settings.
    prefix = rbspx+'_'

    local_root = join_path([default_local_root(),'sdata','rbsp'])
    paths = [rbspx,'ebfield']
    vars = prefix+['b1_gsm','e_uv','b0_gsm','ew','bw_ratio']
    time_vars = ['ut_sec','ut_orbit']

    foreach date, dates do begin
        lprmsg, 'Processing '+time_string(date)+' ...'
        time_range = date+[0,secofday]
        base_name = rbspx+'_ebfield_'+time_string(date,tformat='YYYY_MMDD')+'_v01.cdf'
        year = time_string(date,tformat='YYYY')
        data_file = join_path([local_root,paths,year,base_name])
        if file_test(data_file) eq 0 then continue
        tmp_file = join_path([homedir(),'tmp.cdf'])
        if file_test(tmp_file) then file_delete, tmp_file
        in_id = cdf_open(data_file)
        out_id = cdf_create(tmp_file)
        foreach var, time_vars do begin
            if ~cdf_has_var(var, filename=in_id) then continue
            data = cdf_read_var(var, filename=in_id)
            setting = cdf_read_setting(var, filename=in_id)
            cdf_save_var, var, value=data, filename=out_id
            cdf_save_setting, setting, varname=var, filename=out_id
        endforeach
        foreach var, vars do begin
            if ~cdf_has_var(var, filename=in_id) then continue
            data = float(cdf_read_var(var, filename=in_id))
            setting = cdf_read_setting(var, filename=in_id)
            cdf_save_var, var, value=data, filename=out_id
            cdf_save_setting, setting, varname=var, filename=out_id
        endforeach
        cdf_close, in_id
        cdf_close, out_id
        file_copy, tmp_file, data_file, /overwrite
    endforeach
    if file_test(tmp_file) then file_delete, tmp_file
end

foreach probe, ['b'] do pflux_grant_preprocess_patch2, probe=probe, project=project
end
