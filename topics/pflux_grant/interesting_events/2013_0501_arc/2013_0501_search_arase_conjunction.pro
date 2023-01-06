;+
; Search for RBSP-A or -B conjunction with Arase.
;-

probes = ['a','b']
search_time_range = time_double(['2017-01-01','2019-12'])


log_file = join_path([srootdir(),'2015_0501_search_arase_conjunction.txt'])
if file_test(log_file) eq 0 then begin
    time_step = 60
    common_times = make_bins(search_time_range, time_step)
    rad = constant('rad')
    deg = constant('deg')
    max_dmlt = 1.
    max_dangle = max_dmlt*15
    tab = '    '
    
    ftouch, log_file

    arase_r_var = arase_read_orbit(search_time_range, get_name=1, coord='sm')
    if check_if_update(arase_r_var, search_time_range) then begin
        arase_r_var = arase_read_orbit(search_time_range, coord='sm')
        interp_time, arase_r_var, common_times
    endif
    r_sm = get_var_data(arase_r_var)
    
    arase_mlt = pseudo_mlt(r_sm)*15*rad
    arase_dis = snorm(r_sm)
    arase_cosmlt = cos(arase_mlt)
    arase_sinmlt = sin(arase_mlt)
    
    arase_mlt_var = 'arase_mlt'
    store_data, arase_mlt_var, common_times, arase_mlt
    add_setting, arase_mlt_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'MLT', $
        'unit', 'h' )
    
    arase_dis_var = 'arase_dis'
    store_data, arase_dis_var, common_times, arase_dis
    add_setting, arase_dis_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', '|R|', $
        'unit', 'Re' )


    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        rbsp_r_var = rbsp_read_orbit(search_time_range, get_name=1, coord='sm', probe=probe)
        if check_if_update(rbsp_r_var, search_time_range) then begin
            rbsp_r_var = rbsp_read_orbit(search_time_range, probe=probe, coord='sm')
            interp_time, rbsp_r_var, common_times
        endif
        r_sm = get_var_data(rbsp_r_var)
        
        rbsp_mlt = pseudo_mlt(r_sm)*15*rad
        rbsp_dis = snorm(r_sm)
        rbsp_cosmlt = cos(rbsp_mlt)
        rbsp_sinmlt = sin(rbsp_mlt)
        
        rbsp_mlt_var = prefix+'mlt'
        store_data, rbsp_mlt_var, common_times, rbsp_mlt
        add_setting, rbsp_mlt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'MLT', $
            'unit', 'h' )            
        
        rbsp_dis_var = prefix+'dis'
        store_data, rbsp_dis_var, common_times, rbsp_dis
        add_setting, rbsp_dis_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', '|R|', $
            'unit', 'Re' )

        dangle = acos(rbsp_cosmlt*arase_cosmlt+rbsp_sinmlt*arase_sinmlt)*deg
        index = where(dangle le max_dangle, count)
        if count eq 0 then continue
        time_ranges = common_times[time_to_range(index, time_step=1)]
        durations = (time_ranges[*,1]-time_ranges[*,0])/60
        nrec = n_elements(durations)

        for ii=0,nrec-1 do begin
            index = reform((time_ranges[ii,*]-common_times[0])/time_step)
            rbsp_mlt_range = rbsp_mlt[index]
            arase_mlt_range = arase_mlt[index]
            rbsp_dis_range = rbsp_dis[index]
            arase_dis_range = arase_dis[index]
            mean_dmlt = mean(dangle[index[0]:index[1]],nan=1)/15
            msg = strupcase(probe)+tab+$
                strjoin(reform(time_string(time_ranges[ii,*])),',')+tab+$
                string(durations[ii],format='(F6.2)')+tab+$
                string(mean_dmlt,format='(F6.2)')+tab+$
                strjoin(string(rbsp_mlt_range,format='(F6.2)'),',')+tab+$
                strjoin(string(arase_mlt_range,format='(F6.2)'),',')+tab+$
                strjoin(string(rbsp_dis_range,format='(F4.2)'),',')+tab+$
                strjoin(string(arase_dis_range,format='(F4.2)'),',')
            lprmsg, msg, log_file
        endfor
    endforeach
endif

lines = read_all_lines(log_file)
nline = n_elements(lines)
mlt_range = [-1,1]*5
min_duration = 10*60
min_duration = 0
min_dmlt = 0.5
min_dis = 3.5

flags = fltarr(nline)
for ii=0,nline-1 do begin
    line = lines[ii]
    infos = strsplit(line,', ',extract=1)
    duration = float(infos[3])*60
    if duration le min_duration then flags[ii] = 1
    rbsp_mlt_range = float(infos[5:6])
    arase_mlt_range = float(infos[7:8])
    mlt_ranges = [rbsp_mlt_range,arase_mlt_range]
    if min(mlt_ranges) le mlt_range[0] then flags[ii] = 1
    if max(mlt_ranges) ge mlt_range[1] then flags[ii] = 1
    dmlt = float(infos[4])
    if dmlt gt min_dmlt then flags[ii] = 1
    
    rbsp_dis_range = float(infos[9:10])
    arase_dis_range = float(infos[11:12])
    if min(rbsp_dis_range) le min_dis then flags[ii] = 1
    if min(arase_dis_range) le min_dis then flags[ii] = 1
endfor

index = where(flags eq 0, nline)
if nline eq 0 then stop
lines = lines[index]
for ii=0,nline-1 do begin
    infos = strsplit(lines[ii],', ',extract=1)
    time_range = time_double(infos[1:2])
    data_time_range = time_range+[-1,1]*3600
    probe = strlowcase(infos[0])
    prefix = 'rbsp'+probe+'_'
    
    arase_e_var = arase_read_efield(data_time_range, no_edotb=1, coord='dsi')
    if probe eq 'a' then begin
        rbsp_efw_phasef_read_e_spinfit_diagonal, data_time_range, probe=probe
        e_mgse_var = prefix+'e_spinfit_mgse_v24'
    endif else begin
        rbsp_efw_phasef_read_e_spinfit, data_time_range, probe=probe
        e_mgse_var = prefix+'e_spinfit_mgse_v12'
    endelse
;    e_uvw_var = prefix+'e_uvw'
;    e_uvw = get_var_data(e_uvw_var, times=times)
;    e_mgse = cotran(e_uvw, times, 'uvw2mgse', probe=probe)
;    e_mgse_var = prefix+'e_mgse'
;    store_data, e_mgse_var, times, e_mgse
    add_setting, e_mgse_var, smart=1, dictionary($
        'display_type', 'vector', $
        'unit', 'mV/m', $
        'short_name', 'E', $
        'coord', 'mGSE', $
        'yrange', [-1,1]*20, $
        'coord_labels', constant('xyz') )
    
    arase_r_var = arase_read_orbit(data_time_range, coord='sm')
    rbsp_r_var = rbsp_read_orbit(data_time_range, probe=probe, coord='sm')
    
    sgopen, 0, xsize=6, ysize=6
    tplot, [arase_e_var,e_mgse_var,arase_r_var,rbsp_r_var], trange=data_time_range
    timebar, time_range, linestyle=1, color=sgcolor('red')
    stop
endfor

end