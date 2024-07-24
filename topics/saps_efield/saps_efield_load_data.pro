function saps_efield_load_data, id;, event_info=event_info

;    if n_elements(event_info) ne 0 then if event_info.id eq id then return, event_info

    topic_id = '2024_saps_efield'
    plot_dir = join_path([googledir(),'works',topic_id,'plot',id])
    data_dir = join_path([googledir(),'works',topic_id,'data'])
    version = 'v01'

    if id eq '2015_0104' then begin

    ;---Overall setting.
        event_time_range = time_double(['2015-01-04/10:00','2015-01-04/14:00'])
        rbsp_probes = ['b']
        dmsp_probes = ['f16']

        event_info = dictionary($
            'time_range', event_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'rbsp_probes', rbsp_probes, $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )

    ;---Ground settings.
        ground_file = join_path([data_dir,id+'_ground_data_'+version+'.cdf'])
        ground_time_range = time_double(['2015-01-04/10:00','2015-01-04/14:00'])
        event_info.ground_time_range = ground_time_range

        sites = ['fykn','whit','fsim','atha','gill']
        sites = ['fykn','whit','atha','gill']
        nsite = n_elements(sites)
        min_elevs = [10d,3,8,8]

        mlt_range = [0d,6]
        mlat_range = [55d,75]
        merge_method = 'merge_elev'
        calibration_method = 'moon'
        asi_setting = dictionary($
            'data_file', ground_file, $
            'time_range', ground_time_range, $
            'mlt_range', mlt_range, $
            'mlat_range', mlat_range, $
            'sites', sites, $
            'min_elevs', min_elevs, $
            'merge_method', merge_method, $
            'calibration_method', calibration_method )
        
        mlt_image_var = themis_asf_read_mlt_image($
            ground_time_range, sites=sites, $
            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
        event_info['asi_setting'] = asi_setting

        mlt_range = asi_setting['mlt_range']
        mlat_range = asi_setting['mlat_range']
        mlt_image_rect_var = themis_asf_read_mlt_image_rect($
            ground_time_range, sites=sites, $
            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method, $
            mlt_range=mlt_range, mlat_range=mlat_range)
            
;        mlt_image_var = lets_read_this($
;            func='themis_asf_read_mlt_image', $
;            save_to=ground_file, $
;            ground_time_range, sites=sites, $
;            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
;        event_info['asi_setting'] = asi_setting
;
;        mlt_range = asi_setting['mlt_range']
;        mlat_range = asi_setting['mlat_range']
;        mlt_image_rect_var = lets_read_this($
;            func='themis_asf_read_mlt_image_rect', $
;            save_to=ground_file, $
;            ground_time_range, sites=sites, $
;            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method, $
;            mlt_range=mlt_range, mlat_range=mlat_range)

    ;---RBSP.
        rbsp_time_range = ['2015-01-04/06:00','2015-01-04/16:00']
        rbsp_colors = sgcolor(['magenta','purple'])
        foreach probe, rbsp_probes, probe_id do begin
            rbsp_file = join_path([data_dir,id+'_rbsp'+probe+'_data_'+version+'.cdf'])
            
            ; Remove bad E field.
            bad_e_trs = list()
            if probe eq 'b' then begin
                bad_e_trs.add, ['2015-01-04/11:07:26','2015-01-04/11:12:30']
                bad_e_trs.add, ['2015-01-04/11:34:33','2015-01-04/11:39:24']
                bad_e_trs.add, ['2015-01-04/12:44:00','2015-01-04/13:12:00']
                bad_e_trs.add, ['2015-01-04/13:17:11','2015-01-04/13:19:54']
            endif
            
            sc_info = saps_efield_load_rbsp_data($
                rbsp_time_range, filename=rbsp_file, probe=probe, bad_e_trs=bad_e_trs)
            sc_info['sc_name'] = 'RBSP'
            sc_info['sc_color'] = rbsp_colors[probe_id]
            key = 'rbsp'+probe
            event_info.rbsp[key] = sc_info
            
            

            ; Spinfit E field.
            prefix = 'rbsp'+probe+'_'
            rbsp_efw_phasef_read_e_spinfit, time_double(rbsp_time_range), probe=probe
            e0_mgse_var = prefix+'e_spinfit_mgse_v12'
            options, e0_mgse_var, coord='rbsp_mgse', $
                coord_labels=constant('xyz'), spin_axis='x', probe=probe
            b_var = prefix+'b0_gsm'
            edot0_mgse_var = lets_calc_edotb0(e_var=e0_mgse_var, b_var=b_var, $
                var_info=prefix+'edot0_spinfit_rbsp_mgse', update=update)
            
            b_angle = get_var_setting(edot0_mgse_var, 'b_angle')
            times = get_var_time(edot0_mgse_var)
            b_angle_var = prefix+'edot0_b_angle'
            store_data, b_angle_var, times, b_angle
            b_angle_limit = 15d
            index = where(abs(b_angle) lt b_angle_limit, count)
            if count ne 0 then begin
                e_vec = get_var_data(edot0_mgse_var, times=times)
                e_vec[index,*] = !values.f_nan
                store_data, edot0_mgse_var, times, e_vec
            endif
            
            q_gsm2fac_var = prefix+'q_gsm2fac'
            var = edot0_mgse_var
            in_coord = strlowcase(get_setting(var,'coord'))
            out_var = streplace(var,in_coord,'fac')
            var_fac = lets_cotran([in_coord,'fac'], input=var, output=out_var, q_var=q_gsm2fac_var, update=update)

        endforeach
    
    ;---DMSP.
        dmsp_time_range = time_double(['2015-01-04/12:25','2015-01-04/13:55'])
        ssusi_time = time_double(['2015-01-04/12:42','2015-01-04/13:32'])
        foreach probe, dmsp_probes do begin
            dmsp_file = join_path([data_dir,id+'_dmsp'+probe+'_data_'+version+'.cdf'])
            sc_info = saps_efield_load_dmsp_data($
                dmsp_time_range, filename=dmsp_file, probe=probe)
            sc_info['sc_name'] = 'DMSP'
            sc_info['ssusi_time'] = ssusi_time
            key = 'dmsp'+probe
            event_info.dmsp[key] = sc_info
        endforeach
    
    endif else if id eq '2013_0317' then begin
        
    ;---Overall setting.
        event_time_range = time_double(['2013-03-17/08:00','2013-03-17/10:30'])
        rbsp_probes = ['a','b']
        dmsp_probes = ['f18']

        event_info = dictionary($
            'time_range', event_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'rbsp_probes', rbsp_probes, $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )
    
    ;---Ground settings.
        ground_file = join_path([data_dir,id+'_ground_data_'+version+'.cdf'])
        ground_time_range = time_double(['2013-03-17/05:30','2013-03-17/09:30'])
        event_info.ground_time_range = ground_time_range

;        weygand_setting = dictionary($
;            'data_file', ground_file, $
;            'time_range', ground_time_range )
;        foreach type, ['hor','ver'] do begin
;            j_vars = lets_read_this($
;                func='themis_read_weygand_j', $
;                save_to=ground_file, $
;                ground_time_range, id='j_'+type)
;        endforeach
;        event_info['weygand_setting'] = weygand_setting

        sites = ['kapu','snkq','gill','pina','fsmi','fsim']
        nsite = n_elements(sites)
        min_elevs = fltarr(nsite)+5
        index = where(sites eq 'rank', count)
        if count ne 0 then min_elevs[index] = 10
        index = where(sites eq 'fsmi', count)
        if count ne 0 then min_elevs[index] = 10

        mlt_range = [0d,6]
        mlat_range = [55d,75]
        merge_method = 'merge_elev'
        calibration_method = 'simple'
        asi_setting = dictionary($
            'data_file', ground_file, $
            'time_range', ground_time_range, $
            'mlt_range', mlt_range, $
            'mlat_range', mlat_range, $
            'sites', sites, $
            'min_elevs', min_elevs, $
            'merge_method', merge_method, $
            'calibration_method', calibration_method )
        
        mlt_image_var = lets_read_this($
            func='themis_asf_read_mlt_image', $
            save_to=ground_file, $
            ground_time_range, sites=sites, $
            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
        event_info['asi_setting'] = asi_setting

        mlt_range = asi_setting['mlt_range']
        mlat_range = asi_setting['mlat_range']
        mlt_image_rect_var = lets_read_this($
            func='themis_asf_read_mlt_image_rect', $
            save_to=ground_file, $
            ground_time_range, sites=sites, $
            min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method, $
            mlt_range=mlt_range, mlat_range=mlat_range)

        
    ;---RBSP.
        rbsp_time_range = ['2013-03-17','2013-03-17/12:00']
        rbsp_colors = sgcolor(['magenta','purple'])
        foreach probe, rbsp_probes, probe_id do begin
            rbsp_file = join_path([data_dir,id+'_rbsp'+probe+'_data_'+version+'.cdf'])
            
            sc_info = saps_efield_load_rbsp_data($
                rbsp_time_range, filename=rbsp_file, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
            sc_info['sc_name'] = 'RBSP'
            sc_info['sc_color'] = rbsp_colors[probe_id]
            key = 'rbsp'+probe
            event_info.rbsp[key] = sc_info
        endforeach
    
    ;---DMSP.
        dmsp_time_range = time_double(['2013-03-17/09:20','2013-03-17/09:45'])
        ssusi_time = time_double('2013-03-17/09:32')
        foreach probe, dmsp_probes do begin
            dmsp_file = join_path([data_dir,id+'_dmsp'+probe+'_data_'+version+'.cdf'])
            sc_info = saps_efield_load_dmsp_data($
                dmsp_time_range, filename=dmsp_file, probe=probe)
            sc_info['sc_name'] = 'DMSP'
            sc_info['ssusi_time'] = ssusi_time
            key = 'dmsp'+probe
            event_info.dmsp[key] = sc_info
        endforeach

    endif else if id eq '2013_0601' then begin
        
        ;---Overall setting.
        event_time_range = time_double(['2013-06-01/07:00','2013-06-01/10:30'])
        rbsp_probes = ['a','b']
        dmsp_probes = ['f18']

        event_info = dictionary($
            'time_range', event_time_range, $
            'id', id, $
            'plot_dir', plot_dir, $
            'data_dir', data_dir, $
            'version', version, $
            'rbsp_probes', rbsp_probes, $
            'rbsp', dictionary(), $
            'dmsp', dictionary(), $
            'ground', dictionary() )


    ;---RBSP.
        rbsp_time_range = ['2013-06-01','2013-06-02']
        rbsp_colors = sgcolor(['magenta','purple'])
        foreach probe, rbsp_probes, probe_id do begin
            rbsp_file = join_path([data_dir,id+'_rbsp'+probe+'_data_'+version+'.cdf'])

            sc_info = low_lshell_outflow_load_rbsp_data($
                rbsp_time_range, filename=rbsp_file, probe=probe, bad_e_time_ranges=bad_e_time_ranges)
            sc_info['sc_name'] = 'RBSP'
            sc_info['sc_color'] = rbsp_colors[probe_id]
            key = 'rbsp'+probe
            event_info.rbsp[key] = sc_info
        endforeach

        
    endif

    return, event_info


end


id = '2015_0317'
event_info = low_lshell_outflow_load_data(id)
end