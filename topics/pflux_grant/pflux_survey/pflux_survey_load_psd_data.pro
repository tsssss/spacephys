;+
; Calculate psd for dE, dB, and S.
;-

pro pflux_survey_load_psd_data, probe=probe, project=project

    project = pflux_survey_load_project()
    probe_infos = project['probe_infos']
    rbspx = 'rbsp'+probe
    prefix = rbspx+'_'

    data_file = join_path([project.data_dir,'pflux_survey_psd_data.cdf'])
    if file_test(data_file) eq 0 then begin
        path = fgetpath(data_file)
        if file_test(path) eq 0 then file_mkdir, path
        settings = dictionary($
            'name', 'power spectral density data per orbit')
        cdf_save_setting, settings, filename=data_file
    endif


;---Load orbit data.
    pflux_survey_load_data, 'dis', probe=probe


;---Orbit times.
    the_var = prefix+'orbit_time_range'
    if ~cdf_has_var(the_var, filename=data_file) then begin
        mission_time_range = probe_infos[rbspx].time_range
        data = pflux_grant_calculate_orbit_time_range(mission_time_range, probe=probe)
        settings = dictionary($
            'name', 'orbit time ranges b/w perigees', $
            'unit', 'sec')
        cdf_save_var, the_var, value=data, filename=data_file
        cdf_save_setting, settings, varname=the_var, filename=data_file
    endif
    orbit_time_ranges = cdf_read_var(the_var, filename=data_file)

    the_var = prefix+'orbit_time'
    if ~cdf_has_var(the_var, filename=data_file) then begin
        data = (orbit_time_ranges[*,0]+orbit_time_ranges[*,1])*0.5
        settings = dictionary($
            'name', 'time at orbit center', $
            'unit', 'sec')
        cdf_save_var, the_var, value=data, filename=data_file
        cdf_save_setting, settings, varname=the_var, filename=data_file
    endif
    orbit_times = cdf_read_var(the_var, filename=data_file)
    norbit = n_elements(orbit_times)


;---Scales, same for -A and -B.
    the_var = 'psd_scale'
    if ~cdf_has_var(the_var, filename=data_file) then begin
        settings = dictionary()
        period_range = [1d,1800]
        p2s = (wavelet_info('morlet'))[7]   ; t2s.
        scale_range = period_range*p2s
        settings['period_range'] = period_range
        settings['scale_range'] = scale_range
        dj = 1d/8
        scales = smkgmtrc(scale_range[0], scale_range[1], 2^dj, 'dx')
        nscale = n_elements(scales)
        periods = scales/p2s
        freqs = 1/periods
        settings['nscale'] = nscale
        settings['periods'] = periods
        settings['freqs'] = freqs
        settings['dj'] = dj

        data = scales
        cdf_save_var, the_var, value=data, filename=data_file
        cdf_save_setting, settings, varname=the_var, filename=data_file
    endif
    scales = cdf_read_var(the_var, filename=data_file)
    nscale = n_elements(scales)


;---Load PSD data.
    var_types = ['pf_fac_norm']
    min_dis = 4.

    foreach var_type, var_types do begin
        the_var = prefix+'pf_fac_psd'
        if ~cdf_has_var(the_var, filename=data_file) then begin
            psd_data = fltarr(norbit,nscale)
            psd_time_ranges = fltarr(norbit,2)
            for orbit_id=0, norbit-1 do begin
                orbit_time_range = orbit_time_ranges[orbit_id,*]
                lprmsg, 'Processing '+strjoin(time_string(orbit_time_range),' to ')+' ...'
                lprmsg, systime()

                dis = get_var_data(prefix+'dis', in=orbit_time_range, times=orbit_times)
                time_range = minmax(orbit_times[where(dis ge min_dis)])

                case var_type of
                    'pf_fac_norm': pflux_grant_read_preprocessed_pflux, time_range, probe=probe
                    else: ; do nothing.
                endcase

                tmp_var = prefix+'tmp_data'
                tmp_data = snorm(get_var_data(prefix+'pf_fac_norm', times=times))
                store_data, tmp_var, times, tmp_data
                calc_psd, tmp_var, scales=scales
                get_data, tmp_var+'_psd_cwt', freqs, psds
                ;plot, scales*2, psds, xlog=1, ylog=1
                psd_data[orbit_id,*] = psds
                psd_time_ranges[orbit_id,*] = time_range
            endfor

            time_var = prefix+'psd_time_range'
            settings = dictionary($
                'name', 'time range per orbit PSD is calculaed', $
                'depend_0', prefix+'orbit_time', $
                'unit', 'sec')
            cdf_save_var, time_var, value=psd_time_ranges, filename=data_file
            cdf_save_setting, settings, varname=time_var, filename=data_file

            settings = dictionary($
                'depend_0', prefix+'orbit_time', $
                'depend_1', 'psd_scale')
            cdf_save_var, the_var, value=psd_data, filename=data_file
            cdf_save_setting, settings, varname=the_var, filename=data_file
        endif
        cdf_load_var, the_var, filename=data_file
        scales = cdf_read_var('psd_scale', filename=data_file)
        settings = cdf_read_setting('psd_scale', filename=data_file)
        settings['ylog'] = 1
        settings['ytitle'] = '(Hz)'
        settings['display_type'] = 'spec'
        settings['zlog'] = 1
        settings['ztitle'] = '(mW/m!U2!N)/Hz'
        get_data, the_var, times, data
        store_data, the_var, times, data, settings.freqs
        add_setting, the_var, /smart, settings
    endforeach



end

probes = ['b']
foreach probe, probes do pflux_survey_load_psd_data, probe=probe, project=project
end