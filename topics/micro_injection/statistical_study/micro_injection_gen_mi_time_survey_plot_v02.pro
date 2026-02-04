;+
; The raw mi_times. Need to filter with duration for example.
; Part is from mms_gen_survey_plot_v01
; v02 add ion velocity.
;-


pro micro_injection_gen_mi_time_survey_plot_v02, mission_probe, test=test


    project_info = micro_injection_stat_load_project()
    my_name = get_filename()
    version = get_file_version(my_name)
    plot_root = join_path([diskdir('data'),'survey_plot'])
    plot_root = join_path([plot_root,'microinjection','micro_injection_'+version])

    sample_energys = project_info['sample_energys']
    nsample_energy = n_elements(sample_energys)
    scale_info = project_info['scale_info']
    significant_flux_levels = [1.5,0.3]

    mission = mission_probe[0]
    probe = mission_probe[1]
    mission_probe_str = mission+'_'+probe
    prefix = mission+probe+'_'

    search_trs = micro_injection_load_search_time_range()
    nsearch_tr = n_elements(search_trs[*,0])
    mi_times = micro_injection_read_mi_times(mission_probe)
    min_dis = 2d
    pa_energy_range = [60d,300]    ; keV.
    fig_size = [12,4]
    common_time_step = project_info['common_time_step']

;---Collect all plot_trs.
    plot_trs = list()
    for sid=0,nsearch_tr-1 do begin
        time_range = search_trs[sid,*]+[-1,1]*constant('secofday')
if keyword_set(test) then time_range = time_double(['2015-09-01/00:00','2015-09-02/00:00'])
        r_var = lets_read('orbit', time_range, source=mission_probe, coord='sm')
        r_sm = get_var_data(r_var, times=times)
        dis = snorm(r_sm)
        dis_var = prefix+'dis'
        store_data, dis_var, times, dis
        add_setting, dis_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', 'Re', $
            'short_name', '|R|', $
            'yrange', [min_dis,15], $
            'ytickv', [5,10,15], $
            'yticks', 2, $
            'yminor', 5 )
        index = where(dis ge min_dis, count)
        orbit_trs = times[time_to_range(index, time_step=1)]
        norbit_tr = n_elements(orbit_trs[*,0])
        for oid=0,norbit_tr-1 do begin
            plot_tr = reform(orbit_trs[oid,*])
            index = where_pro(mi_times, '[]', plot_tr, count=count)
            if count eq 0 then continue
            plot_trs.add, plot_tr
        endfor
    endfor


;---Loop through plot_trs.
    foreach plot_tr, plot_trs do begin
        ; Get mi_times for this plot_tr.
        index = where_pro(mi_times, '[]', plot_tr, count=count)
        the_mi_times = mi_times[index]
        mi_trs = time_to_range(the_mi_times, time_step=common_time_step)
        nmi_tr = n_elements(mi_trs[*,0])
        time_range = plot_tr

    ;---Load data.
        ; Electron flux and related settings.
        flux_vars = micro_injection_read_cwt_kev_electron(plot_tr, mission_probe=mission_probe, errmsg=errmsg, log=1)
        spec_vars = flux_vars+'_mor'

        freq_range = [1,10.]*1e-3   ; Hz
        freq_range = minmax(1d/([1d,10]*60))
        significant_psd_levels = 1e-1*[1,significant_flux_levels[1]/significant_flux_levels[0]]
        significant_psd_levels = 1e-1*[1,1]

        energy_strs = list()
        search_info = orderedhash()
        foreach flux_var, flux_vars, vid do begin
            energy_str = get_var_setting(flux_var, 'labels')
            energy_strs.add, energy_str
            spec_var = spec_vars[vid]
            search_info[energy_str] = dictionary($
                'flux_var', flux_var, $
                'spec_var', spec_var, $
                'significant_flux_level', significant_flux_levels[vid], $
                'significant_psd_level', significant_psd_levels[vid], $
                'mi_trs', !null )
            options, [flux_var,spec_var], energy_str=energy_str
        endforeach
        search_settings = dictionary($
            'min_mi_dur', 600d, $
            'sample_energys', sample_energys, $
            'energy_strs', energy_strs, $
            'freq_range', freq_range )
    

        ; Electron PA.
        yrange = [5e4,5e10]
        log_ytickv = make_bins([5,10],1)
        ytickv = 10d^log_ytickv
        yticks = n_elements(ytickv)-1
        ytickn = '10!U'+string(log_ytickv,format='(I0)')
        ytickn[1:*:2] = ' '
        yminor = 9
        pa_var = mms_read_pa_spec_kev(time_range, probe=probe, energy_range=pa_energy_range, id='cdaweb', species='e')
        pa_spec = get_var_data(pa_var, times=times, pas)
        flux_90 = mean(pa_spec[*,5:6],dimension=2,nan=1)
        flux_0 = mean(pa_spec[*,[[0,1],[10,11]]],dimension=2,nan=1)
        pa_flux_var = var_store(prefix+'ele_kev_pa_flux', [[flux_0],[flux_90]], times)
        add_setting, pa_flux_var, smart=1, dictionary($
            'display_type', 'stack', $
            'labels', ['0','90']+' deg', $
            'colors', sgcolor(['red','blue']), $
            'ytitle', '(#/cm!U2!N-s-sr-keV)', $
            'yrange', yrange, $
            'ytickv', ytickv, $
            'yticks', yticks, $
            'yminor', yminor, $
            'ytickname', ytickn, $
            'ylog', 1, $
            'labflag', -1 )
        add_setting, pa_var, smart=1, dictionary($
            'display_type', 'spec', $
            'zrange', yrange, $
            'ztickv', ytickv, $
            'zticks', yticks, $
            'zminor', yminor, $
            'ztickname', ytickn, $
            'zlog', 1, $
            'constant', [45,90,135], $
            'yrange', [0,180], $
            'ytickv', [45,90,135], $
            'yticks', 2, $
            'yminor', 4 )
        aniso_var = var_store(prefix+'ele_kev_anisotropy', $
            flux_90/flux_0, times)
        add_setting, aniso_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', '#', $
            'short_name', 'A', $
            'ylog', 1, $
            'yrange', [0.05,50], $
            'ytickv', [0.1,1,10], $
            'yticks', 2, $
            'yminor', 9, $
            'constant', [0.1,1,10] )
        

        ; Orbit.
        routine = mission+'_read_orbit'
        r_gsm_var = call_function(routine, time_range, probe=probe, coord='gsm')
        options, r_gsm_var, requested_time_range=time_range
        mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
        foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'

        var_labels = prefix+['mlat','dis','mlt']
        nvar_label = n_elements(var_labels)
        options, prefix+'mlat', ytitle='MLat (deg)'
        options, prefix+'dis', ytitle='|R| (Re)'
        options, prefix+'mlt', ytitle='MLT (h)'
        var = prefix+'mlt'
        get_data, var, times, data
        index = where(data le 0, count)
        if count ne 0 then begin
            data[index] += 24
            store_data, var, times, data
        endif

        ; B field related vars.
        fillval = !values.f_nan
        field_time_range = time_range+[-1,1]*30.*60
        b_gsm_var = lets_read_this(func='mms_read_bfield', $
            field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
        if errmsg ne '' then begin
            b_gsm_var = prefix+'b_gsm'
            store_data, b_gsm_var, time_range, fltarr(2,3)+fillval
            add_setting, b_gsm_var, id='bfield', dictionary('coord', 'GSM' )
        endif
        e_gsm_var = lets_read_this(func='mms_read_efield', $
            field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
        if errmsg ne '' then begin
            e_gsm_var = prefix+'e_gsm'
            store_data, e_gsm_var, time_range, fltarr(2,3)+fillval
            add_setting, e_gsm_var, id='efield', dictionary('coord', 'GSM' )
        endif


        ; Model related vars.
        external_models = ['t89','t96','t01','t04s']
        external_models = ['t89']
        internal_models = ['dipole','igrf']


        ; Energy spec from CDAWeb.
        routine = mission+'_read_kev_electron_cdaweb'
        en_spec_var = call_function(routine, time_range, probe=probe)
        
        
        ; Ion velocity.
        u_gsm_var = mms_read_ion_vel(time_range, probe=probe, coord=coord, errmsg=errmsg)
        if errmsg ne '' then begin
            u_gsm_var = prefix+'u_gsm'
            store_data, u_gsm_var, time_range, fltarr(2,3)+fillval
            add_setting, u_gsm_var, id='velocity', dictionary('coord', 'GSM' )
        endif


        ; Low energy electron spec.
        ele_en_spec_var = lets_read_this(func='mms_read_en_spec_ele_cdaweb', $
            time_range, probe=mission_probe, id='thermal', errmsg=errmsg)
        if errmsg ne '' then begin
            ele_en_spec_var = lets_read_this(func='mms_read_en_spec_ele_cdaweb', get_name=1, $
                time_range, probe=mission_probe, id='thermal', errmsg=errmsg)
            foo_energys = [10,2.6e4]
            foo_nenergy = n_elements(foo_energys)
            foo_flux = fltarr(2,foo_nenergy)+fillval
            store_data, ele_en_spec_var, time_range, foo_flux, foo_energys
            unit = '#/cm!U2!N-s-sr-eV'
            add_setting, ele_en_spec_var, smart=1, {$
                requested_time_range: time_range, $
                display_type: 'spec', $
                unit: unit, $
                species: 'e', $
                species_name: 'e-', $
                ytitle: 'Energy!C(eV)', $
                subytitle: 'e-', $
                ylog: 1, $
                zlog: 1, $
                short_name: ''}
        endif
        options, ele_en_spec_var, $
            zrange=[1e3,1e8], zstyle=1, zlog=1, ztickv=[1e4,1e5,1e6,1e7], ztickname='10!U'+['4','5','6','7'], zticks=3, zminor=9, $
            yrange=[1.1e1,2.6e4], ystyle=1, ylog=1, ytickv=[1e2,1e3,1e4], ytickname='10!U'+['2','3','4'], yticks=2, yminor=9


    ;---Prepare for survey plot.

        plot_info = orderedhash()

        plot_info[en_spec_var] = dictionary($
            'routine', 'plot_spec', $
            'panel_label_text', 'e- EN', $
            'ypan', 1, $
            'setting', dictionary() )

        foreach energy_str, search_settings['energy_strs'] do begin
            my_info = search_info[energy_str]
            flux_var = my_info['flux_var']
            spec_var = my_info['spec_var']
            plot_info[flux_var] = dictionary($
                'routine', 'plot_line', $
                'panel_label_text', 'e- flux!C    '+energy_str, $
                'ypan', 0.6, $
                'setting', dictionary() )
            plot_info[spec_var] = dictionary($
                'routine', 'plot_spec', $
                'panel_label_text', 'Wavelet!C    '+energy_str, $
                'ypan', 1, $
                'setting', dictionary() )
        endforeach


        plot_info[pa_var] = dictionary($
            'routine', 'plot_spec', $
            'panel_label_text', 'e- PA', $
            'ypan', 1, $
            'setting', dictionary() )
        
        plot_info[pa_flux_var] = dictionary($
            'routine', 'plot_line', $
            'panel_label_text', 'e- PA flux', $
            'ypan', 0.8, $
            'setting', dictionary() )
        
        plot_info[aniso_var] = dictionary($
            'routine', 'plot_line', $
            'panel_label_text', 'e- aniso', $
            'ypan', 0.8, $
            'setting', dictionary() )


        plot_info[ele_en_spec_var] = dictionary($
            'routine', 'plot_spec', $
            'panel_label_text', 'Ele low', $
            'ypan', 1, $
            'setting', dictionary( ) )            

        plot_info[b_gsm_var] = dictionary($
            'routine', 'plot_linlog', $
            'panel_label_text', 'B GSM', $
            'ypans', 1.2, $
            'setting', dictionary($
                'linear_yrange', [-1,1]*20, $
                'log_yrange', [-1,1]*200, $
                'linear_tick_setting', dictionary($
                    'yticks', 2, $
                    'ytickv', [-1,0,1]*10, $
                    'yminor', 4 ), $
                'log_tick_setting', dictionary($
                    'yticks', 1, $
                    'ytickv', [1,10]*20, $
                    'yminor', 9 ) $
            ) $
        )
        
        plot_info[e_gsm_var] = dictionary($
            'routine', 'plot_linlog', $
            'panel_label_text', 'E GSM', $
            'ypans', 1.2, $
            'setting', dictionary($
                'linear_yrange', [-1,1]*5, $
                'log_yrange', [-1,1]*500, $
                'linear_tick_setting', dictionary($
                    'yticks', 2, $
                    'ytickv', [-1,0,1]*5, $
                    'yminor', 5 ), $
                'log_tick_setting', dictionary($
                    'yticks', 1, $
                    'ytickv', [10,100]*5, $
                    'yminor', 9 ) $
            ) $
        )
        
        plot_info[u_gsm_var] = dictionary($
            'routine', 'plot_linlog', $
            'panel_label_text', 'U GSM', $
            'ypans', 1.2, $
            'setting', dictionary($
                'linear_yrange', [-1,1]*60, $
                'log_yrange', [-1,1]*600, $
                'linear_tick_setting', dictionary($
                    'yticks', 2, $
                    'ytickv', [-1,0,1]*40, $
                    'yminor', 2 ), $
                'log_tick_setting', dictionary($
                    'yticks', 1, $
                    'ytickv', [10,100]*6, $
                    'yminor', 9 ) $
            ) $
        )

        plot_vars = plot_info.keys()
        nplot_var = n_elements(plot_vars)
        
        ; Settings for en_spec_var.
        zrange = [1e-1,1e5]
        log_ztickv = smkarthm(-1,5,1,'dx')
        ztickv = 10d^log_ztickv
        zticks = n_elements(ztickv)-1
        ztickn = '10!U'+string(log_ztickv,format='(I0)')
        ztickn[1:*:2] = ' '
        zminor = 9
        zticklen = -0.5
        yrange = [10,100]*5
        ytickv = yrange
        yticks = n_elements(ytickv)-1
        yminor = 9
        add_setting, en_spec_var, smart=1, dictionary($
            'yrange', yrange, $
            'ytickv', ytickv, $
            'yticks', yticks, $
            'yminor', yminor )
        options, en_spec_var, $
            ytitle='Energy!C(keV)', ysubtitle=' ', $
            zrange=zrange, ztickv=ztickv, zticks=zticks, zminor=zminor, ztickname=ztickn, zticklen=zticklen, $
            constant=sample_energys

        ; Settings for flux_vars.
        yrange = alog10(zrange)
        ystep = 2
        ytickv = smkarthm(yrange[0],yrange[1],ystep,'dx')
        yticks = n_elements(ytickv)-1
        options, flux_vars, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=ystep
        foreach flux_var, flux_vars, vid do begin
            my_info = search_info[energy_strs[vid]]
            significant_flux_level = my_info['significant_flux_level']
            options, flux_var, constant=significant_flux_level
        endforeach

        ; Settings for psd spec_vars.
        ; Change yticks from Hz to mHz.
        spec_plot_vars = spec_vars+'_plot'
        options, spec_plot_vars, zticklen=zticklen
        foreach spec_var, spec_vars, vid do begin
            specs = get_var_data(spec_var, times=times, freqs, settings=settings)
            settings['ytitle'] = 'Freq!C(mHz)'
            settings['yrange'] = settings['yrange']*1e3
            flux_var = flux_vars[vid]
            energy_str = get_var_setting(flux_var, 'labels')
            settings['ztitle'] = 'PSD of Log flux!C'+energy_str
            settings['constant'] = freq_range*1e3
            spec_plot_var = var_store(spec_plot_vars[vid], specs, times, freqs*1e3, settings=settings)
            pid = where(plot_vars eq spec_var)
            plot_vars[pid] = spec_plot_var
            info = plot_info[spec_var]
            plot_info.remove, spec_var
            plot_info[spec_plot_var] = info
        endforeach


    ;---Survey plot.
        base = 'micro_injection_mi_times_survey_plot_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_'+mission+'_'+probe+'_'+version+'.pdf'
        year_str = time_string(plot_tr[0],tformat='YYYY')
        plot_file = join_path([plot_root,'mi_times',year_str,base])
        if keyword_set(test) then plot_file = 0

        margins = [15,3.5+nvar_label,9,2]
        ypans = []
        foreach plot_var, plot_vars, pid do begin
            info = plot_info[plot_var]
            if ~info.haskey('ypan') then info['ypan'] = 1
            ypans = [ypans,info['ypan']]
        endforeach
        pansize = [12,0.8]
        plot_poss = panel_pos(plot_file, nypan=nplot_var, fig_size=fig_size, ypans=ypans, pansize=pansize, margins=margins)
        
        ; Use positions to determine [x,y]ticklen.
        abs_ticklen = 0.3
        foreach plot_var, plot_vars, pid do begin
            my_info = plot_info[plot_var]
            my_info['position'] = plot_poss[*,pid]
            my_info['abs_ticklen'] = abs_ticklen
            ;plot_info[plot_var] = my_info
        endforeach

        sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

        ; Default settings.
        panel_letters = letters(nplot_var)
        foreach plot_var, plot_vars, pid do begin
            my_info = plot_info[plot_var]
            if ~my_info.haskey('ypan') then my_info['ypan'] = 1d
            if ~my_info.haskey('panel_label_text') then my_info['panel_label_text'] = ' '
            if ~my_info.haskey('panel_letter') then my_info['panel_letter'] = panel_letters[pid]
            if ~my_info.haskey('panel_label_msg') then my_info['panel_label_msg'] = my_info['panel_letter']+') '+my_info['panel_label_text']
        endforeach

        tickinterval = 3600*2d
        tplot_options, 'tickinterval', tickinterval
        foreach plot_var, plot_vars, pid do begin
            my_info = plot_info[plot_var]
        
            my_pos = my_info['position']
            plot_routine = my_info['routine']
            plot_setting = my_info['setting']
            plot_setting['position'] = my_pos
            plot_setting['noerase'] = (pid eq 0)? 0: 1
            plot_setting['xtickformat'] = (pid eq nplot_var-1)? '': '(A1)'
            plot_setting['novtitle'] = (pid eq nplot_var-1)? 0: 1
            plot_setting['time_range'] = plot_tr
            plot_setting['tickinterval'] = tickinterval
            plot_setting['panel_label_pos'] = [xchsz*1,my_pos[3]-ychsz*0.7]
            plot_setting['panel_label_msg'] = my_info['panel_label_msg']
            plot_setting['var_labels'] = var_labels
            plot_setting['vlab_margin'] = margins[0]-1
;if plot_var eq aniso_var then stop
            plot_setting = plot_setting.tostruct()        
            tmp = call_function(plot_routine, plot_var, _extra=plot_setting)
        endforeach

        
        
        ; Add significant flux level.
        if keyword_set(test) eq 1 then thick = 2 else thick = 10
        cbpos = get_var_setting(en_spec_var,'zposition')
        zrange = get_var_setting(en_spec_var,'zrange')
        zlog = 1
        set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
        sig_flux_color = sgcolor('red')
        foreach ty, significant_flux_levels do begin
            plots, [0,1], 10d^ty+[0,0], color=sig_flux_color, thick=thick, data=1
        endforeach

        ; Add contours and microinjection time ranges.
        contour_color = sgcolor('brown')
        foreach spec_var, spec_plot_vars, vid do begin
            pid = where(plot_vars eq spec_var, count)
            if count eq 0 then continue
            my_info = search_info[energy_strs[vid]]
            significant_psd_level = my_info['significant_psd_level']
            tpos = plot_poss[*,pid]
            specs = get_var_data(spec_var, freqs, times=times, settings=settings)
            levels = [significant_psd_level]
            yrange = settings['yrange']

            contour, specs, times, freqs, position=tpos, $
                xlog=0, xstyle=5, xrange=xrange, $
                ylog=1, ystyle=5, yrange=yrange, levels=levels, noerase=1, color=contour_color

            ; Add significant level.
            cbpos = get_var_setting(spec_var,'zposition')
            zrange = get_var_setting(spec_var,'zrange')
            zlog = 1
            set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
            ty = significant_psd_levels[vid]
            plots, [0,1], ty+[0,0], color=contour_color, thick=thick, data=1

        endforeach
        
        ; Add mi_times.
        if keyword_set(test) eq 1 then thick = 4 else thick = 20
        pid = nplot_var-1
        tpos = reform(plot_poss[*,pid])
        yrange = [0d,1]
        set_axis, position=tpos, xrange=plot_tr, yrange=yrange
        for ii=0,nmi_tr-1 do begin
            plots, reform(mi_trs[ii,*]), yrange[1]+[0,0], thick=thick, $
                color=sgcolor('red'), data=1
        endfor


        if keyword_set(test) then stop
        sgclose

    endforeach


end


mission_probe = ['mms','4']
test = 0
micro_injection_gen_mi_time_survey_plot_v02, mission_probe, test=test
end