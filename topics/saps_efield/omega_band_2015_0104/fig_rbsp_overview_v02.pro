;+
;-

function fig_rbsp_overview_v02, test=test, event_info=event_info

    version = 'v02'
    id = '2015_0104'

    if n_elements(event_info) eq 0 then event_info = saps_efield_load_data(id)
    time_range = time_double(['2015-01-04/10:50','2015-01-04/12:50'])
    time_range = time_double(['2015-01-04/08:30','2015-01-04/14:30'])
    data_time_range = time_double(['2015-01-04/06:00','2015-01-04/16:00'])
    plot_time_range = time_range
    probes = ['b']
    
    ob_times = time_double('2015-01-04/'+['11:52','12:19:30','12:40','13:04:30'])
    ob_color = sgcolor('red')

    fc_text_time = plot_time_range[0]
    
    
    load_info = dictionary()
    update = 0

    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'

;        o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
;        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=update)

        eu_var = prefix+'e_uvw_comp1'
        spec_var = prefix+'e_uvw_comp1_mor'
        if check_if_update(spec_var) then begin
            rbsp_efw_phasef_read_e_uvw, data_time_range, probe=probe
            vars = stplot_split(prefix+'e_uvw')
            field_var = vars[0]
            time_step = 1d/32
            common_times = make_bins(data_time_range, time_step)
            interp_time, field_var, common_times
            scale_info = {s0:time_step*2, s1:16, dj:1d/8, ns:0d }
            spec_var = stplot_mor_new(field_var, scale_info=scale_info)
            data = get_var_data(spec_var, freqs, times=times, limits=lim)
            freq_range = [0.1,10]
            freq_index = where_pro(freqs, '[]', freq_range, count=nfreq)
            if nfreq eq 0 then begin
                errmsg = 'Invalid frequencies ...'
                return, retval
            endif
            ; Convert to psd spectrogram in X^2/Hz.
            data = data[*,freq_index]
            freqs = freqs[freq_index]
            cwt_info = get_setting(spec_var, 'cwt_info')
            c_unit = 2*cwt_info.c_tau*cwt_info.dt/cwt_info.cdelta
            data *= c_unit  ; convert unit from (X)^2 to (X)^2/Hz
            store_data, spec_var, times, data, freqs
        endif
        
        
        
        e_spec_var = prefix+'e_spec2'
        out_var = e_spec_var

        spec_vars = list()
        spec_vars.add, rbsp_read_wave_spec_mhz(data_time_range, probe=probe)
        spec_vars.add, rbsp_read_wave_spec_khz(data_time_range, probe=probe, id='e')
        spec_vars.add, spec_var
        freq_ranges = [[1e4,5e4],[10,1e4],[0.1,10]]
        
        
        yrange = [0.1,5e4]
        dfreq = 1.15
        freqs = smkgmtrc(yrange[0],yrange[1],dfreq, 'dx')
        nfreq = n_elements(freqs)
        xrange = time_range
        time_step = 1d
        common_times = make_bins(xrange,time_step)
        ntime = n_elements(common_times)
        specs = fltarr(ntime,nfreq)
        foreach var, spec_vars, var_id do begin
            if var eq '' then continue
            data = get_var_data(var, vals, at=common_times, limits=lim)
            freq_range = freq_ranges[*,var_id]
            index = where_pro(freqs, '[)', freq_range, count=count)
            if count eq 0 then message, 'Inconsistency ...'
            specs[*,index] = transpose(sinterpol(transpose(data),vals,freqs[index]))
        endforeach
        unit = lim.unit
        zrange = [1e-8,1e1]

        log_ytickv = make_bins(minmax(alog10(yrange)),1,inner=1)
        ytickv = 10d^log_ytickv
        ytickname = get_short_log_tickname(log_ytickv)
        yminor = 9

        store_data, out_var, common_times, specs, freqs
        add_setting, out_var, smart=1, dictionary($
            'requested_time_range', time_range, $
            'no_interp', 1, $
            'display_type', 'spec', $
            'unit', unit, $
            'ytitle', 'Freq (Hz)', $
            'yrange', yrange, $
            'ylog', 1, $
            'ytickv', ytickv, $
            'ytickname', ytickname, $
            'yminor', yminor, $
            'zlog', 1, $
            'zrange', zrange, $
            'short_name', 'E' )
        
        
        
           
        
        ; Cyclotron freqs.
        fc_vars = list()
        foreach species, ['e','o','he','p'] do fc_vars.add, rbsp_read_gyro_freq(time_range, probe=probe, species=species, get_name=1)
        options, prefix+'fce', labels='f!Dc,e'
        options, prefix+'fcp', labels='f!Dc,H'
        options, prefix+'fco', labels='f!Dc,O'
        options, prefix+'fche', labels='f!Dc,He'

        var = prefix+'flh'
        fcp = get_var_data(prefix+'fcp', times=times)
        store_data, var, times, fcp*43, limits={labels:'f!DLH!N'}
        fc_vars.add, var
        fc_vars = fc_vars.toarray()
        fc_colors = get_color(n_elements(fc_vars))
        foreach var, fc_vars, ii do options, var, 'colors', fc_colors[ii]

        e_spec_combo = e_spec_var+'_combo'
        store_data, e_spec_combo, data=[e_spec_var,fc_vars]
        options, e_spec_combo, 'yrange', get_setting(e_spec_var,'yrange')
        options, e_spec_combo, 'labels', ''
        options, e_spec_combo, 'ztitle', 'E ((mV/m)!U2!N/Hz)'


        b_spec_var = prefix+'b_spec_hz'
        b_spec_combo = b_spec_var+'_combo'
        store_data, b_spec_combo, data=[b_spec_var,fc_vars]
        options, b_spec_combo, 'yrange', get_setting(b_spec_var,'yrange')
        options, b_spec_combo, 'labels', ''
        

        ; Model field and tilt angle.
        b_vars = list()
        b_gsm_var = prefix+'b_gsm'
        b_vars.add, b_gsm_var


        external_models = ['t89','t96','t01','t04s']
        internal_model = 'igrf'
        bmod_vars = list()
        foreach external_model, external_models do begin
            suffix = '_'+internal_model+'_'+external_model
            bmod_var = prefix+'bmod_gsm'+suffix
            bmod_vars.add, bmod_var
            b_vars.add, bmod_var
        endforeach
        
        bmag_vars = list()
        foreach var, b_vars do begin
            get_data, var, times, b_vec, limits=lim
            bmag = snorm(b_vec)
            bmag_var = var+'_mag'
            store_data, bmag_var, times, bmag
            add_setting, bmag_var, smart=1, dictionary($
                'display_type', 'scalar', $
                'short_name', '|B|', $
                'unit', 'nT' )
            bmag_vars.add, bmag_var
        endforeach
        bmag_combo_var = prefix+'bmag_combo'
        times = get_var_time(b_gsm_var)
        foreach var, bmag_vars do interp_time, var, times
        bmag_combo_var = stplot_merge(bmag_vars[1:*], output=bmag_combo_var)
        add_setting, bmag_combo_var, smart=1, dictionary($
            'display_type', 'stack', $
            'short_name', '|B|!Dmod!N-|B|!Dobs', $
            'unit', 'nT', $
            'constant', 0, $
            'labels', strupcase(external_models), $
            'colors', sgcolor(['blue','orange','green','purple']))
        bmag_var = bmag_vars[0]
        options, bmag_var, 'ytitle', '|B| (nT)'
        bmag_obs = get_var_data(bmag_var, times=times)
        bmag_mod = get_var_data(bmag_combo_var, times=times)
        foreach tmp, external_models, ii do bmag_mod[*,ii] -= bmag_obs
        store_data, bmag_combo_var, times, bmag_mod
        set_ytick, bmag_combo_var, yrange=[-25,65], ytickv=[-25,0,25,50], yminor=5
        

        b_theta_vars = list()
        foreach var, b_vars do begin
            get_data, var, times, b_vec, limits=lim
            coord = strlowcase(lim.coord)
            if coord ne 'sm' then begin
                b_vec = cotran(b_vec, times, coord+'2sm')
            endif
            theta_var = var+'_theta'
            theta = atan(b_vec[*,2]/snorm(b_vec))*constant('deg')
            store_data, theta_var , times, theta
            add_setting, theta_var, smart=1, dictionary($
                'display_type', 'scalar', $
                'short_name', 'Tilt', $
                'unit', 'deg' )
            if var eq b_gsm_var then begin
                b_sm = b_vec
                b_sm_var = prefix+'b_sm'
                lim.coord = 'sm'
                store_data, b_sm_var, times, b_sm, limits=lim
            endif
            b_theta_vars.add, theta_var
        endforeach

        theta_combo_var = prefix+'b_theta_combo'
        times = get_var_time(b_gsm_var)
        foreach var, b_theta_vars do interp_time, var, times
        theta_combo_var = stplot_merge(b_theta_vars, output=theta_combo_var)
        add_setting, theta_combo_var, smart=1, dictionary($
            'display_type', 'stack', $
            'short_name', tex2str('theta'), $
            'unit', 'deg', $
            'constant', 0, $
            'labels', ['Obs',strupcase(external_models)], $
            'colors', sgcolor(['black','blue','orange','green','purple']) )
        set_ytick, theta_combo_var, yrange=[30,45], ytickv=[30,40], yminor=5
;        theta_combo_var = stplot_merge(b_theta_vars[1:*], output=theta_combo_var)
;        add_setting, theta_combo_var, smart=1, dictionary($
;            'display_type', 'stack', $
;            'short_name', tex2str('theta'), $
;            'unit', 'deg', $
;            'constant', 0, $
;            'labels', strupcase(external_models), $
;            'colors', sgcolor(['blue','orange','green','purple']) )
;        b_theta_var = b_theta_vars[0]
;        theta_obs = get_var_data(b_theta_var, times=times)
;        theta_mod = get_var_data(theta_combo_var, times=times)
;        foreach tmp, external_models, ii do theta_mod[*,ii] -= theta_obs
;        store_data, theta_combo_var, times, theta_mod
;        set_ytick, theta_combo_var, yrange=[-4,12], ytickv=[0,5,10], yminor=5
        
        
        
        ; magnetic pressure.
        pmag_var = rbsp_read_magnetic_pressure(b_var=b_gsm_var)
        options, pmag_var, 'ylog', 1
        ; thermal_pressure.
        pthermal_var = prefix+'p_thermal'
        foreach species, ['e','p','o'] do begin
            vars = rbsp_read_hope_moments(time_range, probe=probe, species=species)
        endforeach
        
        species = 'e'
        p_thermal = calc_thermal_pressure($
            get_var_data(prefix+species+'_n',times=times), $
            get_var_data(prefix+species+'_t',at=times))
        foreach species, ['p','o'] do begin
            p_thermal += calc_thermal_pressure($
                get_var_data(prefix+species+'_n',at=times), $
                get_var_data(prefix+species+'_t',at=times))
        endforeach
        store_data, pthermal_var, times, p_thermal
        add_setting, pthermal_var, smart=1, dictionary($
            'ylog', 1, $
            'display_type', 'scalar', $
            'short_name', 'P!Dth', $
            'unit', 'nPa' )
        
        
        p_thermal = get_var_data(pthermal_var, times=times)
        p_mag = get_var_data(pmag_var, at=times)
        beta = p_thermal/p_mag
        beta_var = prefix+'beta'
        store_data, beta_var, times, beta
        add_setting, beta_var, smart=1, dictionary($
            'ylog', 1, $
            'yrange', [1e-4,1e-1]*7, $
            'constant', [1e-2,1e-1], $
            'display_type', 'scalar', $
            'short_name', tex2str('beta'), $
            'unit', '#' )
        
        
        ; FMLat.
        fmlat_north_var = prefix+'fmlat'+suffix+'_north'
        fmlat_south_var = prefix+'fmlat'+suffix+'_south'
        fmlat_var = prefix+'fmlat'
        fmlat_var = stplot_merge([fmlat_north_var,fmlat_south_var], output=fmlat_var)
        get_data, fmlat_var, times, data
        store_data, fmlat_var, times, abs(data)
        add_setting, fmlat_var, smart=1, dictionary($
            'display_type', 'stack', $
            'ytitle', '(deg)', $
            'yrange', [40,70], $
            'ytickv', [50,60], $
            'yticks', 1, $
            'yminor', 2, $
            'constant', [50,60], $
            'labels', ['North','South'], $
            'colors', sgcolor(['red','blue']) )



        ; dis.
        r_gsm_var = prefix+'r_gsm'
        dis_var = prefix+'dis'
        get_data, r_gsm_var, times, data
        store_data, dis_var, times, snorm(data)
        add_setting, dis_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', '|R|', $
            'unit', 'Re' )
        var = dis_var
        options, var, 'yrange', [1,6]
        options, var, 'ytickv', [2,4,6]
        options, var, 'yticks', 2
        options, var, 'yminor', 2
        options, var, 'constant', [2,4]


        ; MLat.
        mlat_var = rbsp_read_mlat(time_range, probe=probe)
        options, mlat_var, 'yrange', [-1,1]*20
        options, mlat_var, 'constant', [-1,0,1]*10
        options, mlat_var, 'ytickv', [-1,0,1]*10
        options, mlat_var, 'yticks', 2
        options, mlat_var, 'yminor', 2


        ; MLT.
        mlt_var = prefix+'mlt'
        yrange = [-1,1]*12
        ytickv = [-1,0,1]*6
        yticks = n_elements(ytickv)-1
        yminor = 3
        constants = [-1,0,1]*6
        options, mlt_var, 'yrange', yrange
        options, mlt_var, 'constant', constants
        options, mlt_var, 'ytickv', ytickv
        options, mlt_var, 'yticks', yticks
        options, mlt_var, 'yminor', yminor

;        ; Tilt.
;        var = bmod_gsm_var+'_theta'
;        options, var, 'yrange', [0,90]
;        options, var, 'ytickv', [40,80]
;        options, var, 'constant', [40,80]
;        options, var, 'yticks', 1
;        options, var, 'yminor', 4
;        var = dtheta_var
;        yrange = [-45,10]
;        ystep = 20
;        ytickv = make_bins(yrange,ystep, inner=1)
;        yticks = n_elements(ytickv)-1
;        options, var, 'ytickv', ytickv
;        options, var, 'yticks', yticks
;        options, var, 'yrange', yrange
;        options, var, 'yminor', 2
;        options, var, 'constant', ytickv

        ; E
        e_var = prefix+'edot0_fac'
        e_vars = stplot_split(e_var)
        e_vars = e_vars[1:2]
        var = e_vars
        
;        foreach var, e_vars, ii do begin
;            get_data, prefix+'edot0_fac', times, e_srvy
;            e_sf = get_var_data(prefix+'edot0_spinfit_fac', at=times)
;            data = [[e_srvy[*,ii+1]],[e_sf[*,ii+1]]]
;            store_data, var, times, data
;            options, var, 'labels', ['survey','spinfit']
;            options, var, 'colors', sgcolor(['red','green'])
;            options, var, 'labflag', -1
;            
;;            index = where_pro(times, '[]', time_double(['2015-01-04/12:47','2015-01-04/12:51']), count=count)
;;            if count ne 0 then begin
;;                data[index,*] = !values.f_nan
;;                store_data, var, times, data
;;            endif
;            stop
;        endforeach
        
        

        
        
        dens_var = prefix+'density'
        options, dens_var, 'constant', [1,10,100,1000]
        options, dens_var, 'yrange', [0.5,3e3]
        e_en_var = prefix+'e_en_spec'
        p_en_var = prefix+'p_en_spec'
        o_en_var = prefix+'o_en_spec'
        b_theta_var = prefix+'b_gsm_theta'
        b_var = prefix+'b_sm'
        pthermal_var = prefix+'p_thermal'
        
        
        
        
        model_suffix = '_igrf_t04s_north'
        model_suffix2 = '_igrf_t96_north'
        mlt_image_var = 'thg_asf_mlt_image_rect'
        mlt_images = get_var_data(mlt_image_var, times=times, settings=settings)
        ntime = n_elements(times)
        
        ; KEO
        keo_var = 'thg_asf_keo'
        ytitle = 'MLat!C(deg)'

        mlat_range = [60,67]
        mlt_range = [1,4]
        mlt_bins = settings['mlt_bins']
        mlat_bins = settings['mlat_bins']
        index = where_pro(mlt_bins, '[]', mlt_range, count=nbin)
        ;bins = mlt_bins[index]
        ;keo = fltarr(ntime,nbin)
        keo = total(mlt_images[*,index,*],2)/nbin
        
        fmlt_var = prefix+'fmlt'+model_suffix
        fmlt = get_var_data(fmlt_var, at=times)
        fmlt_del = [-1,1]*0.2
        nbin = n_elements(mlat_bins)
        keo = fltarr(ntime,nbin)
        foreach time, times, time_id do begin
            index = where_pro(mlt_bins,'[]',fmlt[time_id]+fmlt_del,count=count)
            if count eq 0 then continue
            keo[time_id,*] = total(mlt_images[time_id,index,*],2)/count
        endforeach

        store_data, keo_var, times, keo, mlat_bins, limits={ytitle:ytitle,spec:1,zlog:0,color_table:49}
        set_ytick, keo_var, yrange=[58,72], ytickv=[60,65,70], yminor=5
        zlim, keo_var, 0,3000, 0
        options, keo_var, ztickv=[0,1,2]*1e3, zticklen=-0.5, $
            zticks=2, zminor=2, ztitle='Count (#)'
        
        ewo_var = 'thg_asf_ewo'
        index = where_pro(mlat_bins, '[]', mlat_range, count=nbin)
        ewo = total(mlt_images[*,*,index],3)/nbin
        
        fmlat_var = prefix+'fmlat'+model_suffix
        fmlat = abs(get_var_data(fmlat_var, at=times))
        fmlat_del1 = [0.5,3]
        fmlat_del2 = [-1,-0.5]
        fmlat_del = [-0.5,3]
        nbin = n_elements(mlt_bins)
        ewo = fltarr(ntime,nbin)
        foreach time, times, time_id do begin
            index1 = where_pro(mlat_bins,'[]',fmlat[time_id]+fmlat_del1,count=count1)
            index2 = where_pro(mlat_bins,'[]',fmlat[time_id]+fmlat_del2,count=count2)
            index = [index1,index2]
            count = count1+count2
            index = where_pro(mlat_bins,'[]',fmlat[time_id]+fmlat_del,count=count)
            if count eq 0 then continue
            ewo[time_id,*]  = total(mlt_images[time_id,*,index],3)/count
        endforeach
        
        store_data, ewo_var, times, ewo, mlt_bins, limits={ytitle:'MLT (h)', spec:1, color_table:49}
        set_ytick, ewo_var, yrange=[0.8,2.8], ytickv=[1,2], yminor=5
        zlim, ewo_var, 0,3000, 0
        options, ewo_var, ztickv=[0,1,2]*1e3, zticklen=-0.5, $
            zticks=2, zminor=2, ztitle='Count (#)', ytitle='MLT!C(h)'
        
        db_var = prefix+'b1_fac'
        update = 1
        

        e_sv = get_var_data(prefix+'edot0_fac', times=times, in=plot_time_range)
        copy_data, prefix+'edot0_spinfit_fac', prefix+'edot0_spinfit_fac_copy'
        interp_time, prefix+'edot0_spinfit_fac_copy', times
        e_sf = get_var_data(prefix+'edot0_spinfit_fac_copy')
        foreach var, e_vars, ii do begin
            data = [[e_sv[*,ii+1]],[e_sf[*,ii+1]]]
            store_data, var, times, data, limits={$
                colors:sgcolor(['black','red']),labels:['survey','spinfit']}
        endforeach
        var = e_vars
        options, var, 'yrange', [-1,1]*5
        options, var, 'ytickv', [-1,0,1]*3
        options, var, 'yticks', 2
        options, var, 'yminor', 3
        options, var, 'ytitle', '(mV/m)'
        options, var, 'constant', [-1,0,1]*3
        options, var, 'labflag', -1

        e_spinfit_var = prefix+'edot0_spinfit_fac'
        var = e_spinfit_var
        set_ytick, var, yrange=[-1,1]*3.5, ytickv=[-1,0,1]*2, yminor=2
        
        plot_vars = [theta_combo_var,bmag_combo_var,dens_var,e_en_var,p_en_var,$
            o_en_var,e_spec_combo,keo_var,eu_var,e_spinfit_var]
        nvar = n_elements(plot_vars)
        fig_letters = letters(nvar)
        fac_labels = ['||',tex2str('perp')+','+['out','west']]
        fig_labels = fig_letters+') '+['Tilt','|B|!DObs-Mod!N','N','e-','H+','O+','E spec','Keo','Eu','E!Dspinfit']
        ypans = fltarr(nvar)+1.
        index = where(plot_vars eq e_spec_combo, count)
        if count ne 0 then ypans[index] = 1.5
        index = where(plot_vars eq theta_combo_var, count)
        if count ne 0 then ypans[index] = 1.1
        index = where(plot_vars eq bmag_combo_var, count)
        if count ne 0 then ypans[index] = 0.9
        index = where(plot_vars eq beta_var, count)
        if count ne 0 then ypans[index] = 0.8
        index = where(plot_vars eq mlat_var, count)
        if count ne 0 then ypans[index] = 0.6
        index = where(plot_vars eq mlt_var, count)
        if count ne 0 then ypans[index] = 0.6
        index = where(plot_vars eq dis_var, count)
        if count ne 0 then ypans[index] = 0.6
        index = where(plot_vars eq fmlat_var, count)
        if count ne 0 then ypans[index] = 0.9
        index = where(plot_vars eq keo_var, count)
        if count ne 0 then ypans[index] = 1.3
        index = where(plot_vars eq ewo_var, count)
        if count ne 0 then ypans[index] = 1.3
                
        the_vars = [e_en_var,p_en_var,o_en_var]
        options, the_vars, 'ytitle', 'Energy!C(eV)'
        options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
            ytickv=[100,1000,10000], yticks=2, yminor=9, $
            yrange=[20,5.1e4]

        the_vars = e_en_var 
        ztickn = '10!U'+['4','5','6','7','8','9','10']
        ztickn[1:*:2] = ' '       
        options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
            ytickv=[100,1000,10000], yticks=2, yminor=9, $
            ztickv=10d^[4,5,6,7,8,9,10], zticks=6, ztickname=ztickn
            
        options, b_spec_var, ztickv=10d^[-3,-2,-1,0,1], zticks=4, $
            ztickname=['10!U-3','10!U-2','0.1','1','10']

        zticklen = -0.5
        the_vars = [e_en_var,p_en_var,o_en_var,b_spec_var,e_spec_var]
        options, the_vars, zticklen=zticklen, zminor=9
        
        
        
        the_vars = [b_spec_combo,e_spec_combo]
        options, the_vars, 'ytitle', 'Freq!C(Hz)'
        
        options, mlt_var, yrange=[-1,1]*7
        options, mlat_var, yrange=[-1,1]*15
        
        pansize = [5,0.8]
        margins = [12,7,8,1]
        
        if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
        plot_file = join_path([plot_dir,'fig_rbsp'+probe+'_overview_'+id+'_'+version+'.pdf'])
        print, plot_file
        if keyword_set(test) then plot_file = 0
        poss = panel_pos(plot_file, ypans=ypans, fig_size=fig_size, $
            nypan=nvar, pansize=pansize, panid=[0,1], margins=margins)
        sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, inch=1

        xticklen_chsz = -0.25   ; in ychsz.
        yticklen_chsz = -0.40   ; in xchsz.
        for ii=0,nvar-1 do begin
            tpos = poss[*,ii]
            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
            options, plot_vars[ii], 'xticklen', xticklen
            options, plot_vars[ii], 'yticklen', yticklen
            
            spec = get_setting(plot_vars[ii],'spec',exist)
            if exist and spec eq 1 or plot_vars[ii] eq e_spec_combo then begin
                cbpos = tpos
                cbpos[0] = cbpos[2]+xchsz*0.8
                cbpos[2] = cbpos[0]+xchsz*0.8
                options, plot_vars[ii], 'zposition', cbpos
            endif
        endfor

        lshell_var = prefix+'lshell'
        var_labels = [mlt_var,lshell_var,bmag_var]
        options, mlt_var, 'ytitle', 'MLT (h)'
        options, lshell_var, 'ytitle', 'L (#)'
        
        
      
        ;margins = [12,6,10,1]
        ;poss = sgcalcpos(nvar, margins=margins, ypans=ypans)
        tplot_options, 'tickinterval', 30*60
        tplot, plot_vars, trange=plot_time_range, var_label=var_labels, position=poss, vlab_margin=10, noerase=1
        for pid=0,nvar-1 do begin
            tpos = poss[*,pid]
            tx = tpos[0]-xchsz*10
            ty = tpos[3]-ychsz*0.7
            msg = fig_labels[pid]
            xyouts, tx,ty,msg, normal=1
        endfor
        
;        pid = 0
;        tpos = poss[*,pid]
;        tx = tpos[0]+xchsz*0.5
;        ty = tpos[3]-ychsz*1
;        msg = strupcase('rbsp-'+probe)
;        xyouts, tx,ty,msg, normal=1;, color=sgcolor('white')
        
        
        ;bar_times = make_bins(plot_time_range, 600, inner=1)
        ;timebar, bar_times, linestyle=1, color=sgcolor('silver')
        
        
        timebar, ob_times, linestyle=2, color=ob_color
        


        foreach spec_var, [e_spec_combo] do begin
            pid = where(plot_vars eq spec_var, count)
            if count eq 0 then continue
            tpos = poss[*,pid]
            spec_settings = get_var_setting(spec_var)
            xrange = plot_time_range
            yrange = spec_settings['yrange']
            ylog = 1
            plot, xrange, yrange, $
                ylog=ylog, xrange=xrange, yrange=yrange, $
                xstyle=5, ystyle=5, $
                position=tpos, noerase=1, nodata=1
            
            nfreq_var = n_elements(fc_vars)
            colors = get_color(nfreq_var)
            foreach var, fc_vars, var_id do begin
                color = colors[var_id]
                yys = get_var_data(var, times=xxs, settings=settings)
                tx = fc_text_time
                ty = interpol(yys, xxs, tx)
                if product(ty-yrange) ge 0 then continue

                msg = settings.labels
                color = settings.colors
                oplot, xxs, yys, color=color
                tmp = convert_coord(tx,ty, data=1, to_normal=1)
                tx = tmp[0]-xchsz*1.5
                ty = tmp[1]-ychsz*0.35
                tx = tmp[0]+xchsz*1
                ty = tmp[1]+ychsz*0.35
                xyouts, tx,ty,normal=1, msg, charsize=label_size, color=color
            endforeach
        endforeach
        
        
        pid = where(plot_vars eq ewo_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            yrange = get_setting(ewo_var, 'yrange')
            set_axis, position=tpos, xrange=plot_time_range, yrange=yrange
            fmlt_var = prefix+'fmlt'+model_suffix
            fmlts = get_var_data(fmlt_var, times=times, in=plot_time_range)
            oplot, times, fmlts, color=sgcolor('red')
            tmp = convert_coord(times[0],fmlts[0], data=1, to_normal=1)
            tx = tmp[0]+ychsz*0.5
            ty = tmp[1]+ychsz*0.5
            msg = strupcase('rbsp-'+probe)
            xyouts, tx,ty,msg, normal=1, color=sgcolor('red')
        endif
        
        
        pid = where(plot_vars eq keo_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            yrange = get_setting(keo_var, 'yrange')
            set_axis, position=tpos, xrange=plot_time_range, yrange=yrange
            fmlat_var = prefix+'fmlat'+model_suffix
            fmlats = abs(get_var_data(fmlat_var, times=times, in=plot_time_range))
            oplot, times, fmlats, color=sgcolor('red')
            tmp = convert_coord(times[0],fmlats[0], data=1, to_normal=1)
            tx = tmp[0]+ychsz*0.5
            ty = tmp[1]+ychsz*0.5
            msg = strupcase('rbsp-'+probe)
            xyouts, tx,ty,msg, normal=1, color=sgcolor('red')
            fmlat_var = prefix+'fmlat'+model_suffix2
            fmlats = abs(get_var_data(fmlat_var, times=times, in=plot_time_range))
            oplot, times, fmlats, color=sgcolor('red')
        endif
        
        
        if keyword_set(test) then stop
        sgclose
    endforeach
    
    return, plot_file

end

print, fig_rbsp_overview_v02(test=0, event_info=event_info)
end