
function fig_rbsp_overview_v01, plot_dir, event_info=event_info

    version = 'v01'
    id = '2015_0317'
    event_info = low_lshell_outflow_load_data(id)

    time_range = time_double(['2015-03-17','2015-03-18/12:00'])


    ;foreach probe, rbsp_probes do begin
    foreach probe, ['a','b'] do begin
        prefix = 'rbsp'+probe+'_'

        dens_var = rbsp_read_density(time_range, probe=probe, id='emfisis', suffix='')
        e_en_var = rbsp_read_en_spec(time_range, probe=probe, species='e', errmsg=errmsg)
        if errmsg ne '' then return, retval
        p_en_var = rbsp_read_en_spec(time_range, probe=probe, species='p', errmsg=errmsg)
        if errmsg ne '' then return, retval
        var_info = rbsp_read_en_spec_combo(time_range, probe=probe, species='o', errmsg=errmsg)
        if errmsg ne '' then return, retval
        o_en_vars = [var_info.para, var_info.anti]
        o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', errmsg=errmsg, energy_range=[50,5e4], update=1)

        
        ; E field spec.
        spec_var = rbsp_read_efield_spec(time_range, probe=probe, update=update)
        ; Cyclotron freqs.
        fc_vars = list()
        foreach species, ['e','o','he','p'] do fc_vars.add, rbsp_read_gyro_freq(time_range, probe=probe, species=species)
        var = prefix+'fce_half'
        fce = get_var_data(prefix+'fce', times=times)
        store_data, var, times, fce*0.5
        fc_vars.add, var
        var = prefix+'flh'
        fcp = get_var_data(prefix+'fcp', times=times)
        store_data, var, times, fcp*43
        fc_vars.add, var
        fc_vars = fc_vars.toarray()
        fc_colors = get_color(n_elements(fc_vars))
        foreach var, fc_vars, ii do options, var, 'colors', fc_colors[ii]

        spec_combo_var = spec_var+'_combo'
        store_data, spec_combo_var, data=[spec_var,fc_vars]
        options, spec_combo_var, 'yrange', get_setting(spec_var,'yrange')
        options, spec_combo_var, 'labels', ''


        ; Model field and tilt angle.
        external_model = 't89'
        internal_model = 'dipole'
        igrf = (internal_model eq 'igrf')? 1: 0
        suffix = '_'+internal_model+'_'+external_model
        bmod_gsm_var = prefix+'bmod_gsm'+suffix
        b_gsm_var = prefix+'b_gsm'
        foreach var, [b_gsm_var,bmod_gsm_var] do begin
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
        endforeach
        dtheta_var = prefix+'b_dtheta'
        b_theta = get_var_data(b_gsm_var+'_theta', times=times)
        bmod_theta = get_var_data(bmod_gsm_var+'_theta', at=times)
        dtheta = b_theta-bmod_theta
        store_data, dtheta_var, times, dtheta
        add_setting, dtheta_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'dTilt', $
            'unit', 'deg' )
        
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

        ; Tilt.
        var = bmod_gsm_var+'_theta'
        options, var, 'yrange', [0,90]
        options, var, 'ytickv', [40,80]
        options, var, 'constant', [40,80]
        options, var, 'yticks', 1
        options, var, 'yminor', 4
        var = dtheta_var
        yrange = [-45,10]
        ystep = 20
        ytickv = make_bins(yrange,ystep, inner=1)
        yticks = n_elements(ytickv)-1
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yrange', yrange
        options, var, 'yminor', 2
        options, var, 'constant', ytickv

        ; E
        e_var = prefix+'edot0_fac'
        var = e_var
        options, var, 'yrange', [-1,1]*20
        options, var, 'ytickv', [-1,0,1]*10
        options, var, 'yticks', 2
        options, var, 'yminor', 4
        options, var, 'constant', [-1,0,1]*10
        
        
        vars = [e_var,dens_var,e_en_var,o_en_var,o_pa_var,spec_combo_var,bmod_gsm_var+'_theta',dtheta_var,fmlat_var,mlat_var,dis_var]
        ypans = [1,1,1,1,1,1,0.6,0.8,0.6,0.6,0.6]
        if n_elements(xpansize) eq 0 then xpansize = 8
        pansize = [xpansize,1]
        nvar = n_elements(vars)
        margins = [10,8,8,1]
        
        if n_elements(plot_dir) eq 0 then plot_dir = event_info.plot_dir
        plot_file = join_path([plot_dir,'fig_rbsp_overview_'+prefix+version+'.pdf'])
        if keyword_set(test) then plot_file = 0
        poss = panel_pos(plot_file, ypans=ypans, fig_size=fig_size, nypan=nvar, pansize=pansize, panid=[0,1], margins=margins)
        sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, inch=1

        xticklen_chsz = -0.25   ; in ychsz.
        yticklen_chsz = -0.40   ; in xchsz.
        for ii=0,nvar-1 do begin
            tpos = poss[*,ii]
            xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
            yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
            options, vars[ii], 'xticklen', xticklen
            options, vars[ii], 'yticklen', yticklen
        endfor
        
        label_vars = prefix+['fmlat','mlt','dis']
        tplot, vars, trange=time_range, position=poss, var_label=label_vars
        if keyword_set(test) then stop
        sgclose
    endforeach
    

end

print, fig_rbsp_overview_v01()
end