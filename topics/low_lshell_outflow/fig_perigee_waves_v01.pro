
function fig_perigee_waves_v01, plot_file, event_info=event_info

    time_range = time_double(['2015-03-17','2015-03-18/12:00'])
    rbsp_probes = ['a','b']
    themis_probes = ['d','e']
    mms_probes = ['1']

    event_info = dictionary($
        'time_range', time_range, $
        'rbsp_probes', rbsp_probes, $
        'themis_probes', themis_probes, $
        ;'pp_ct', 67, $
        'pp_ct', 63, $
        'pp_psym', 8, $
        'pp_symsize', 0.75 )

    internal_model = 'igrf'
    external_model = 't89'
    igrf = internal_model eq 'igrf'
    t89_par = 2
    
    foreach probe, ['b'] do begin
        prefix = 'rbsp'+probe+'_'

        ; Density.
        dens_var = rbsp_read_density(time_range, probe=probe, id='emfisis', suffix='')
        ;p_info = rbsp_read_hope_moments(time_range, probe=probe, species='p')
        ;o_info = rbsp_read_hope_moments(time_range, probe=probe, species='o')
        p_dens_var = rbsp_read_density_hope(time_range, probe=probe, species='p')
        o_dens_var = rbsp_read_density_hope(time_range, probe=probe, species='o')
        

        ; Position.
        r_gsm_var = rbsp_read_orbit(time_range, probe=probe, coord='gsm')
        mlt_var = rbsp_read_mlt(time_range, probe=probe)
        lshell_var = rbsp_read_lshell(time_range, probe=probe)
        mlat_var = rbsp_read_mlat(time_range, probe=probe)
        options, mlat_var, 'yrange', [-1,1]*20
        options, mlat_var, 'constant', [-1,0,1]*10
        options, mlat_var, 'ytickv', [-1,0,1]*10
        options, mlat_var, 'yticks', 2
        options, mlat_var, 'yminor', 2
        
        ; E and B field.
        e_mgse_var = rbsp_read_efield(time_range, probe=probe, coord='mgse')
        b_gsm_var = rbsp_read_bfield(time_range, probe=probe, coord='gsm', errmsg=errmsg, resolution='hires')

        ; B model.
        bmod_gsm_var = prefix+'bmod_gsm_'+external_model+'_'+internal_model
        if check_if_update(bmod_gsm_var, time_range) then begin
            bmod_gsm_var = geopack_read_bfield(r_var=r_gsm_var, models=external_model, $
                suffix='_'+internal_model, igrf=igrf, t89_par=t89_par, coord='gsm')
            options, bmod_gsm_var, 'requested_time_range', time_range
        endif
        
        ; B-B_model.
        b_gsm = get_var_data(b_gsm_var, times=times)
        bmod_gsm = get_var_data(bmod_gsm_var, at=times)
        b_mag = snorm(b_gsm)
        bmod_mag = snorm(bmod_gsm)
        db_mag = b_mag-bmod_mag
        db_mag_var = prefix+'db_mag'
        store_data, db_mag_var, times, db_mag
        add_setting, db_mag_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'yrange', [-1,1]*200, $
            'unit', 'nT', $
            'short_name', 'dB_mag' )
        b_sm = cotran_pro(b_gsm, times, 'gsm2sm')
        bmod_sm = cotran_pro(bmod_gsm, times, 'gsm2sm')
        b_tilt = asin(b_sm[*,2]/b_mag)*constant('deg')
        bmod_tilt = asin(bmod_sm[*,2]/bmod_mag)*constant('deg')
        db_tilt = b_tilt-bmod_tilt
        db_tilt_var = prefix+'db_tilt'
        store_data, db_tilt_var, times, db_tilt
        add_setting, db_tilt_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'yrange', [-50,10], $
            'unit', 'deg', $
            'short_name', 'dB_tilt' )
        b_azim = atan(b_sm[*,1],b_sm[*,0])*constant('deg')
        bmod_azim = atan(bmod_sm[*,1],bmod_sm[*,0])*constant('deg')
        db_azim = b_azim-bmod_azim
        db_azim_var = prefix+'db_azim'
        store_data, db_azim_var, times, db_azim
        add_setting, db_azim_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'yrange', [-1,1]*40, $
            'unit', 'deg', $
            'short_name', 'dB_azim' )
        stop

        ; Trace to ionosphere.
        bf_gsm_var = prefix+'bf_gsm_t89_dipole_north'
        if check_if_update(bf_gsm_var, time_range) then begin
            vinfo_north = geopack_trace_to_ionosphere(r_gsm_var, models='t89', $
                igrf=0, north=1, refine=1, suffix='_dipole_north')
            options, bf_gsm_var, 'requested_time_range', time_range
        endif
        bf_gsm_var = prefix+'bf_gsm_t89_dipole_south'
        if check_if_update(bf_gsm_var, time_range) then begin
            vinfo_south = geopack_trace_to_ionosphere(r_gsm_var, models='t89', $
                igrf=0, south=1, refine=1, suffix='_dipole_south')
            options, bf_gsm_var, 'requested_time_range', time_range
        endif

        ; O spec.
        o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o')
        
        ; E/B ratio.
        
        stop
        
        ; Give me a time_range, a window, a suffix.
        if probe eq 'a' then begin
            tr = time_double(['2015-03-17/22:45:50','2015-03-17/22:46:10'])
            tr = time_double(['2015-03-17/22:45','2015-03-17/22:47'])
            window = 3d
            ;tr = time_double(['2015-03-18/07:57','2015-03-18/08:05'])
            tr = time_double(['2015-03-17/22:35','2015-03-17/22:40'])   ; Alfven waves?
            window = 300d
        endif else begin
            tr = time_double(['2015-03-17/21:30','2015-03-17/21:40'])   ; Alfven waves?
            window = 300d
        endelse

        suffix = '_1'
        tr1 = tr+[-1,1]*window
        
        ; dB.
        e_mgse = get_var_data(e_mgse_var, in=tr1, times=common_times, limits=e_lim)
        time_step = sdatarate(common_times)
        b_gsm = get_var_data(b_gsm_var, at=common_times, limits=b_lim)
        bmod_gsm = get_var_data(bmod_gsm_var, at=common_times)
        db_gsm_var = prefix+'db_gsm'+suffix
        db_gsm = b_gsm-bmod_gsm
        width = window/time_step
        ndim = 3
        for ii=0,ndim-1 do db_gsm[*,ii] -= smooth(db_gsm[*,ii],width, nan=1)
        store_data, db_gsm_var, common_times, db_gsm, limits=b_lim
        b0_gsm_var = prefix+'b0_gsm'+suffix
        b0_gsm = b_gsm-db_gsm
        store_data, b0_gsm_var, common_times, b0_gsm, limits=b_lim
        

        ; E dot B.
        edot0_mgse_var = prefix+'edot0_mgse'+suffix
        edot0_mgse = e_mgse
        b0_mgse = cotran(b0_gsm, common_times, 'gsm2mgse', probe=probe)
        edot0_mgse[*,0] = -total(e_mgse[*,1:2]*b0_mgse[*,1:2],2)/b0_mgse[*,0]
        for ii=0,ndim-1 do edot0_mgse[*,ii] -= smooth(edot0_mgse[*,ii],width, nan=1)
        e_lim.coord = 'mgse'
        store_data, edot0_mgse_var, common_times, edot0_mgse, limits=e_lim
        edot0_gsm_var = prefix+'edot0_gsm'+suffix
        edot0_gsm = cotran(edot0_mgse, common_times, 'mgse2gsm', probe=probe)
        e_lim.coord = 'gsm'
        store_data, edot0_gsm_var, common_times, edot0_gsm, limits=e_lim
        
        edot0_angle_var = prefix+'edot0_angle'+suffix
        edot0_angle = asin(b0_mgse[*,0]/snorm(b0_mgse))
        store_data, edot0_angle_var, common_times, edot0_angle*constant('deg')
        add_setting, edot0_angle_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'unit', 'deg', $
            'short_name', tex2str('theta')+'!DEdot0!N')
        
        
        ; FAC.
        fac_labels = ['||',tex2str('perp')+','+['west','out']]
        define_fac, b0_gsm_var, r_gsm_var, time_var=b0_gsm_var
        de_fac_var = prefix+'de_fac'+suffix
        to_fac, edot0_gsm_var, to=de_fac_var
        db_fac_var = prefix+'db_fac'+suffix
        to_fac, db_gsm_var, to=db_fac_var
        
        
        ; Spec.
        de_fac_comp_var = stplot_index(de_fac_var, 1)
        db_fac_comp_var = stplot_index(db_fac_var, 2)
        de_mor_var = stplot_mor_new(de_fac_comp_var)
        db_mor_var = stplot_mor_new(db_fac_comp_var)
        zlim, de_mor_var, [1e-1,1e3]
        zlim, db_mor_var, [1e-2,1e2]


        ; Pflux.
        pf_fac_var = prefix+'pf_fac'+suffix
;        cpoynt = 1d/(400d*!dpi) ; from mV/m x nT -> mW/m^2.
;        de_fac = get_var_data(de_fac_var)
;        db_fac = get_var_data(db_fac_var)
;        pf_fac = cpoynt*vec_cross(de_fac, db_fac)
;        store_data, pf_fac_var, common_times, pf_fac
        stplot_calc_pflux_mor, de_fac_var, db_fac_var, pf_fac_var, scaleinfo=scale_info
        add_setting, pf_fac_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'S', $
            'unit', 'mW/m!U2!N', $
            'coord', 'FAC', $
            'coord_labels', fac_labels )

        ; Map.
        bf_gsm = get_var_data(bf_gsm_var, at=common_times)
        cmap_var = prefix+'cmap'+suffix
        cmap = snorm(bf_gsm)/snorm(b0_gsm)
        store_data, cmap_var, common_times, cmap

        pf_fac_map_var = prefix+'pf_fac_map'+suffix
        pf_fac_map = get_var_data(pf_fac_var)
        for ii=0,ndim-1 do pf_fac_map[*,ii] *= cmap
        store_data, pf_fac_map_var, common_times, pf_fac_map
        add_setting, pf_fac_map_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'S', $
            'unit', 'mW/m!U2!N', $
            'coord', 'FAC', $
            'coord_labels', fac_labels )

        
        ; O+ and H+ ratio.
        o_dens = get_var_data(prefix+'o_density', at=common_times)
        p_dens = get_var_data(prefix+'p_density', at=common_times)
        o_ratio = o_dens/(o_dens+p_dens)
        p_ratio = 1-o_ratio
        o_mass = 16d
        p_mass = 1d
        avg_mass = o_mass*o_ratio+p_mass*p_ratio
        
        ; Alfven speed.
        va_var = prefix+'va'+suffix
        va0 = 1e-9/sqrt(1e6*!dpi*4e-7*1.67e-27)*1e-3    ; km/s, B in nT, n in cc, m in m_p.
        bmag = snorm(b0_gsm)
        num_dens = get_var_data(dens_var, at=common_times)
        
        va = va0*bmag/sqrt(num_dens*avg_mass)
        store_data, va_var, common_times, va
        add_setting, va_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'V!DA!N', $
            'unit', 'km/s' )
        
        ; Gyro-frequency.
        f_g0 = 1.6e-19*1e-9/1.67e-27/2/!dpi   ; in Hz.
        fci_var = prefix+'fci'+suffix
        fci = f_g0*bmag/avg_mass
        store_data, fci_var, common_times, fci
        add_setting, fci_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'f!Dci!N', $
            'unit', 'Hz' )
        fco_var = prefix+'fco'+suffix
        fco = f_g0*bmag/o_mass
        store_data, fco_var, common_times, fco
        add_setting, fco_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'f!Dco!N', $
            'unit', 'Hz' )
        fcp_var = prefix+'fcp'+suffix
        fcp = f_g0*bmag/p_mass
        store_data, fcp_var, common_times, fcp
        add_setting, fcp_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'f!Dcp!N', $
            'unit', 'Hz' )
            
        plot_vars = [pf_fac_map_var,de_fac_var,db_fac_var,de_mor_var,db_mor_var,$
            o_pa_var,va_var,lshell_var]
        nplot_var = n_elements(plot_vars)
        poss = sgcalcpos(nplot_var)
        tplot, plot_vars, trange=tr, position=poss
        
        foreach var, [db_mor_var,de_mor_var] do begin
            pid = where(plot_vars eq var)
            tpos = poss[*,pid]
            yrange = get_setting(var, 'yrange')
            xrange = tr
            ylog = get_setting(var, 'ylog')
            xlog = 0
            plot, xrange, yrange, $
                xstyle=5, xrange=xrange, xlog=xlog, $
                ystyle=5, yrange=yrange, ylog=ylog, $
                position=tpos, nodata=1, noerase=1
            f_vars = [fcp_var,fco_var,fci_var]
            nf_var = n_elements(f_vars)
            colors = get_color(nf_var)
            foreach f_var, f_vars, var_id do begin
                fc = get_var_data(f_var, in=xrange, times=uts)
                oplot, uts, fc, color=colors[var_id]
            endforeach
        endforeach

        

        stop
        e_en_var = rbsp_read_en_spec(time_range, probe=probe, species='e', errmsg=errmsg)
        if errmsg ne '' then return, retval
        p_en_var = rbsp_read_en_spec(time_range, probe=probe, species='p', errmsg=errmsg)
        if errmsg ne '' then return, retval
        var_info = rbsp_read_en_spec_combo(time_range, probe=probe, species='o', errmsg=errmsg)
        if errmsg ne '' then return, retval
        o_en_vars = [var_info.para, var_info.anti]
        
        stop
        
        fac_labels = ['b','w','o']

        ; Model field.
        external_model = 't89'
        internal_model = 'dipole'
        igrf = (internal_model eq 'igrf')? 1: 0
        t89_par = 2
        bmod_gsm_var = geopack_read_bfield(time_range, r_var=r_gsm_var, models=external_model, igrf=igrf, suffix='_'+internal_model, t89_par=t89_par)
        vinfo_north = geopack_trace_to_ionosphere(r_gsm_var, models='t89', igrf=0, north=1, refine=1, suffix='_'+internal_model+'_north')
        vinfo_south = geopack_trace_to_ionosphere(r_gsm_var, models='t89', igrf=0, south=1, refine=1, suffix='_'+internal_model+'_south')

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
        fmlat_north_var = vinfo_north.fmlat
        fmlat_south_var = vinfo_south.fmlat
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
        var = e_var
        options, var, 'yrange', [-1,1]*20
        options, var, 'ytickv', [-1,0,1]*10
        options, var, 'yticks', 2
        options, var, 'yminor', 4
        options, var, 'constant', [-1,0,1]*10
        
        vars = [e_var,dens_var,e_en_var,p_en_var,o_en_vars,bmod_gsm_var+'_theta',dtheta_var,fmlat_var,mlat_var,dis_var]
        ypans = [1,1,1,1,1,1,0.6,0.8,0.6,0.6,0.6]
        if n_elements(xpansize) eq 0 then xpansize = 8
        pansize = [xpansize,1]
        nvar = n_elements(vars)
        margins = [10,4,8,1]
        
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
        
        tplot, vars, trange=time_range, position=poss
        stop
    endforeach
    

end

print, fig_perigee_waves_v01()
end