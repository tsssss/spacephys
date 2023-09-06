pro _2014_0828_10_load_weygand, the_info, reload=reload, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload) then begin
        file_delete, the_info['file'], /allow_nonexistent
        del_data, the_info['var']
    endif

    ; Settings.
    time_range = the_info.time_range
    mlon_range = the_info.mlon_range
    mlat_range = the_info.mlat_range
    mlt_range = the_info.mlt_range
    ewo_mlat_range = the_info.ewo_mlat_range
    mlon_binsize = the_info.mlon_binsize
    mlat_binsize = the_info.mlat_binsize
    mlt_binsize = mlon_binsize/15.
    j_var = 'thg_j_ver'
    mlt_ewo_var = 'thg_j_up_ewo'


    ; Update.
    update_data_file = 0
    load = 0
    save_vars = the_info.var
    foreach tvar, save_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        if file_test(the_info.file) eq 1 then begin
            tplot_restore, filename=the_info.file
        endif else begin
            mlat_bins = make_bins(mlat_range,mlat_binsize)
            mlon_bins = make_bins(mlon_range,mlon_binsize)
            nmlon_bin = n_elements(mlon_bins)
            nmlat_bin = n_elements(mlat_bins)
            new_image_size = [nmlon_bin,nmlat_bin]

            if check_if_update(j_var, time_range) then themis_read_weygand, time_range
            get_data, j_var, times, j_orig
            glonbins = get_setting(j_var, 'glonbins')
            glatbins = get_setting(j_var, 'glatbins')
            glonbinsize = total(glonbins[[0,1]]*[-1,1])
            glatbinsize = total(glatbins[[0,1]]*[-1,1])
            nglonbin = n_elements(glonbins)
            nglatbin = n_elements(glatbins)
            old_image_size = [nglonbin,nglatbin]

            ; Get the mlon/mlat for glon/glat bins.
            pixel_glons = fltarr(nglonbin,nglatbin)
            pixel_glats = fltarr(nglonbin,nglatbin)
            for ii=0,nglatbin-1 do pixel_glons[*,ii] = glonbins
            for ii=0,nglonbin-1 do pixel_glats[ii,*] = glatbins
            apexfile = join_path([homedir(),'Projects','idl','spacephys','aurora','image','support','mlatlon.1997a.xdr'])
            geo2apex, pixel_glats, pixel_glons, pixel_mlats, pixel_mlons

            ; Map to uniform mlon/mlat bins.
            mlon_bin_min = mlon_range[0]
            mlat_bin_min = mlat_range[0]
            i0_bins = round((pixel_mlons-mlon_bin_min)/mlon_binsize)
            j0_bins = round((pixel_mlats-mlat_bin_min)/mlat_binsize)

            i1_range = [0,nmlon_bin-1]
            j1_range = [0,nmlat_bin-1]

            i_bins = make_bins(i1_range, 1)
            j_bins = make_bins(j1_range, 1)
            ni_bin = nmlon_bin
            nj_bin = nmlat_bin

            index_map_from_old = list()
            index_map_to_new = list()
            for ii=0, nmlon_bin-1 do begin
                mlon_range = mlon_bins[ii]+[-1,1]*mlon_binsize*0.5
                for jj=0, nmlat_bin-1 do begin
                    mlat_range = mlat_bins[jj]+[-1,1]*mlat_binsize*0.5
                    index = where($
                        pixel_mlons ge mlon_range[0] and $
                        pixel_mlons lt mlon_range[1] and $
                        pixel_mlats ge mlat_range[0] and $
                        pixel_mlats lt mlat_range[1], count)
                    if count eq 0 then continue
                    index_map_from_old.add, index
                    index_map_to_new.add, ii+jj*nmlon_bin
                endfor
            endfor

            ntime = n_elements(times)
            j_new = fltarr([ntime,new_image_size])
            for ii=0,ntime-1 do begin
                img_old = reform(j_orig[ii,*,*])
                img_new = fltarr(new_image_size)
                foreach pixel_new, index_map_to_new, pixel_id do begin
                    img_new[pixel_new] = mean(img_old[index_map_from_old[pixel_id]])
                endforeach
                j_new[ii,*,*] = img_new
            endfor

        ;---EWOgram.
            mlat_range = [60,70]
	    stop
	    ; Need to change to the packaged routine.
            mlat_index = where_pro(mlat_bins, '[]', mlat_range)
            ewo = total(-j_new[*,*,mlat_index], 3)/n_elements(mlat_index)
            ewo = fltarr(ntime,nmlon_bin)
            foreach time, times, ii do begin
                foreach mlon, mlon_bins, jj do begin
                    tmp = reform(-j_new[ii,jj,mlat_index])
                    index = where(tmp lt 0, count)
                    if count ne 0 then tmp[index] = 0
                    ewo[ii,jj] = mean(tmp)   ; total -> zrange [0.2.5e5], max -> zrange [0,2.5e5], mean -> zrange [0,1.5e5], total -> zrange [3e5]
                endforeach
            endforeach
            store_data, ewo_var, times, ewo, mlon_bins, limits={$
                spec: 1, $
                no_interp: 1, $
                ytitle: 'MLon (deg)', $
                ystyle: 1, $
                yrange: reverse(minmax(mlon_bins)), $
                ztitle: 'Upward current (A)', $
                zlog: 0 , $
                zrange: [0,.5e5], $
                yticklen: -0.02, $
                xticklen: -0.02 }

        ;---Convert EWOgram from MLon to MLT.
            mlt_bins = make_bins(mlt_range, mlt_binsize)
            nmlt_bin = n_elements(mlt_bins)
            get_data, ewo_var, times, ewo, mlon_bins
            ntime = n_elements(times)
            mlt_ewo = fltarr(ntime,nmlt_bin)
            for ii=0,ntime-1 do begin
                the_mlts = mlon2mlt(mlon_bins,times[ii])
                dmlt = the_mlts[1:-1]-the_mlts[0:-2]
                index = where(abs(dmlt) gt 12, count)
                if count ne 0 then begin
                    if dmlt[index] ge 0 then begin
                        the_mlts[index+1:*] -= 24
                    endif else begin
                        the_mlts[index+1:*] += 24
                    endelse
                endif

                mlt_ewo[ii,*] = interpol(ewo[ii,*],the_mlts,mlt_bins)
                index = where(mlt_bins le min(the_mlts) or mlt_bins ge max(the_mlts), count)
                if count ne 0 then mlt_ewo[ii,index] = 0
            endfor
            ystep = 3
            ytickv = make_bins(mlt_range, ystep)
            yticks = n_elements(ytickv)-1
            yminor = ystep
            store_data, mlt_ewo_var, times, mlt_ewo, mlt_bins, limits={$
                spec: 1, $
                no_interp: 1, $
                ytitle: 'MLT (hr)', $
                ystyle: 1, $
                yrange: mlt_range, $
                ytickv: ytickv, $
                yticks: yticks, $
                yminor: yminor, $
                ztitle: 'Upward current (A)', $
                zlog: 0 , $
                zrange: [0,.5e5], $
                yticklen: -0.02, $
                xticklen: -0.02 }


            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['weygand'] = the_info
        tplot_save, save_vars, filename=the_info.file

        info_var = event_info.var
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info.file
    endif



end


pro _2014_0828_10_calc_pflux, the_info, reload=reload, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload) then begin
        file_delete, the_info['file'], /allow_nonexistent
        del_data, the_info['var']
    endif

    ; Settings.
    scale_info = the_info.scale_info
    probes = the_info.probes
    time_range = the_info.psd_time_range
    fac = ['||','east','north']
    rgb = constant('rgb')
    pf_unit = 'mW/m!U2!N'
    ct = 66
    nprobe = n_elements(probes)

    ; Update.
    update_data_file = 0
    load = 0
    save_vars = the_info.var
    foreach tvar, save_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        if file_test(the_info.file) eq 1 then begin
            tplot_restore, filename=the_info.file
        endif else begin
            foreach probe, probes do begin
                prefix = probe+'_'

                ; Uniform time.
                vars = prefix+['b0_gsm','db_gsm','r_gsm']
                foreach var, vars do interp_time, var, to=prefix+'e_gsm'

                ; Convert E/B fields to FAC.
                define_fac, prefix+'b0_gsm', prefix+'r_gsm'
                vars = prefix+['db_gsm','e_gsm']
                foreach var, vars do to_fac, var

                vars = prefix+['db','e']+'_fac'
                short_names = ['dB','dE']
                foreach var, vars, ii do begin
                    add_setting, var, /smart, {$
                        display_type: 'vector', $
                        short_name: short_names[ii], $
                        coord: 'FAC', $
                        coord_labels: fac, $
                        colors: rgb}
                endforeach

                ; Calculate pflux.
                stplot_calc_pflux_mor, prefix+'e_fac', prefix+'db_fac', prefix+'pf_fac', scaleinfo=scale_info

                scales = wavelet_make_scale(s0=scale_info.s0, sj=scale_info.s1, dj=scale_info.dj)
                vars = prefix+['db','e']
                foreach var, vars do begin
                    get_data, var+'_gsm', times, data
                    index = where_pro(times, time_range)
                    store_data, var+'_mag', times[index], snorm(data[index,*])
                    calc_psd, var+'_mag', scales=scales
                endforeach

                tvar = prefix+'pf_fac_mor_spec_1'
                get_data, tvar, uts, dat
                index = where_pro(uts, time_range, count=nrec)
                spsd = total(dat[index,*], 1)/nrec

                get_data, prefix+'db_mag_psd_cwt', freqs, bpsd
                get_data, prefix+'e_mag_psd_cwt', freqs, epsd
                eb_ratio = sqrt(epsd/bpsd)*1e3  ; in km/s.
                add_setting, prefix+'pf_fac', /smart, {$
                    display_type: 'vector', $
                    unit: 'mW/m!U2!N', $
                    short_name: 'S', $
                    coord: 'FAC', $
                    coord_labels: fac, $
                    colors: rgb, $
                    freq: freqs, $
                    freq_unit: 'Hz', $
                    constant: 0., $
                    e_psd: epsd, $
                    e_psd_unit: '(mV/m)!U2!N/Hz', $
                    b_psd: bpsd, $
                    b_psd_unit: '(nT)!U2!N/Hz', $
                    s_psd: spsd, $
                    eb_ratio: eb_ratio, $
                    eb_ratio_unit: 'km/s'}
            endforeach
            update_data_file = 1
        endelse
    endif

    ; Update file and others.
    if update_data_file then begin
        event_info['pflux'] = the_info
        tplot_save, save_vars, filename=the_info.file

        info_var = event_info.var
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info.file
    endif

end


pro _2014_0828_10_calc_drift_period, drift_period_info, reload=reload_drift_period, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_drift_period) then begin
        file_delete, drift_period_info['file'], /allow_nonexistent
        del_data, drift_period_info['var']
    endif


    all_species = drift_period_info.species
    test_times = drift_period_info.test_times
    save_vars = drift_period_info.var

    update_data_file = 0
    load = 0
    foreach tvar, save_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        if file_test(drift_period_info['file']) eq 1 then tplot_restore, filename=drift_period_info['file'] else begin

            ; Loop over species.
            foreach species, all_species do begin
                lprmsg, 'Processing '+species+' ...'
                the_info = drift_period_info[species]
                probes = the_info.probes
                energy_range = the_info.energy_range
                flux_vars = probes+'_kev_'+species+'_flux'
                drift_vars = probes+'_'+species+'_drift_period'
                foreach probe, probes, probe_id do begin
                    lprmsg, 'Processing '+probe+' ...'
                    prefix = probe+'_'
                    drift_period_var = drift_vars[probe_id]
;if tnames(drift_period_var) ne '' then continue
                    get_data, flux_vars[probe_id], times, fluxes, energies
                    ; Filter energy.
                    index = where_pro(energies,'[]', energy_range)
                    energies = energies[index]
                    ; Select times to calc drift period.
                    times = test_times
                    ntime = n_elements(times)
                    nenergy = n_elements(energies)
                    drift_periods = fltarr(ntime,nenergy)
                    bounce_periods = fltarr(ntime,nenergy)
                    lshells = fltarr(ntime,nenergy)
                    ; Calc drift period.
                    pos_var = prefix+'r_gsm'
                    r_gsms = get_var_data(pos_var, at=times)
                    foreach time, times, ii do begin
                        foreach energy, energies, jj do begin
                            drift_periods[ii,jj] = calc_drift_period(rgsm=r_gsms[ii,*], species=species, $
                                time=time, energy=energy, bounce_period=bounce_period, lshell=lshell)
                            bounce_periods[ii,jj] = bounce_period
                            lshells[ii,jj] = lshell
                        endforeach
                    endforeach

                    store_data, drift_period_var, times, drift_periods, energies
                    add_setting, drift_period_var, /smart, {$
                        display_type: 'list', $
                        bounce_period: bounce_periods, $
                        lshells: lshells, $
                        color_table: 52, $
                        unit: 'sec', $
                        value_unit: 'keV', $
                        short_name: 'T drift'}
                endforeach
            endforeach
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['calc_drift_period'] = drift_period_info
        tplot_save, save_vars, filename=drift_period_info['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif





end

pro _2014_0828_10_calc_gmag_info, gmag_info, event_info=event_info

    gmag = event_info['gmag_info']

    ; Get the basic info.
    sites = gmag['sites']
    foreach site, sites do begin
        gmag_info[site] = dictionary(themis_read_mag_metadata(sites=site))
    endforeach

    ; Calculate the time FAC arrives.
    mag_var = gmag['var']
    get_data, mag_var, times, mag_data

    data_rate = sdatarate(times)
    base_time = gmag_info['base_time']
    slope_width = gmag_info['slope_width']
    time_range = gmag_info['time_range']
    foreach site, sites, ii do begin
        f0 = mag_data[*,ii]
        ; Find the time when value start to decrease.
        ; By using the time when its derivative excedes 3 sigma of the background level.
        dn0 = slope_width/data_rate
        uts = make_bins(time_range, slope_width)
        nut = n_elements(uts)
        f1 = interpol(f0,times, uts)
        slopes = (f1[1:nut-1]-f1[0:nut-2])/slope_width
        uts = (uts[1:nut-1]+uts[0:nut-2])*0.5
        ; Determine the backgrond level of the slope.
        index = where(uts le base_time)
        slope_background_level = stddev(slopes[index])
        slope_threshold = -slope_background_level*3
        if min(slopes) ge slope_threshold then slope_threshold = min(slopes)*0.5
        ; Find the start time of the search.
        index = where(slopes le slope_threshold and uts ge base_time, count)
        if count eq 0 then stop     ; no site stops.
        start_time = uts[index[0]]
        (gmag_info[site])['start_time'] = start_time[0]
        ; Find the end time of the search.
        index = where(uts ge start_time)
        f1_min = min(f1[index])
        index = where(uts ge start_time and f1 eq f1_min)
        end_time = uts[index]
        (gmag_info[site])['end_time'] = end_time[0]
        ; Find the time when FAC arrives.
        index = where_pro(times, [start_time,end_time])
        f1 = f0[index]
        t1 = times[index]
        value_threshold = (gmag_info['fac_method'] eq 'absolute')? f1[0]-gmag_info['fac_threshold']: f1[0]+(f1_min-f1[0])*gmag_info['fac_threshold']
        index = where_pro(f1, 'le', value_threshold, count=count)
        fac_time1 = t1[index[0]]   ; should always find a result.
        index = where_pro(f1, 'ge', value_threshold, count=count)
        fac_time2 = t1[index[count-1]]   ; should always find a result.
        fac_time = (abs(f1[0]-f1_min) le gmag_info['fac_db_threshold'])? fac_time1: fac_time2
        (gmag_info[site])['fac_time'] = fac_time[0]

        test = 0
        if keyword_set(test) then begin
            test_var = 'test_var'
            store_data, test_var, times, f0, limits={labels:site}
            store_data, test_var+'_df', uts, slopes, limits={labels:site}
            tplot, test_var+['','_df']
            timebar, start_time, color=sgcolor('red')
            timebar, end_time, color=sgcolor('red')
            timebar, fac_time, color=sgcolor('blue')
            print, time_string(start_time)
            print, time_string(end_time)
            print, time_string(fac_time)
            stop
        endif
    endforeach

    event_info['gmag_info'] = gmag_info
    info_var = event_info['var']
    store_data, info_var, 0, event_info
    tplot_save, info_var, filename=event_info['file']
end


pro _2014_0828_10_load_event_info, event_info=event_info
    ; Haven't decide if this is useful.
end

pro _2014_0828_10_load_gmag, gmag, reload=reload_gmag, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_gmag) then begin
        file_delete, gmag['file'], /allow_nonexistent
        store_data, gmag['var'], /delete
    endif

    update_data_file = 0
    load = 0
    foreach tvar, gmag['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(gmag['file']) eq 1 then tplot_restore, filename=gmag['file'] else begin
            themis_read_mag, gmag['time_range'], /sort_by_mlon, $
                mlat_range=gmag['mlat_range'], mlon_range=gmag['mlon_range'];, component=comp

            ; Add sites to gmag.
            gmag_var = gmag['var']
            get_data, gmag_var, tmp, tmp, sites
            gmag['sites'] = strlowcase(sites)
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['gmag'] = gmag
        tplot_save, gmag['var'], filename=gmag['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif
end


pro _2014_0828_10_load_aurora, aurora, reload=reload_aurora, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_aurora) then begin
        file_delete, aurora['file'], /allow_nonexistent
        store_data, aurora['var'], /delete
    endif

    update_data_file = 0
    load = 0
    foreach tvar, aurora['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(aurora['file']) eq 1 then tplot_restore, filename=aurora['file'] else begin
            sites = aurora['sites']
            nsite = n_elements(sites)
            site_infos = replicate(aurora[sites[0]].tostruct(),nsite)
            for ii=0, nsite-1 do site_infos[ii] = aurora[sites[ii]].tostruct()
            themis_read_mlonimg, aurora['time_range'], sites=sites, site_infos=site_infos, mlon_range=aurora['mlon_range'], mlat_range=aurora['mlat_range'], merge_method=aurora['merge_method']
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['aurora'] = aurora
        tplot_save, aurora['var'], filename=aurora['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif
end


pro _2014_0828_10_load_efield, efield, reload=reload_efield, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_efield) then begin
        file_delete, efield['file'], /allow_nonexistent
        store_data, efield['var'], /delete
    endif

    probe_abbrev = event_info['probe_abbrev']
    update_data_file = 0

    ; Do nothing if data are found in memory.
    load = 0
    foreach tvar, efield['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(efield['file']) eq 1 then tplot_restore, filename=efield['file'] else begin
            all_probes = efield['probe']
            time_range = efield['time_range']
            foreach key, all_probes.keys() do begin
                probes = all_probes[key]
                prefix = probe_abbrev[key]
                foreach probe, probes do begin
                    ; Load xxx_exx_gsm.
                    routine = key+'_read_efield'
                    call_procedure, routine, time_range, probe=probe
                    pre0 = prefix+probe+'_'

                    case key of
                        'rbsp': begin
                            get_data, pre0+'e_gse', times, egse
                            egsm = cotran(egse, times, 'gse2gsm')
                            store_data, pre0+'e0_gsm', times, egsm
                            add_setting, pre0+'e0_gsm', /smart, {$
                                display_type: 'vector', $
                                unit: 'mV/m', $
                                short_name: 'E0', $
                                coord: 'GSM', $
                                coord_labels: ['x','y','z'], $
                                colors: constant('rgb')}
                            end
                        'themis': ; do nothing.
                    endcase

                    ; Rename variable.
                    evar = pre0+'e_gsm'
                    case key of
                        'rbsp': rename_var, pre0+'e0_gsm', to=evar
                        'themis': ; do nothing.
                        else: stop
                    endcase

                    options, evar, 'colors', event_info['rgb']
                endforeach
            endforeach
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['efield'] = efield
        tplot_save, efield['var'], filename=efield['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif
end


pro _2014_0828_10_load_bfield, bfield, reload=reload_bfield, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_bfield) then begin
        file_delete, bfield['file'], /allow_nonexistent
        store_data, bfield['var'], /delete
    endif

    probe_abbrev = event_info['probe_abbrev']
    update_data_file = 0

    ; Do nothing if data are found in memory.
    load = 0
    foreach tvar, bfield['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(bfield['file']) eq 1 then tplot_restore, filename=bfield['file'] else begin
            all_probes = bfield['probe']
            time_range = bfield['time_range']
            model_chosen = bfield['model_chosen']
            foreach key, all_probes.keys() do begin
                probes = all_probes[key]
                prefix = probe_abbrev[key]
                foreach probe, probes do begin
                    ; Load xxx_b_gsm.
                    routine = key+'_read_bfield'
                    call_procedure, routine, time_range, probe=probe
                    pre0 = prefix+probe+'_'

                    ; Separate dB and B0.
                    calc_db, pre0+'b_gsm', pre0+'bmod_gsm_'+model_chosen, db_name=pre0+'db_gsm', b0_name=pre0+'b0_gsm'

                    ; Correct the mapping coefficient.
                    get_data, pre0+'c0map_'+model_chosen, times, c0map
                    get_data, pre0+'b0_gsm', uts, b0gsm
                    b0gsm = sinterpol(b0gsm, uts, times)
                    b0mag = snorm(b0gsm)
                    bmodmag = snorm(get_var_data(pre0+'bmod_gsm_'+model_chosen))
                    cmap = c0map*bmodmag/b0mag
                    cmap_var = pre0+'cmap'
                    store_data, cmap_var, times, cmap
                    add_setting, cmap_var, /smart, {$
                        display_type:'scalar', $
                        unit:'#', $
                        short_name:'C!Dmap!N'}

                    ; Calc |B|-|B|model.
                    get_data, pre0+'b_gsm', times, bgsm
                    bmag = snorm(bgsm)
                    get_data, pre0+'bmod_gsm_'+model_chosen, uts, bmodgsm
                    bmodgsm = sinterpol(bmodgsm, uts, times)
                    bmodmag = snorm(bmodgsm)
                    dbmag = bmag-bmodmag
                    sdespike, times, dbmag
                    tvar = pre0+'dbmag'
                    store_data, tvar, times, dbmag
                    add_setting, tvar, /smart, {$
                        display_type: 'scalar', $
                        unit: 'nT', $
                        short_name: '|B|-|B!Dmod!N|'}

                    vars = pre0+['b_gsm','db_gsm','b0_gsm']
                    options, vars, 'colors', event_info['rgb']
                endforeach
            endforeach
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['bfield'] = bfield
        tplot_save, bfield['var'], filename=bfield['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif
end


pro _2014_0828_10_load_model_info, model_info, reload=reload_model_info, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_model_info) then begin
        file_delete, model_info['file'], /allow_nonexistent
        store_data, model_info['var'], /delete
    endif

    probe_abbrev = event_info['probe_abbrev']
    update_data_file = 0

    ; Do nothing if data are found in memory.
    load = 0
    foreach tvar, model_info['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(model_info['file']) eq 1 then tplot_restore, filename=model_info['file'] else begin
            trace_dir = model_info['trace_dir']
            foreach model, model_info['models'] do begin
                foreach key, model_info['probe'].keys() do begin
                    probes = (model_info['probe'])[key]
                    prefix = probe_abbrev[key]
                    foreach probe, probes do begin
                        rvar = prefix+probe+'_r_gsm'
                        if tnames(rvar) eq '' then message, 'Load R first ...'
                        lprmsg, prefix+probe+' '+model+' ...'
                        read_geopack_info, rvar, model=model, direction=trace_dir

                        pre0 = prefix+probe+'_'
                        vars = pre0+['bmod_gsm','fpt_gsm']+'_'+model
                        options, vars, 'colors', event_info['rgb']
                    endforeach
                endforeach
            endforeach
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['model_info'] = model_info
        tplot_save, model_info['var'], filename=model_info['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif
end


pro _2014_0828_10_load_injection, injection, reload=reload_injection, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_injection) then begin
        file_delete, injection['file'], /allow_nonexistent
        store_data, injection['var'], /delete
    endif

    probe_abbrev = event_info['probe_abbrev']
    update_data_file = 0

    ; Do nothing if data are found in memory.
    load = 0
    foreach tvar, injection['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(injection['file']) eq 1 then tplot_restore, filename=injection['file'] else begin
            all_probes = injection['probe']
            time_range = injection['time_range']
            foreach key, all_probes.keys() do begin
                probes = all_probes[key]
                prefix = probe_abbrev[key]
                foreach species, injection['species'].keys() do begin
                    vars = prefix+probes+'_'+'kev_'+species+'_flux'
                    routine = key+'_read_kev_'+(injection['species'])[species]
                    foreach probe, probes, ii do begin
                        call_procedure, routine, time_range, probe=probe, energy=injection['energy_range'], pitch_angle=injection['pitch_angle']
                    endforeach
                endforeach
            endforeach
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['injection'] = injection
        tplot_save, injection['var'], filename=injection['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif


end


pro _2014_0828_10_load_ephemeris, ephemeris, reload=reload_ephemeris, event_info=event_info

    ; Remove data from disk and memory, to force reloading data.
    if keyword_set(reload_ephemeris) then begin
        file_delete, ephemeris['file'], /allow_nonexistent
        store_data, ephemeris['var'], /delete
    endif

    probe_abbrev = event_info['probe_abbrev']
    update_data_file = 0

    ; Do nothing if data are found in memory.
    load = 0
    foreach tvar, ephemeris['var'] do if tnames(tvar) eq '' then load = 1
    if load then begin
        ; Try to load data from disk.
        if file_test(ephemeris['file']) eq 1 then tplot_restore, filename=ephemeris['file'] else begin
            all_probes = ephemeris['probe']
            time_range = ephemeris['time_range']
            foreach key, all_probes.keys() do begin
                probes = all_probes[key]
                prefix = probe_abbrev[key]
                foreach probe, probes do begin
                    ; Load xxx_r_gsm.
                    routine = key+'_read_orbit'
                    call_procedure, routine, time_range, probe=probe
                    if key eq 'rbsp' then begin
                        get_data, prefix+probe+'_r_gse', times, r_gse
                        store_data, prefix+probe+'_r_gsm', times, cotran(r_gse,times,'gse2gsm')
                        add_setting, prefix+probe+'_r_gsm', /smart, dictionary($
                            'display_type', 'vector', $
                            'unit', 'Re', $
                            'short_name', 'R', $
                            'coord', 'GSM', $
                            'coord_labels', ['x','y','z'] )
                    endif

                    ; Fix r_gsm for 1991-080.
                    if probe eq '1991-080' then begin
                        lanl_read_orbit, make_bins(time_range,86400d), probe=probe
                        bad_time = time_double(['2014-08-28/05:30','2014-08-28/16:30'])
                        var = probe+'_r_gsm'
                        fix_lanl_orbit, var, bad_time=bad_time, to=var
                        get_data, var, times, rgsm
                        index = where_pro(times, time_range)
                        store_data, var, times[index], rgsm[index,*]
                    endif

                    ; Get other ephemeris data.
                    deg = (event_info['constant'])['deg']
                    pre0 = prefix+probe+'_'

                    get_data, pre0+'r_gsm', times, rgsm
                    options, pre0+'r_gsm', 'colors', event_info['rgb']
                    rmag = cotran(rgsm, times, 'gsm2mag')
                    mlat = asin(rmag[*,2]/snorm(rmag))*deg
                    mlon = atan(rmag[*,1],rmag[*,0])*deg
                    mlt = mlon2mlt(mlon, times)

                    ; Save them to memory.
                    mlat_var = pre0+'mlat'
                    store_data, mlat_var, times, mlat
                    add_setting, mlat_var, /smart, {$
                        display_type: 'scalar', $
                        unit: 'deg', $
                        short_name: 'MLat'}

                    mlon_var = pre0+'mlon'
                    store_data, mlon_var, times, mlon
                    add_setting, mlon_var, /smart, {$
                        display_type: 'scalar', $
                        unit: 'deg', $
                        short_name: 'MLon'}

                    mlt_var = pre0+'mlt'
                    store_data, mlt_var, times, mlt
                    add_setting, mlt_var, /smart, {$
                        display_type: 'scalar', $
                        unit: 'deg', $
                        short_name: 'MLT'}
                endforeach
            endforeach
            update_data_file = 1
        endelse
    endif

    if update_data_file then begin
        event_info['ephemeris'] = ephemeris
        tplot_save, ephemeris['var'], filename=ephemeris['file']

        info_var = event_info['var']
        store_data, info_var, 0, event_info
        tplot_save, info_var, filename=event_info['file']
    endif
end


pro _2014_0828_10_load_data, $
    event_info=event_info, $                ; output.
    reload_event_info=reload_event_info, $  ; Load event info.
    reload_ephemeris=reload_ephemeris, $    ; Load position in GSM.
    reload_injection=reload_injection, $    ; Need to set probes, energy range, pitch angle.
    reload_drift_period=reload_drift_period, $  ; Need injection and position.
    reload_model_info=reload_model_info, $  ; Need position. Load footpoint, model B, mapping coefficient.
    reload_bfield=reload_bfield, $          ; Need model info. Load B in GSM, dB and B0.
    reload_efield=reload_efield, $          ; Need some settings.
    reload_pflux=reload_pflux, $            ; Need E and B loaded.
    reload_aurora=reload_aurora, $          ; Need some settings.
    reload_gmag=reload_gmag, $              ; Need some settings.
    recalc_gmag_info=recalc_gmag_info, $    ; Need mag data.
    reload_weygand=reload_weygand, $        ; Need some settings.
    reload_all=reload_all


;---Settings.
    event_info = dictionary()

    event_id = '2014_0828_10'
    event_info['id'] = event_id
    info_var = event_id+'_event_info'
    event_info['var'] = info_var

    root_dir = sparentdir(srootdir())
    event_info['root_dir'] = root_dir
    data_dir = join_path([root_dir,'data'])
    if file_test(data_dir,/directory) eq 0 then file_mkdir, data_dir
    event_info['data_dir'] = data_dir
    plot_dir = join_path([root_dir,'plot'])
    if file_test(plot_dir,/directory) eq 0 then file_mkdir, plot_dir
    event_info['plot_dir'] = plot_dir
    event_info['file'] = join_path([data_dir,event_id+'_event_info.tplot'])
    if file_test(event_info['file']) eq 1 then begin
        tplot_restore, filename=event_info['file']
        event_info = get_var_data(info_var)
    endif

    ; time range for loading E/B fields and positions.
    time_range_long = time_double(['2014-08-28/09:25','2014-08-28/11:25']) ; E/B fields are bad before 09:25 for THD.

    ; Probes.
    all_probes = dictionary()
    all_probes['rbsp'] = ['b']
    all_probes['lanl'] = ['LANL-01A','1991-080','1994-084','LANL-97A','LANL-04A','LANL-02A']
    all_probes['goes'] = ['13','15']
    all_probes['themis'] = ['a','d','e']
    event_info['all_probes'] = all_probes

    probe_abbrev = dictionary()
    probe_abbrev['rbsp'] = 'rbsp'
    probe_abbrev['lanl'] = ''
    probe_abbrev['goes'] = 'g'
    probe_abbrev['themis'] = 'th'
    event_info['probe_abbrev'] = probe_abbrev

    constant = dictionary()
    constant['deg'] = 180d/!dpi
    constant['rad'] = !dpi/180d
    constant['re'] = 6378d
    event_info['constant'] = constant

    event_info['rgb'] = sgcolor(['red','green','blue'])


;---Plot settings.
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero'
    tplot_options, 'ystyle', 1
    tplot_options, 'yticklen', -0.02
    tplot_options, 'xticklen', -0.03


;---Load data file.
    if keyword_set(reload_all) then begin
        reload_ephemeris = 1
        reload_bfield = 1
        reload_efield = 1
        reload_injection = 1
        reload_model_info = 1
        reload_aurora = 1
        reload_gmag = 1
        recalc_gmag_info = 1
        file_delete, event_info['file'], /allow_nonexistent
        store_data, '*', /delete
    endif



;---Load injection data.
    injection = dictionary()
    injection['time_range'] = time_double(['2014-08-28/09:30','2014-08-28:11:00'])
    injection['probe'] = all_probes[['rbsp','lanl','goes','themis']]    ; slicing a dictionary returns a dictionary
    injection['energy_range'] = [0d,1000]   ; keV.
    injection['pitch_angle'] = 90d  ; deg.
    injection['species'] = dictionary('e','electron','h','proton')
    injection['file'] = join_path([data_dir,event_id+'_injection.tplot'])
    injection['var'] = []
    foreach key, injection['probe'].keys() do begin
        probes = (injection['probe'])[key]
        prefix = probe_abbrev[key]
        foreach species, injection['species'].keys() do begin
            injection['var'] = [injection['var'],prefix+probes+'_'+'kev_'+species+'_flux']
        endforeach
    endforeach

    _2014_0828_10_load_injection, injection, reload=reload_injection, event_info=event_info


;---Load ephemeris, r_gsm, mlt, mlon, mlat.
    ephemeris = dictionary()
    ephemeris['time_range'] = time_range_long
    ephemeris['probe'] = all_probes[['rbsp','lanl','goes','themis']]
    ephemeris['file'] = join_path([data_dir,event_id+'_ephemeris.tplot'])
    ephemeris['var'] = []
    foreach key, ephemeris['probe'].keys() do begin
        probes = (injection['probe'])[key]
        prefix = probe_abbrev[key]
        foreach probe, probes do begin
            ephemeris['var'] = [ephemeris['var'],prefix+probe+'_'+['r_gsm','mlat','mlon','mlt']]
        endforeach
    endforeach

    _2014_0828_10_load_ephemeris, ephemeris, reload=reload_ephemeris, event_info=event_info


;---Calc drift period.
    drift_period_info = dictionary($
        'species', ['h','e'], $
        'test_times', time_double(['2014-08-28/10:00','2014-08-28/11:00']), $
        'e', dictionary($
            'energy_range', [50.,500], $
            'probes', ['thd','g15','1991-080','g13','rbspb','1994-084']), $
        'h', dictionary($
            'energy_range', [70.,500], $
            'probes', ['thd','LANL-01A','LANL-02A','LANL-04A','LANL-97A']), $
        'file', join_path([data_dir,event_id+'_drift_period.tplot']))
    save_vars = list()
    foreach species, drift_period_info.species do begin
        the_info = drift_period_info[species]
        probes = the_info.probes
        drift_vars = probes+'_'+species+'_drift_period'
        save_vars.add, drift_vars, /extract
    endforeach
    drift_period_info['var'] = save_vars.toarray()

    _2014_0828_10_calc_drift_period, drift_period_info, reload=reload_drift_period, event_info=event_info


;---Load model info.
    model_info = dictionary()
    model_info['time_range'] = time_range_long
    model_info['probe'] = all_probes[['rbsp','lanl','goes','themis']]
    model_info['file'] = join_path([data_dir,event_id+'_model_info.tplot'])
    model_info['models'] = ['t89','t96','t01','t04s']
    model_info['trace_dir'] = -1
    model_info['var'] = []
    tvars = ['bmod_gsm','fpt_gsm','fmlat','fmlon','fmlt','c0map']
    foreach key, model_info['probe'].keys() do begin
        probes = (injection['probe'])[key]
        prefix = probe_abbrev[key]
        foreach probe, probes do begin
            foreach model, model_info['models'] do begin
                model_info['var'] = [model_info['var'],prefix+probe+'_'+tvars+'_'+model]
            endforeach
        endforeach
    endforeach

    _2014_0828_10_load_model_info, model_info, reload=reload_model_info, event_info=event_info


;---Load B field.
    bfield = dictionary()
    bfield['time_range'] = time_range_long
    bfield['probe'] = all_probes[['rbsp','goes','themis']]
    bfield['file'] = join_path([data_dir,event_id+'_bfield.tplot'])
    bfield['model_chosen'] = 't89'
    bfield['var'] = []
    tvars = ['b_gsm','db_gsm','b0_gsm','cmap','dbmag']
    foreach key, bfield['probe'].keys() do begin
        probes = (bfield['probe'])[key]
        prefix = probe_abbrev[key]
        foreach probe, probes do begin
            bfield['var'] = [bfield['var'],prefix+probe+'_'+tvars]
        endforeach
    endforeach

    _2014_0828_10_load_bfield, bfield, reload=reload_bfield, event_info=event_info


;---Load E field.
    efield = dictionary()
    efield['time_range'] = time_range_long
    efield['probe'] = all_probes[['rbsp','themis']]
    efield['file'] = join_path([data_dir,event_id+'_efield.tplot'])
    efield['var'] = []
    foreach key, efield['probe'].keys() do begin
        probes = (efield['probe'])[key]
        prefix = probe_abbrev[key]
        efield['var'] = [efield['var'],prefix+probes+'_e_gsm']
    endforeach

    _2014_0828_10_load_efield, efield, reload=reload_efield, event_info=event_info


;---Load pflux.
    pflux = dictionary()
    pflux['scale_info'] = {s0:6d, s1:600, dj:1d/8, ns:0d}
    pflux['probes'] = ['rbspb','th'+['a','d','e']]
    vars = list()
    foreach probe, pflux.probes do begin
        vars.add, probe+'_'+['e_fac','db_fac','pf_fac','e_mag_psd_cwt','db_mag_psd_cwt','pf_fac_mor_spec_1'], /extract
    endforeach
    pflux['var'] = vars.toarray()
    pflux['file'] = join_path([data_dir,event_id+'_pflux.tplot'])
    pflux['psd_time_range'] = time_double(['2014-08-28/10:00','2014-08-28/10:30'])

    _2014_0828_10_calc_pflux, pflux, reload=reload_pflux, event_info=event_info


;---Load aurora.
    aurora = dictionary()
    aurora['time_range'] = time_double(['2014-08-28/10:00','2014-08-28/10:30'])
    aurora['sites'] = ['whit','fsim']
    min_elevs = [5,20]
    foreach site, aurora['sites'], ii do begin
        aurora[site] = dictionary(themis_read_mlonimg_default_site_info(site),/extract)
        (aurora[site])['min_elev'] = min_elevs[ii]
    endforeach
    aurora['mlon_range'] = [-100,-55]
    aurora['mlat_range'] = [60,75]
    aurora['merge_method'] = 'max_elev'
    aurora['file'] = join_path([data_dir,event_id+'_aurora.tplot'])
    aurora['var'] = 'thg_mlonimg'

    _2014_0828_10_load_aurora, aurora, reload=reload_aurora, event_info=event_info


;---Load gmag data.
    gmag = dictionary()
    gmag['time_range'] = time_double(['2014-08-28/09:30','2014-08-28/11:00'])
    gmag['mlon_range'] = [-170,10]
    gmag['mlat_range'] = [55,80]
    gmag['component'] = ['h','d','z'] ; ['h','d','z']   north, east, down.
    gmag['file'] = join_path([data_dir,event_id+'_gmag.tplot'])
    gmag['var'] = 'thg_db'+gmag['component']

    _2014_0828_10_load_gmag, gmag, reload=reload_gmag, event_info=event_info


;---Load weygand.
    weygand = dictionary()
    weygand['time_range'] = time_double(['2014-08-28/08:00','2014-08-28/24:00'])
    weygand['mlon_range'] = [50.,-100]
    weygand['mlat_range'] = [60.,80]
    weygand['ewo_mlat_range'] = [60.,70]
    weygand['mlt_range'] = [-3.,9]
    weygand['mlon_binsize'] = 4
    weygand['mlat_binsize'] = 2
    weygand['var'] = ['thg_j_ver','thg_j_up_ewo']
    weygand['file'] = join_path([data_dir,event_id+'_weygand.tplot'])

    _2014_0828_10_load_weygand, weygand, reload=reload_weygand, event_info=event_info



;;---Calculate gmag info.
;    gmag_info = dictionary()
;    gmag_info['base_time'] = time_double('2014-08-28/10:08')
;    gmag_info['time_range'] = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
;    gmag_info['fac_method'] = 'absolute'
;    gmag_info['fac_threshold'] = 50 ; nT.
;    gmag_info['fac_db_threshold'] = 100 ; nT. For sites have small dB, use the first point. Otherwise, use the last point.
;    gmag_info['slope_width'] = 30.  ; sec.
;    info_var = event_info['var']
;    event_info = get_var_data(info_var)
;    recalc = ~event_info.haskey('gmag_info')
;    if keyword_set(recalc_gmag_info) then recalc = 1
;    if recalc then _2014_0828_10_calc_gmag_info, gmag_info, event_info=event_info


;---Cleanup.
    event_info = get_var_data(info_var)
    lprmsg, 'Done ...'
end

_2014_0828_10_load_data, event_info=event_info, $
    reload_bfield=0, $          ; Need model info. Load B in GSM, dB and B0.
    reload_model_info=0, $
    reload_efield=0, $          ; Need some settings.
    reload_pflux=0, $
    reload_aurora=0, $          ; Need some settings.
    reload_gmag=0, $              ; Need some settings.
    recalc_gmag_info=0, $
    reload_drift_period=0, $
    reload_weygand=0

end
