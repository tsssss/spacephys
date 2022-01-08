;+
; Map a probe to equatorial plane.
;-


pro map_probe_to_equatorial_plane_t89, r_var

    get_data, r_var, times, r_vec, limits=lim
    coord = get_setting(r_var, 'coord')
    r_gsm = cotran(r_vec, times, strlowcase(coord)+'2gsm')
    rf_gsm = r_gsm

    time_range = minmax(times+[-1,1]*3600)
    time_range = time_range-(time_range mod 3600)+[0,1]*3600

    r0 = 4.
    ntime = n_elements(times)
    ndim = 3
    r_eq_gsm = fltarr(ntime,ndim)
    foreach time, times, time_id do begin
        tilt = geopack_recalc(time)
        rr = r_gsm[time_id,*]
        flines = rr
        par = 2
        foreach dir, [1,-1] do begin
            geopack_trace, rr[0],rr[1],rr[2], dir, par, xf,yf,zf, fline=fline, r0=r0, t89=1
            flines = (dir eq -1)? [flines,fline]: [flines,reverse(fline,1)]
        endforeach
        flines = flines[2:-2,*]     ;
        tmp = max(snorm(flines), index)
        ; In some cases, flines diverge from the middle. This means the mapping isn't working.
        if index eq 0 or index eq n_elements(flines[*,0])-1 then begin
            r_eq_gsm[time_id,*] = !values.f_nan
        endif else begin
            dfline = deriv(smooth(snorm(flines),10))
            nodes = find_node(dfline)
            ; nontrivial nodes often appear when <6 Re.
            if n_elements(nodes) gt 1 and max(snorm(flines)) ge 6 then begin
                stop
                r_eq_gsm[time_id,*] = !values.f_nan
            endif else begin
                r_eq_gsm[time_id,*] = flines[index,*]
            endelse
        endelse
    endforeach
    r_eq_vec = cotran(r_eq_gsm, times, 'gsm2'+strlowcase(coord))
    store_data, r_var+'_eq', times, r_eq_vec, limits=lim

end

pro map_probe_to_equatorial_plane_t96, r_var

    get_data, r_var, times, r_vec, limits=lim
    coord = get_setting(r_var, 'coord')
    r_gsm = cotran(r_vec, times, strlowcase(coord)+'2gsm')
    rf_gsm = r_gsm

    time_range = minmax(times+[-1,1]*3600)
    time_range = time_range-(time_range mod 3600)+[0,1]*3600
    if check_if_update('t96_var', time_range) then sgeopack_par, time_range, 't96'
    
    r0 = 4.
    ntime = n_elements(times)
    ndim = 3
    r_eq_gsm = fltarr(ntime,ndim)
    foreach time, times, time_id do begin
        tilt = geopack_recalc(time)
        rr = r_gsm[time_id,*]
        flines = rr
        par = get_var_data('t96_par', at=time)
        foreach dir, [1,-1] do begin
            geopack_trace, rr[0],rr[1],rr[2], dir, par, xf,yf,zf, fline=fline, r0=r0, t96=1
            flines = (dir eq -1)? [flines,fline]: [flines,reverse(fline,1)]
        endforeach
        flines = flines[2:-2,*]     ; 
        tmp = max(snorm(flines), index)
        ; In some cases, flines diverge from the middle. This means the mapping isn't working.
        if index eq 0 or index eq n_elements(flines[*,0])-1 then begin
            r_eq_gsm[time_id,*] = !values.f_nan
        endif else begin
            dfline = deriv(smooth(snorm(flines),10))
            nodes = find_node(dfline)
            ; nontrivial nodes often appear when <6 Re.
            if n_elements(nodes) gt 1 and max(snorm(flines)) ge 6 then begin
                stop
                r_eq_gsm[time_id,*] = !values.f_nan                
            endif else begin
                r_eq_gsm[time_id,*] = flines[index,*]
            endelse
        endelse
    endforeach
    r_eq_vec = cotran(r_eq_gsm, times, 'gsm2'+strlowcase(coord))
    store_data, r_var+'_eq', times, r_eq_vec, limits=lim

end



    test_mapping = 0
    if keyword_set(test_mapping) then begin
        time_range = time_double(['2014-08-28/11:42','2014-08-28/11:43'])
        probe = 'b'
        rbsp_read_orbit, time_range, probe=probe, coord='gsm'
        prefix = 'rbsp'+probe+'_'
        r_var = prefix+'r_gsm'
        r_gsm = get_var_data(r_var, times=times, limits=lim)
        r_sm = cotran(r_gsm, times, 'gsm2sm')
        lim.coord = 'SM'
        store_data, prefix+'r_sm', times, r_sm, limits=lim
        r_var = prefix+'r_sm'
        map_probe_to_equatorial_plane, r_var
    endif

    model = 't89'
    test = 1
    if n_elements(mlt) eq 0 then begin
        if n_elements(project) eq 0 then project = azim_df_load_project()
        events = azim_df_find_dfgroup(project=project)
    
        all_df_r_sm = list()
        all_df_eq_r_sm = list()
        foreach event, events do begin
            foreach df, event.df_list do begin
                r_var = df.probe+'_r_sm'
                store_data, r_var, df.obs_time, transpose(df.obs_r_sm)
                add_setting, r_var, dictionary('coord', 'SM')
                call_procedure, 'map_probe_to_equatorial_plane_'+model, r_var
                r_eq_sm = get_var_data(r_var+'_eq')
                df['eq_r_sm'] = reform(r_eq_sm)
                df['eq_mlt'] = pseudo_mlt(r_eq_sm)
                all_df_r_sm.add, df.obs_r_sm
                all_df_eq_r_sm.add, df.eq_r_sm
            endforeach
        endforeach
        
        all_df_r_sm = all_df_r_sm.toarray()
        all_df_eq_r_sm = all_df_eq_r_sm.toarray()
        
        rxy = snorm(all_df_r_sm[*,0:1])
        eq_rxy = snorm(all_df_eq_r_sm[*,0:1])
        mlt = pseudo_mlt(all_df_r_sm)
        eq_mlt = pseudo_mlt(all_df_eq_r_sm)
    end
    
    plot_file = join_path([homedir(),'check_uncertainty_of_equatorial_mapping_'+model+'.pdf'])
    if keyword_set(test) then plot_file = 0
    symsize = 0.5
    poss = panel_pos(plot_file, pansize=[2,2], xpans=[1,1], fig_size=fig_size)
    sgopen, plot_file, xsize=fig_size[0], ysize=fig_size[1]
    
    xrange = [-1,1]*10
    yrange = [-1,1]*10
    xtitle = 'MLT (hr)'
    ytitle = 'Equatorial footpoint!C'+strupcase(model)+' MLT (hr)'
    tpos = poss[*,0]
    plot, mlt, eq_mlt, psym=1, symsize=symsize, $
        xstyle=1, xrange=xrange, xtitle=xtitle, $
        ystyle=1, yrange=yrange, ytitle=ytitle, $
        position=tpos, /noerase
    plots, xrange, yrange, linestyle=1, color=sgcolor('silver')
    

    xrange = [0,1]*40
    yrange = [0,1]*40
    xtitle = 'Rxy (Re)'
    ytitle = 'Equatorial footpoint!C'+strupcase(model)+' Rxy (Re)'
    tpos = poss[*,1]
    plot, rxy, eq_rxy, psym=1, symsize=symsize, $
        xstyle=1, xrange=xrange, xtitle=xtitle, $
        ystyle=1, yrange=yrange, ytitle=ytitle, $
        position=tpos, /noerase
    plots, xrange, yrange, linestyle=1, color=sgcolor('silver')

    if keyword_set(test) then stop
    sgclose
end
