;+
; Calculate the # of spacecraft available in pre-midnight and post-midnight
; Calculate the # of good triad in pre-midnight and post-midnight
;-

function normalize_data, zzs, mlt, smooth_width=smooth_width, zrange=range, ct=ct, reverse_ct=reverse_ct

    if keyword_set(smooth_width) then zzs -= smooth(zzs,smooth_width, /edge_truncate, /nan)
    mean = mean(zzs,/nan)
    zzs = (zzs-mean)*sqrt(abs(mlt)+1)+mean

    linear_range = 1
    signs = sign(zzs)
    data = alog10(abs(zzs)>linear_range)*signs
    the_range = [-1,1]*alog10(max(abs(range)))
    index_color = bytscl(data, min=the_range[0], max=the_range[1])
    if keyword_set(reverse_ct) then index_color = 256-index_color

    if n_elements(ct) eq 0 then ct = 70
    loadct, ct, rgb_table=rgb0
    rgb0 = float(rgb0)
    true_colors = 256L*(256L*rgb0[index_color,2]+rgb0[index_color,1])+rgb0[index_color,0]

    return, true_colors
end

pro azim_df_analysis_storm, project=project

    small_angle_limit = 15. ; deg.
    triad_settings = dictionary($
        'small_angle_limit', small_angle_limit, $
        'angle_range', [small_angle_limit,180-small_angle_limit])

    ; Cut off data outside the range.
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9
    if n_elements(xsm_range) ne 2 then xsm_range = [-20,4]
    if n_elements(dis_range) ne 2 then dis_range = minmax(abs(xsm_range))


    if n_elements(project) eq 0 then project = azim_df_load_project()
    if ~project.haskey('done_filter_storm') then azim_df_filter_storm, project=project
    candidate = project.candidate
    storm_ids = candidate.keys()
    storm_ids = storm_ids.sort()

test = 0
    if keyword_set(test) then storm_ids = '2014_0828_01'


    filter_storm_setting = project.filter_storm_setting
    mlt_limit = filter_storm_setting.mlt_limit
    regions = list()
    foreach name, ['pre','post'] do begin
        the_mlt_range = (name eq 'pre')? [-1,0]: [0,1]
        the_mlt_range *= mlt_limit
        regions.add, dictionary('name', name+'_midn', $
            'mlt', the_mlt_range)
    endforeach


    time_step = 90.    ; sec.
    foreach event_id, storm_ids do begin
        tinfo = candidate[event_id]
        time_range = tinfo.time_range

        the_key = 'data_file_suffix'
        if ~tinfo.haskey(the_key) then tinfo[the_key] = event_id+'_data.tplot'
        data_file = join_path([project.data_dir,tinfo[the_key]])
        if file_test(data_file) eq 0 then azim_df_load_data, project=project, event_id=event_id
        del_data, '*'
        tplot_restore, filename=data_file

        the_key = 'available_probes'
        if ~tinfo.haskey(the_key) then begin
            del_data, '*'
            tplot_restore, filename=data_file
            mission_probes = tnames('*_mlt')
            foreach probe, mission_probes, jj do mission_probes[jj] = strmid(probe,0,strpos(probe,'_'))
            tinfo[the_key] = mission_probes
        endif else mission_probes = tinfo[the_key]

        nmission_probe = n_elements(mission_probes)
        probe_colors = smkarthm(0.,254,nmission_probe,'n')
        ct = 40
        for jj=0, nmission_probe-1 do probe_colors[jj] = sgcolor(probe_colors[jj], ct=ct)


    ;---Detrend.
        detrend_width = 80*60/10
        foreach mission_probe, mission_probes, jj do begin
            the_var = mission_probe+'_db_tilt'
            get_data, the_var, times, tilt
            tilt -= smooth(tilt, detrend_width, /nan, /edge_truncate)
            store_data, the_var, times, tilt
            options, the_var, 'labels', strupcase(mission_probe)
        endforeach


    ;---Downsample.
        common_times = make_bins(time_range, time_step)
        vars = tnames('*')
        foreach var, vars do begin
            get_data, var, times
            if n_elements(times) eq 0 then continue
            interp_time, var, common_times
        endforeach


    ;---Collect instantaneous status on how many s/c in one region:
    ;   1. Within magnetopause
    ;   2. Have data
        ntime = n_elements(common_times)
        nregion = n_elements(regions)
        data_flags = bytarr(ntime,nmission_probe,nregion)
        key = 'mlt'
        foreach mission_probe, mission_probes, jj do begin
            mlt = get_var_data(mission_probe+'_'+key, at=common_times)
            foreach region, regions, kk do begin
                index = lazy_where(mlt, '[]', region[key], count=count)
                if count ne 0 then data_flags[index,jj,kk] = 1
            endforeach
            rsm = get_var_data(mission_probe+'_r_sm', at=common_times)
            rgsm = cotran(rsm, common_times, 'sm2gsm')
            magn_flags = check_if_in_magn(rgsm)
            index = where(magn_flags eq 0, count)
            if count ne 0 then data_flags[index,jj,*] = 0
            dis = snorm(rgsm)
            index = where(dis lt min(dis_range), count)
            if count ne 0 then data_flags[index,jj,*] = 0

            tilt = get_var_data(mission_probe+'_db_tilt', at=common_times)
            index = where(finite(tilt,/nan), count)
            if count ne 0 then data_flags[index,jj,*] = 0
        endforeach

        foreach region, regions, jj do begin
            the_var = region['name']+'_data_flag'
            store_data, the_var, common_times, reform(data_flags[*,*,jj])
            add_setting, the_var, /smart, {$
                yrange: [-0.2,1.2], $
                ytitle: strjoin(strsplit(the_var,'_',/extract),' '), $
                labels: strupcase(mission_probes), $
                colors: probe_colors }
        endforeach


    ;---Collect instantaneous status on how many good triads
        foreach mission_probe, mission_probes do begin
            the_var = mission_probe+'_r_sm'
            rsm = get_var_data(the_var, at=common_times)
            store_data, mission_probe+'_x_sm', common_times, rsm[*,0]
            store_data, mission_probe+'_y_sm', common_times, rsm[*,1]
        endforeach

        ndim = 3
        angle_range = triad_settings.angle_range
        foreach region, regions, jj do begin
            num_good_triads = intarr(ntime)
            the_var = region['name']+'_data_flag'
            data_flags = get_var_data(the_var, at=common_times)
            foreach time, common_times, kk do begin
                index = where(data_flags[kk,*] eq 1, count)
                if count le ndim then continue
                probes = mission_probes[index]
                triads = choose_from(probes, ndim)
                ntriad = n_elements(triads)
                the_count = intarr(ntriad)
                xxs = fltarr(ndim)
                yys = fltarr(ndim)
                angles = fltarr(ndim)
                foreach triad, triads, ll do begin
                    for mm=0,ndim-1 do xxs[mm] = get_var_data(triad[mm]+'_x_sm', at=time)
                    for mm=0,ndim-1 do yys[mm] = get_var_data(triad[mm]+'_y_sm', at=time)
                    for mm=0,ndim-1 do begin
                        xxs = shift(xxs,1)
                        yys = shift(yys,1)
                        angles[mm] = sang($
                            [xxs[1]-xxs[0],yys[1]-yys[0]], $
                            [xxs[2]-xxs[0],yys[2]-yys[0]], /degree)
                    endfor
                    index = lazy_where(angles, '[]', angle_range, count=count)
                    if count eq ndim then the_count[ll] = 1 ; A good triad.
                endforeach
                num_good_triads[kk] = total(the_count)
            endforeach
            the_var = region['name']+'_triad_flag'
            store_data, the_var, common_times, num_good_triads
            add_setting, the_var, /smart, {$
                yrange: minmax(num_good_triads)+[-1,1]*0.2, $
                ytitle: strjoin(strsplit(the_var,'_',/extract),' ') }
        endforeach

    ;---Prepare to plot.
        if keyword_set(test) then plot_file = 0 else plot_file = join_path([project.plot_dir,'selected_storms','fig_overview_'+event_id+'.pdf'])

    ;---Figure out the size of the plot.
        secofday = 86400d
        duration = total(time_range*[-1,1])
        nday = duration/secofday
        ysize = 6.
        aspect_ratio = 2
        xsize = ysize*aspect_ratio*nday ; normalize by duration.
        sgopen, plot_file, xsize=xsize, ysize=ysize

        npanel = 4
        poss = sgcalcpos(npanel, ypans=[1,2,2,1], tmargin=2, rmargin=8, lmargin=12, xchsz=xchsz, ychsz=ychsz)
        full_ychsz = 0.7
        half_ychsz = 0.35
        label_xshift = 2
        label_yshift = full_ychsz
        label_size = 0.7
        constant_linestyle = 1

    ;---The common xrange.
        xrange = time_range
        xticklen_chsz = -0.15
        yticklen_chsz = -0.30
        xminor = 8  ; hour.
        xtickv = make_bins(xrange, 3600*xminor, /inner) ; make time line up at hours.
        xticks = n_elements(xtickv)-1
        xtickn = strarr(xticks+1)
        for jj=0, xticks do begin
            the_time = xtickv[jj]
            xtickn[jj] = time_string(the_time,tformat='hh:mm')
            date = time_string(the_time,tformat='YYYY-MM-DD')
            if jj eq 0 then begin
                xtickn[jj] += '!C'+date
                continue
            endif
            if the_time mod secofday ne 0 then continue
            xtickn[jj] += '!C'+date
        endfor

    ;---Dst and AE.
        tpos = poss[*,0]
        xtickformat='(A1)'
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        vars = ['dst','ae']
        if tnames(vars[0]) eq '' then omni_read_index, time_range
        options, vars, 'yticklen', yticklen
        options, vars, 'ystyle', 1
        options, vars, 'yminor', 5
        ystep = [50,500]
        constant = [0,500]
        ytitles = ['Dst','AE']+' (nT)'
        foreach var, vars, jj do begin
            get_data, var, times, data
            ytickv = make_bins(data, ystep[jj])
            yticks = n_elements(ytickv)-1
            options, var, 'ytitle', ytitles[jj]
            options, var, 'ytickv', ytickv
            options, var, 'yticks', yticks
            options, var, 'yrange', minmax(ytickv)
            options, var, 'constant', constant[jj]
        endforeach

        get_data, 'dst', xxs, yys, limits=lim
        yrange = lim.yrange
        ytickv = lim.ytickv
        yticks = lim.yticks
        yminor = lim.yminor
        ytitle = lim.ytitle
        constant = lim.constant
        dst_color = sgcolor('black')
        plot, xrange, yrange, position=tpos, $
            xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
            /nodata, /noerase
        oplot, xxs, yys, color=dst_color
        oplot, xrange, constant+[0,0], linestyle=constant_linestyle
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        get_data, 'ae', xxs, yys, limits=lim
        yrange = lim.yrange
        ytickv = lim.ytickv
        yticks = lim.yticks
        yminor = lim.yminor
        ytitle = lim.ytitle
        constant = lim.constant
        ae_color = sgcolor('red')
        axis, yaxis=1, /save, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            color=ae_color
        oplot, xxs, yys, color=ae_color
        oplot, xrange, constant+[0,0], linestyle=constant_linestyle, color=color

        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'a. Dst/'
        ;ty = tpos[3]-label_yshift*ychsz-full_ychsz*ychsz
        xyouts, tx,ty,/normal, '           AE', color=ae_color



    ;---Common z-range, colorbar, symbol.
        data_var_suffix = 'db_tilt'
        ztitle = 'Log!D10!N[Tilt angle - T89 (deg)]'
        tilt_range = [-1,1]*64
        zrange = tilt_range
        ztickv = [-64,-16,-4,0,4,16,64]
        ztickn = string(ztickv,format='(I0)')
        ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
        zrange = alog10(abs(zrange))*[-1,1]
        zticks = n_elements(ztickv)-1

        ; Symbol.
        psym = 8
        symsize = 0.5
        usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
        ct = 70
        reverse_ct = 1

        ; About transform data.
        smooth_window = 80*60   ; sec.
        prefixes = tnames('*_mlt')
        foreach prefix, prefixes, jj do prefixes[jj] = get_prefix(prefixes[jj])

    ;---EWOgram.
        tpos = poss[*,1]
        xtickformat = '(A1)'
        ytitle = 'MLT (hr)'
        yrange = mlt_range
        yminor = 3
        ytickv = make_bins(yrange,yminor, /inner)
        yticks = n_elements(ytickv)-1
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        constant = 0

        plot, xrange, yrange, position=tpos, $
            xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
            /nodata, /noerase

        pos_var_suffix = 'mlt'
        available_prefixes = prefixes
        flags = list()
        foreach prefix, prefixes do begin
            pos_var = prefix+pos_var_suffix
            data_var = prefix+data_var_suffix
            get_data, data_var, xxs, zzs
            get_data, pos_var, xxs, yys
            if n_elements(zzs) ne n_elements(xxs) then begin
                errmsg = handle_error(data_var+': inconsistent data and time ...')
                stop
            endif
            mlt = get_var_data(prefix+'mlt')
            smooth_width = 0
            zzs = normalize_data(zzs, mlt, smooth_width=smooth_width, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

            ; Remove data out of range.
            get_data, prefix+'r_sm', times, r_sm
            xsm = r_sm[*,0]
            dis = snorm(r_sm)
            ; Exclude data outside the MLT range.
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then xxs[index] = !values.f_nan
            ; Exclude data outside the x-range.
            index = lazy_where(xsm, '][', xsm_range, count=count)
            if count ne 0 then xxs[index] = !values.f_nan
            ; Exclude data ouside the distance range.
            index = lazy_where(dis, '][', dis_range, count=count)
            if count ne 0 then xxs[index] = !values.f_nan
            ; Exclude data outside magnetopause.
            magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
            index = where(magn_flags eq 0, count)
            if count ne 0 then xxs[index] = !values.f_nan

            index = where(finite(xxs), count)
            the_flag = (count gt 0)
            flags.add, the_flag
            if count eq 0 then continue

            nxx = count
            xxs = xxs[index]
            yys = yys[index]
            zzs = zzs[index]
            for jj=0, nxx-1 do plots, xxs[jj], yys[jj], color=zzs[jj], psym=psym, symsize=symsize

        endforeach
        oplot, xrange, constant+[0,0], linestyle=constant_linestyle

        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'b. UT-'+strmid(ytitle,0,strpos(ytitle,' '))

        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
        sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks

        flags = flags.toarray()
        index = where(flags eq 1, count)
        available_prefixes = available_prefixes[index]

    ;---Keogram.
        tpos = poss[*,2]
        xtickformat = '(A1)'
        ytitle = 'X!USM!N (Re)'
        yrange = xsm_range
        yminor = 5
        ytickv = make_bins(yrange,yminor, /inner)
        yticks = n_elements(ytickv)-1
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        plot, xrange, yrange, position=tpos, $
            xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
            /nodata, /noerase

        pos_var_suffix = 'r_sm'
        foreach prefix, prefixes do begin
            pos_var = prefix+pos_var_suffix
            data_var = prefix+data_var_suffix
            get_data, data_var, xxs, zzs
            get_data, pos_var, xxs, yys & yys = yys[*,0]
            if n_elements(zzs) ne n_elements(xxs) then begin
                errmsg = handle_error(data_var+': inconsistent data and time ...')
                stop
            endif
            mlt = get_var_data(prefix+'mlt')
            smooth_width = smooth_window/time_step
            zzs = normalize_data(zzs, mlt, smooth_width=smooth_width, ct=ct, reverse_ct=reverse_ct, zrange=tilt_range)

            ; Remove data out of range.
            get_data, prefix+'r_sm', times, r_sm
            xsm = r_sm[*,0]
            dis = snorm(r_sm)
            ; Exclude data outside the MLT range.
            index = lazy_where(mlt, '][', mlt_range, count=count)
            if count ne 0 then xxs[index] = !values.f_nan
            ; Exclude data outside the x-range.
            index = lazy_where(xsm, '][', xsm_range, count=count)
            if count ne 0 then xxs[index] = !values.f_nan
            ; Exclude data ouside the distance range.
            index = lazy_where(dis, '][', dis_range, count=count)
            if count ne 0 then xxs[index] = !values.f_nan
            ; Exclude data outside magnetopause.
            magn_flags = check_if_in_magn(cotran(r_sm, times, 'sm2gsm'))
            index = where(magn_flags eq 0, count)
            if count ne 0 then xxs[index] = !values.f_nan

            index = where(finite(xxs), count)
            if count eq 0 then continue

            nxx = count
            xxs = xxs[index]
            yys = yys[index]
            zzs = zzs[index]
            for jj=0, nxx-1 do plots, xxs[jj], yys[jj], color=zzs[jj], psym=psym, symsize=symsize
        endforeach

        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'c. UT-'+strmid(ytitle,0,strpos(ytitle,' '))

        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
        sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


    ;---Print available_probes.
        msg = ['Available spacecraft: ']
        prefixes = sort_uniq(available_prefixes)
        foreach prefix, prefixes do begin
            mission_probe = strmid(prefix,0,strpos(prefix,'_'))
            mission_info = resolve_probe(mission_probe)
            msg = [msg,strupcase(mission_info.short_name)]
        endforeach
        if n_elements(msg) eq 1 then msg = [msg,'none']
        msg = msg[0]+strjoin(msg[1:*],',')
        tpos = poss[*,0]
        tx = tpos[0]
        ty = tpos[3]+ychsz*half_ychsz
        xyouts, tx,ty,/normal, msg


    ;---Add triad count.
        vars = tnames('*_triad_flag')
        labels = vars
        foreach var, labels, jj do labels[jj] = strupcase(strmid(var,0,strpos(var,'_')))
        colors = sgcolor(['red','blue'])
        the_var = event_id+'_triad_flag'
        stplot_merge, vars, newname=the_var, labels=labels, colors=colors
        get_data, the_var, xxs, yys
        tpos = poss[*,3]
        xtickformat = ''
        ytitle = '(#)'
        yrange = minmax(yys)
        yminor = 5
        ytickv = make_bins(yrange,yminor, /inner)
        yticks = n_elements(ytickv)-1
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        plot, xrange, yrange, position=tpos, $
            xstyle=9, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=9, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        dy = (tpos[3]-tpos[1])/(n_elements(labels)+1)
        tys = smkarthm(tpos[3]-dy,tpos[1]+dy,-dy,'dx')-half_ychsz*ychsz
        tx = tpos[2]+xchsz
        foreach label, labels, jj do begin
            oplot, xxs, yys[*,jj], color=colors[jj]
            xyouts, tx, tys[jj], /normal, labels[jj], color=colors[jj]
        endforeach
        
        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'd. No. triad'

        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase


        if keyword_set(test) then stop
        sgclose


    endforeach
    update_project, project

end
