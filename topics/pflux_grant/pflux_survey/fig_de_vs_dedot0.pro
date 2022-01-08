;+
; To compare dE and dEdot0 when the latter can be determined.
;-

test = 0

    project = pflux_survey_load_project()
    probes = ['a','b']
    plot_file = join_path([project.plot_dir,'figures','fig_de_vs_dedot0_vs_mlt.pdf'])
    if keyword_set(test) then plot_file = test
    mission_time_range = project.time_range

    de_vars = ['de','dedot0']
    fac_labels = ['b','w','o']
    fac_full_labels = ['parallel','westward','outward']
    ndim = 3
    de_step = 16*60.    ; To select data every minute.
    secofday = constant('secofday')

    ; Settings for edot0_flag.
    min_bw_ratio = 0.2
    edot0_flag_pad_time = 1800.


;---Check if data for plotting are saved.
    save_vars = list()
    data_file = join_path([project.data_dir,'pflux_survey_de_vs_dedot0.tplot'])
    ;tmp = pflux_grant_load_project()
    ;data_file2 = join_path([tmp.data_dir,'pflux_survey_de_vs_dedot0.tplot'])
    if file_test(data_file) eq 1 then tplot_restore, filename=data_file
    file_updated = file_test(data_file) eq 0

    foreach probe, probes do begin
        rbspx = 'rbsp'+probe
        prefix = rbspx+'_'

    ;---Load orbit data.
        r_gsm_var = prefix+'r_gsm'
        save_vars.add, r_gsm_var
        if check_if_update(r_gsm_var) then begin
            rbsp_read_orbit, mission_time_range, probe=probe
            file_updated = 1
        endif

    ;---Load E field.
        get_data, r_gsm_var, common_times
        ntime = n_elements(common_times)
        secofday = constant('secofday')
        nday = total(mission_time_range*[-1,1])/secofday
        days = smkarthm(mission_time_range[0],secofday,nday, 'x0')
        foreach de_var, de_vars do begin
            de_var = prefix+de_var+'_fac'
            de_save_var = de_var+'_interp'
            save_vars.add, de_save_var
            if ~check_if_update(de_save_var) then continue
            file_updated = 1
            de_data = fltarr(ntime,ndim)
            for day_id=0, nday-1 do begin
                time_range = days[day_id]+[0,secofday]
                pflux_grant_read_preprocessed_ebfield, time_range, probe=probe
                index = lazy_where(common_times,'[)', time_range)
                de_data[index,*] = (get_var_data(de_var))[0:*:de_step,*]
            endfor
            store_data, de_save_var, common_times, de_data
            short_name = (de_var eq 'de')? 'E': 'E!S!Udot0!R!N'
            add_setting, de_save_var, /smart, {$
                display_type: 'vector', $
                unit: 'mV/m', $
                short_name: short_name, $
                coord: '', $
                coord_labels: fac_labels}
        endforeach

    ;---Load Bw ratio.
        bw_var = prefix+'bw_ratio'
        bw_save_var = bw_var+'_interp'
        save_vars.add, bw_save_var
        if check_if_update(bw_save_var) then begin
            bw_data = fltarr(ntime)
            for day_id=0,nday-1 do begin
                time_range = days[day_id]+[0,secofday]
                pflux_grant_read_preprocessed_ebfield, time_range, probe=probe, id='bw_ratio'
                index = lazy_where(common_times,'[)', time_range)
                bw_data[index,*] = (get_var_data(bw_var))[0:*:de_step]
            endfor
            file_updated = 1

            store_data, bw_save_var, common_times, bw_data
            add_setting, bw_save_var, /smart, {$
                display_type: 'scalar', $
                unit: '#', $
                short_name: 'Bw/|B|'}
        endif

    ;---Index for good Edot0 data.
        edot0_flag_var = prefix+'edot0_flag'
        save_vars.add, edot0_flag_var
        if check_if_update(edot0_flag_var) then begin
            bw_ratio = get_var_data(bw_save_var)
            edot0_flag = abs(bw_ratio) le min_bw_ratio  ; 1 for bad data.
            index = where(edot0_flag eq 1, count)
            time_step = total(common_times[0:1]*[-1,1])
            drec = edot0_flag_pad_time/time_step
            ntime = n_elements(common_times)
            for ii=0,count-1 do begin
                i0 = (index[ii]-drec)>0
                i1 = (index[ii]+drec)<(ntime-1)
                edot0_flag[i0:i1] = 1
            endfor

            file_updated = 1
            store_data, edot0_flag_var, common_times, edot0_flag
            add_setting, edot0_flag_var, /smart, {$
                display_type: 'scalar', $
                unit: '#', $
                short_name: 'Edot0 flag', $
                min_bw_ratio: min_bw_ratio, $
                pad_time: edot0_flag_pad_time, $
                fieldnam: '1 for bad data'}
        endif
    endforeach

    save_vars = save_vars.toarray()
    if file_updated then tplot_save, save_vars, filename=data_file


;---Remove data when bw_ratio does not meet criteria.
    de_range = list()
    foreach probe, probes do begin
        prefix = 'rbsp'+probe+'_'
        get_data, prefix+'edot0_flag', common_times, dedot0_flag
        dedot0_fac = get_var_data(prefix+'dedot0_fac_interp')
        de_fac = get_var_data(prefix+'de_fac_interp')
        r_gsm = get_var_data(prefix+'r_gsm')
        de_range.add, minmax(de_fac)
        de_range.add, minmax(dedot0_fac)

        if ~check_if_update(prefix+'de_fac') then continue
        index = where(dedot0_flag eq 0)
        times = common_times[index]
        dedot0_fac = dedot0_fac[index,*]
        de_fac = de_fac[index,*]
        r_sm = cotran(r_gsm[index,*], common_times, 'gsm2sm')
        mlt = pseudo_mlt(r_sm)
        store_data, prefix+'de_fac', times, de_fac
        store_data, prefix+'dedot0_fac', times, dedot0_fac
        store_data, prefix+'mlt', times, mlt
    endforeach
    de_range = minmax(de_range.toarray())
    de_tick_step = 50.
    de_range = [-1,1]*ceil(max(abs(de_range))/de_tick_step)*de_tick_step
    de_tickv = make_bins(de_range,de_tick_step*2,/inner)
    de_range = [-1,1]*180
    de_tickv = make_bins(de_range,de_tick_step*2,/inner)
    de_ticks = n_elements(de_tickv)-1
    de_minor = 5
    de_unit = 'mV/m'
    de_log = 0
    short_name = 'E'
    de_labels = fac_full_labels
    de_ticklen = -0.02
    fig_letters = letters(ndim)


;---MLT settings.
    mlt_range = [-6,6]
    mlt_step = 1
    mlt_tick_step = 6
    mlt_centers = make_bins(mlt_range,mlt_step)
    mlt_vertices = [[mlt_centers-mlt_step*0.5],[mlt_centers+mlt_step*0.5]]
    nmlt_bin = n_elements(mlt_centers)

    ct = 6  ; circular color table.
    ncolor = n_elements(mlt_centers)
    color_min = 0+60
    color_max = 255-60
    index_colors = smkarthm(color_min, color_max, ncolor, 'n')
    colors = index_colors
    foreach color, index_colors, ii do colors[ii] = sgcolor(color, ct=ct)


;---Symbol settings.
    psyms = [6,7,8]
    npoint = 40
    tmp = findgen(npoint)/(npoint-1)*2*!dpi
    circ_x = cos(tmp)
    circ_y = sin(tmp)
    usersym, circ_x, circ_y

;---The plot.
    npanel = 4
    panel_size = 1.8
    xpad = 1.5
    sgopen, 0, xsize=1,ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    margins = [8,5,7,2]
    fig_xsize = panel_size*npanel+(xpad*(npanel-1)+total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_size+total(margins[[1,3]])*abs_ychsz
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz

    poss = sgcalcpos(1,npanel, xpad=xpad, ypad=0, margins=margins)
    poss[[0,2],npanel-1] += xchsz*(margins[2]-1.5)

;---Colorbar.
;    cbpos = poss[*,npanel-1]
;    cbpos[0] = cbpos[2]+xchsz*xpad
;    cbpos[2] = cbpos[0]+xchsz*1
;    ztickv = make_bins(mlt_range, mlt_tick_step)
;    zticks = n_elements(ztickv)
;    ztitle = 'MLT (hr)'
;    sgcolorbar, index_colors, ct=ct, position=cbpos, $
;        zrange=mlt_range, zticks=zticks, ztickv=ztickv, ztitle=ztitle

    cbpos = poss[*,npanel-1]
    cbpos[3] = cbpos[1]+ychsz*0.5
    ztick_step = mlt_tick_step*0.5
    ztickv = make_bins(mlt_range, ztick_step)
    zticks = n_elements(ztickv)
    zrange = mlt_range+[-1,1]*mlt_step*0.5
    zminor = ztick_step
    zticklen = de_ticklen*(cbpos[2]-cbpos[0])/(cbpos[3]-cbpos[1])*2
    ztitle = 'MLT (hr)'
    plot, zrange, [0,1], $
        xstyle=9, xticks=zticks, xtickv=ztickv, xminor=zminor, xticklen=zticklen, xtitle=ztitle, $
        ystyle=5, $
        position=cbpos, /nodata, /noerase
    sgcolorbar, index_colors, ct=ct, position=cbpos, /horizontal, $
        zrange=mlt_range, zticks=1, zminor=0, ztitle=' ', ztickformat='(A1)', zticklen=0


;---Panel settings.
    psym = 6
    symsize = 0.2
    xtitle = short_name+'!D0!N ('+de_unit+')'
    line_color = sgcolor('silver')
    fac_slopes = fltarr(nmlt_bin,ndim)

    for panel_id=0,ndim-1 do begin
        tpos = poss[*,panel_id]
        tpos[2] = tpos[0]+(tpos[3]-tpos[1])*fig_ysize/fig_xsize

        if panel_id ne 0 then ytitle = ' ' else ytitle = short_name+'!Ddot0!N ('+de_unit+')'
        if panel_id ne 0 then ytickformat='(A1)' else ytickformat=''

        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        fig_label = fig_letters[panel_id]+'. FAC '+short_name+'!D'+fac_labels[panel_id]+'!N '+de_labels[panel_id]
        xyouts, tx,ty,/normal, fig_label

        plot, de_range,de_range, /iso, $
            xstyle=1, xrange=de_range, xlog=de_log, xtitle=xtitle, xticklen=de_ticklen, xtickv=de_tickv, xticks=de_ticks, xminor=de_minor, $
            ystyle=1, yrange=de_range, ylog=de_log, ytitle=ytitle, yticklen=de_ticklen, ytickv=de_tickv, yticks=de_ticks, yminor=de_minor, $
            position=tpos, /noerase, ytickformat=ytickformat, /nodata
        plots, de_range,de_range, linestyle=1, color=line_color
        plots, [0,0],de_range, linestyle=1, color=line_color
        plots, de_range,[0,0], linestyle=1, color=line_color

        slopes = fltarr(nmlt_bin)
        foreach mlt_center, mlt_centers, mlt_id do begin
            slope = list()
            foreach probe, probes do begin
                if keyword_set(test) then if test eq 2 then if probe eq 'a' then continue
                if keyword_set(test) then if test eq 3 then if probe eq 'b' then continue
                prefix = 'rbsp'+probe+'_'
                xdata = (get_var_data(prefix+'de_fac'))[*,panel_id]
                ydata = (get_var_data(prefix+'dedot0_fac'))[*,panel_id]
                mlts = get_var_data(prefix+'mlt')

                index = lazy_where(mlts,'[)', mlt_vertices[mlt_id,*], count=count)
                if count eq 0 then continue
                xxs = xdata[index]
                yys = ydata[index]
                index = where(finite(xxs+yys))
                xxs = xxs[index]
                yys = yys[index]

                fit_result = linfit(xxs,yys)
                slope.add, fit_result[1]
                step = 5
                plots, xxs[0:*:step],yys[0:*:step], color=colors[mlt_id], psym=psyms[panel_id], symsize=symsize
            endforeach
            slope = mean(slope.toarray())
            slopes[mlt_id] = slope
            ;print, mlt_centers[mlt_id], slope
            ;oplot, de_range,de_range*slope, linestyle=1, color=colors[mlt_id]
        endforeach
        fac_slopes[*,panel_id] = slopes
    endfor

;---The coef.
    tpos = poss[*,npanel-1]
    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    fig_label = 'd. 3D E calibration per MLT bin'
    xyouts, tx,ty,/normal, fig_label

    tpos[1] = cbpos[3]+ychsz*0.2
    yrange = minmax(fac_slopes)
    yrange = [floor(yrange[0]),ceil(yrange[1])]
    ytitle = 'Avg E!Ddot0 !N/E!D0!N'
    ytick_step = 2
    ytickv = make_bins(yrange,ytick_step,/inner)
    yticks = n_elements(ytickv)-1
    yminor = 4

    xrange = zrange
    xtickv = ztickv
    xticks = zticks
    xminor = zminor
    xtitle = ' '

    plot, xrange,yrange, $
        xstyle=1, xrange=xrange, xlog=0, xtitle=xtitle, xticklen=de_ticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, yrange=yrange, ylog=0, ytitle=ytitle, yticklen=de_ticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
       position=tpos, /noerase, /nodata, xtickformat='(A1)'
    plots, xrange,[1,1], linestyle=1, color=line_color

    linestyles=[0,0,0]
    symsize = 0.5
    for ii=0,ndim-1 do begin
        yys = fac_slopes[*,ii]
        xxs = mlt_centers
        plots, xxs, yys, linestyle=linestyles[ii]
        for jj=0,nmlt_bin-1 do plots, xxs[jj],yys[jj], color=colors[jj], symsize=symsize, psym=psyms[ii]
    endfor

    tx0 = tpos[0]+xchsz*2
    ty0 = tpos[3]-ychsz*1
    label_size = 0.8
    for ii=0,ndim-1 do begin
        txs = tx0+[-1,1]*xchsz
        tys = ty0-ii*ychsz*0.9
        plots, txs,tys,/normal
        plots, tx0,tys[0],/normal, psym=psyms[ii], symsize=symsize
        tx = txs[1]+xchsz*1
        ty = tys-ychsz*0.2
        xyouts, tx,ty,/normal, fac_full_labels[ii], charsize=label_size
    endfor

    if keyword_set(test) then stop
    sgclose


end
