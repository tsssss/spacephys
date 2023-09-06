;+
; Plot theta and DF for given time range.
;-

pro check_df_in_time_range, time_range, probes=probes


    mlt_range = [-1,1]*9
    rxy_range = [4.,30]


    min_scaled_height = 8

    all_probes = ['rbsp'+letters('b'),'th'+letters('e'),'g'+['13','14','15'],'mms1']
    if n_elements(probes) eq 0 then probes = all_probes
    df_list = list()
    foreach probe, probes do begin
        dfs = azim_df_read_df(time_range, probe=probe, mlt_range=mlt_range, rxy_range=rxy_range)
        if dfs.length eq 0 then continue
        df_list.add, dfs, /extract
    endforeach

    probe_list = list()
    foreach df, df_list do if n_elements(probe_list.where(df.probe)) eq 0 then probe_list.add, df.probe

    project = azim_df_load_project()
    foreach probe, probe_list do begin
        prefix = probe+'_'
        if check_if_update(prefix+'r_sm', time_range) then azim_df_read_data, 'r_sm', time_range=time_range, probe=probe, project=project
        if check_if_update(prefix+'theta', time_range) then azim_df_read_data, 'theta', time_range=time_range, probe=probe, project=project
        if check_if_update(prefix+'pseudo_mlt', time_range) then begin
            r_sm = get_var_data(prefix+'r_sm', times=times)
            mlt = azim_df_calc_pseudo_mlt(r_sm)
            store_data, prefix+'pseudo_mlt', times, mlt
        endif
    endforeach

    if df_list.length eq 0 then stop
    if probe_list.length eq 0 then stop


;---Plot settings.
    xtickformat = ''
    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.
    full_ychsz = constant('full_ychsz')
    half_ychsz = full_ychsz*0.5
    lineskip = constant('lineskip')
    label_size = 0.7
    label_xshift = 10
    label_yshift = full_ychsz
    secofday = constant('secofday')
    bar_thick = keyword_set(test)? 0.5: 4

    rad = constant('rad')
    deg = constant('deg')


    sgopen, 0, xsize=6, ysize=4, /inch
    pos_list = list()
    pos_list.add, sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, rmargin=10, lmargin=12)

    tpos = pos_list[0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ;---The common x-axis setting.
        xrange = time_range
        xstep = constant('secofhour')
        if total(time_range*[-1,1]) le xstep then xstep = 600.
        xticks = floor(total(xrange*[-1,1])/xstep)
        xtickv = smkarthm(xrange[0], xstep, xticks+1, 'x0')
        xminor = 6
        xtickn = strarr(xticks+1)
        for ii=0, xticks do begin
            the_time = xtickv[ii]
            xtickn[ii] = time_string(the_time,tformat='hh:mm')
            date = time_string(the_time,tformat='YYYY-MM-DD')
            if ii eq 0 then begin
                xtickn[ii] += '!C'+date
                continue
            endif
            if the_time mod secofday ne 0 then continue
            xtickn[ii] += '!C'+date
        endfor

    ;---Common z-range, colorbar, symbol.
        theta_var = 'theta'
        ztitle = 'Log!D10!N[Detrended tilt angle (deg)]'
        theta_range = [-1,1]*64
        zrange = theta_range
        ztickv = [-64,-16,-4,0,4,16,64]
        ztickn = string(ztickv,format='(I0)')
        ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
        zrange = alog10(abs(zrange))*[-1,1]
        zticks = n_elements(ztickv)-1
        spec_psym = 8
        spec_symsize = 0.3
        usersym, [1,1,-1,-1,1]*0.5, [-1,1,1,-1,-1]*1, /fill
        spec_ct = 70


    ;---Common x data.
        time_step = 10.
        xxs = make_bins(time_range, time_step)
        nxx = n_elements(xxs)

    ;---MLT-UT.
        pos_var = 'pseudo_mlt'
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        ytitle = 'MLT (hr)'
        yminor = 3
        yrange = []
        foreach probe, probe_list do yrange = [yrange,minmax(get_var_data(probe+'_pseudo_mlt'))]
        yrange = [floor(min(yrange)),ceil(max(yrange))]
        ytickv = make_bins(yrange,yminor, /inner)
        yticks = n_elements(ytickv)-1
        constants = 0


        ; Set up coord.
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, /nodata, /noerase

        ; Plot data.
        foreach probe, probe_list do begin
            prefix = probe+'_'
            yys = get_var_data(prefix+pos_var, at=xxs)
            zzs = get_var_data(prefix+theta_var, at=xxs)
            mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
            zzs = azim_df_scale_theta(zzs, mlt) ; apply MLT scaling.
            zzs = azim_df_normalize_theta(zzs, zrange=theta_range, ct=spec_ct, /reverse_ct)

            ; Remove data outside ROI.
            index = where_pro(mlt, '][', mlt_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            rsm = get_var_data(prefix+'r_sm', at=xxs)
            rxy = snorm(rsm[*,0:1])
            index = where_pro(rxy, '][', rxy_range, count=count)
            if count ne 0 then yys[index] = !values.f_nan
            index = where(finite(zzs,/nan), count)
            if count ne 0 then yys[index] = !values.f_nan

            ; Plot data.
            for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=spec_psym, symsize=spec_symsize
        endforeach

        ; Add DF.
        foreach df, df_list do begin
            ty = df.obs_mlt
            tx = df.obs_time
            plots, tx,ty,/data, psym=1, symsize=0.5
        endforeach

        ; Add large DF.
        foreach df, df_list do begin
            if round(df.scaled_height) le min_scaled_height then continue
            ty = df.obs_mlt
            tx = df.obs_time
            plots, tx,ty,/data, psym=6, symsize=0.5, color=sgcolor('red')
        endforeach

        ; Add notations.
        foreach constant, constants do oplot, xrange, constant+[0,0], linestyle=1

        ; Add axes.
        plot, xrange, yrange, position=tpos, $
            xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
            ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
            /nodata, /noerase

        ; Add label.
        tx = tpos[0]-label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, 'c. UT-MLT'

        cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
        sgcolorbar, 256-findgen(256), ct=spec_ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


    foreach df, df_list do azim_df_vertex_write, df

    stop
end


time_range = time_double(['2014-08-28/09:30','2014-08-28/11:30'])
;time_range = time_double(['2016-10-13/11:30','2016-10-13/14:00'])
;time_range = time_double(['2007-11-20/16:50','2007-11-20/17:30'])
time_range = time_double(['2008-02-28/03:00','2008-02-28/04:30'])
time_range = time_double(['2008-02-29/08:00','2008-02-29/09:00'])

;time_range = time_double(['2013-08-31/02:00','2013-08-31/03:00'])
;time_range = time_double(['2013-10-09/06:00','2013-10-09/08:00'])
;time_range = time_double(['2013-10-14/02:30','2013-10-14/04:00'])


;time_range = time_double(['2009-03-19/08:00','2009-03-19/08:50'])
;time_range = time_double(['2009-03-15/08:40','2009-03-15/09:10'])
;time_range = time_double(['2009-02-27/07:40','2009-02-27/08:10'])
;time_range = time_double(['2009-03-05/03:00','2009-03-05/03:30'])
;time_range = time_double(['2009-02-27/07:40','2009-02-27/08:10'])
check_df_in_time_range, time_range
end
