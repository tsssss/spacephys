;+
; Load Dst/AE, EWOgram for DP and current.
;-


test = 0


    long_time_range = time_double(['2014-08-28/08:00','2014-08-28/24:00'])

    if n_elements(project) eq 0 then project = azim_df_load_project()
    azim_df_load_basic_data, project=project
    _2014_0828_10_load_data
    mlt_range = [-3,9]
    probe_xys = [-1,1]
    asi_time_range = time_double(['2014-08-28/10:05','2014-08-28/10:20'])
    asi_mlt_range = [-1.5,2]

    line_linestyle = 2
    full_ychsz = 0.7
    half_ychsz = 0.35
    label_xshift = 1
    label_yshift = full_ychsz
    label_size = 0.7
    constant_linestyle = 1
    sorted_probes = ['thd','g15','the','tha','g13','rbspb']
    rxy_range = [4.,30]
    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'case2_fig_long_term.pdf'])

    ; The lines.
    line_color = sgcolor('grey')
    line_xx = list()
    line_yy = list()
    line_xx.add, time_double('2014-08-28/'+['10:05','11:10'])
    line_yy.add, [0.0,total(line_xx[0]*[-1,1])/60*2.1/15]
;    line_xx.add, time_double('2014-08-28/'+['16:30','17:10'])
;    line_yy.add, [0.0,9.0]
;    line_xx.add, time_double('2014-08-28/'+['15:40','16:30'])
;    line_yy.add, [0.0,9.0]


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    panel_ysize = 1.2    ; inch.
    ypans = [0.6,1,1]
    nypanel = n_elements(ypans)
    ypad = 0.4
    panel_x_ratio = 0.5 ; inch per hour.
    panel_xsize = abs(total(long_time_range*[-1,1]))/3600*panel_x_ratio
    margins = [10.,4,8,1]
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_ysize = panel_ysize*total(ypans)+(total(margins[[1,3]])+ypad*(nypanel-1))*abs_ychsz


;---Load data.
    if check_if_update('ae',long_time_range) then omni_read_index, long_time_range


;---Common x-axis.
    xrange = long_time_range
    xticklen_chsz = -0.2
    yticklen_chsz = -0.4
    xminor = 6 ; hr.
    xstep = 3600.
    xtickv = make_bins(xrange, xstep, /inner) ; make time line up.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    secofday = 86400d
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


    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, /inch, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypanel, ypans=ypans, ypad=ypad, margins=margins)


    ; Label.
    fig_labels = letters(nypanel)+'. '+['AE','Dipolari-!C     zation','Upward!C    current']
    for ii=0,nypanel-1 do begin
        tpos = poss[*,ii]
        tx = label_xshift*xchsz
        ty = tpos[3]-label_yshift*ychsz
        xyouts, tx,ty,/normal, fig_labels[ii]
    endfor



;---AE.
    the_var = 'ae'
    tpos = poss[*,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    xtickformat = '(A1)'
    ystep = 500.
    constants = [500]
    yys = get_var_data(the_var, in=long_time_range, times=xxs)
    yrange = minmax([yys,constants])
    ytickv = make_bins(yrange, ystep)
    yticks = n_elements(ytickv)-1
    yrange = minmax(ytickv)
    yminor = 2
    ytitle = '(nT)'

    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, /nodata, /noerase

    ; Add data.
    oplot, xxs, yys

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = [500,1000]
    x_constants = make_bins(long_time_range, 3600)
    x_constants = x_constants[1:-2]
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add axes.
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


;---Panel c, EWOgram of current.
    tpos = poss[*,2]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz

    j_ewo_var = 'thg_j_up_ewo'
    j_ewo_var1 = j_ewo_var+'1'
    copy_data, j_ewo_var, j_ewo_var1
    get_data, j_ewo_var1, xx,yy,zz
    yy *= 1e-3
    store_data, j_ewo_var1, xx,yy,zz
    options, j_ewo_var1, 'ztitle', 'Upward current (kA)'
    options, j_ewo_var1, 'zrange', [0,50]
    options, j_ewo_var1, 'zcharsize', 0.8
    options, j_ewo_var1, 'zposition', cbpos
    j_ewo_color_table = 62
    options, j_ewo_var1, 'xticklen', xticklen
    options, j_ewo_var1, 'yticklen', yticklen
    options, j_ewo_var1, 'color_table', j_ewo_color_table
    tplot, j_ewo_var1, trange=long_time_range, position=tpos, /noerase

    yrange = mlt_range
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, noerase=1, nodata=1

    ; Add grid.
    grid_color = sgcolor('gray')
    grid_linestyle = 1
    y_constants = [0,3,6]
    x_constants = make_bins(long_time_range, 3600)
    x_constants = x_constants[1:-2]
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add lines.
    for ii=0, line_xx.length-1 do begin
        oplot, line_xx[ii], line_yy[ii], color=line_color, linestyle=line_linestyle
    endfor


    ; Add zoom-in box.
;    oplot, asi_time_range, -2+[0,0], color=grid_color
;    oplot, time_double('2014-08-28/10:00')+[0,0], asi_mlt_range, color=grid_color
    xxs = asi_time_range
    yys = asi_mlt_range
    plots, xxs[[0,1,1,0,0]], yys[[0,0,1,1,0]], color=grid_color

    ; Add box.
    xtickformat = '(A1)'
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase


;---Panel b. UT-MLT long.
    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    xtickformat = '(A1)'
    ytitle = 'MLT (hr)'
    yrange = mlt_range
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    constant = 0

    data_var_suffix = 'scaled_theta'
    ztitle = 'Scaled '+tex2str('theta')+' (deg)'
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

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=5, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, ytitle=ytitle, $
        /nodata, /noerase

    pos_var_suffix = 'pseudo_mlt'
    flags = list()
    foreach probe, reverse(sorted_probes) do begin
        prefix = probe+'_'
        pos_var = prefix+pos_var_suffix
        data_var = prefix+data_var_suffix
        zzs = get_var_data(data_var, in=long_time_range, times=xxs)
        yys = get_var_data(pos_var, at=xxs)
        mlt = get_var_data(prefix+'pseudo_mlt', at=xxs)
        zzs = azim_df_normalize_theta(zzs, zrange=tilt_range, ct=ct, /reverse_ct)

        ; Remove data out of range.
        r_sm = get_var_data(prefix+'r_sm', at=xxs)
        rxy = snorm(r_sm[*,0:1])
        ; Exclude data outside the MLT range.
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data ouside the distance range.
        index = where_pro(rxy, '][', rxy_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        ; Exclude data outside magnetopause.
        magn_flags = check_if_in_magn(r_sm)
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
        for ii=0, nxx-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    ; Add grid.
    foreach constant, x_constants do $
        plots, constant+[0,0], yrange, color=grid_color, linestyle=grid_linestyle
    foreach constant, y_constants do $
        plots, xrange, constant+[0,0], color=grid_color, linestyle=grid_linestyle

    ; Add lines.
    for ii=0, line_xx.length-1 do begin
        oplot, line_xx[ii], line_yy[ii], color=line_color, linestyle=line_linestyle
    endfor
    ; Add angular speed.
    tx = time_double('2014-08-28/10:55')
    ty = 7.7
    xyouts, tx,ty,/data, '2.1 deg/min', alignment=1, charsize=0.8;, color=line_color

    ; Add box.
    xtickformat = '(A1)'
    plot, xrange, yrange, position=tpos, $
        xstyle=1, xrange=xrange, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, xtickformat=xtickformat, $
        ystyle=1, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, ytitle=ytitle, $
        /nodata, /noerase

    ; Draw color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks


;---EWOgram for aurora.
    rough_pos = poss[*,2]
    rough_pos[0] = (rough_pos[2]-rough_pos[0])*5/8+rough_pos[0]
    poss = sgcalcpos(1,2, region=rough_pos, xpad=1, margins=[4,1.8,6,0.5])

    line_xs = list()
    line_ys = list()
    line_xs.add, time_double('2014-08-28/'+['10:12:10','10:12:48'])
    line_ys.add, [0.12,-0.54]
    line_xs.add, time_double('2014-08-28/'+['10:12:10','10:12:48'])
    line_ys.add, [0.12,0.99]


    tpos = poss[*,1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])*0.5
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])*0.5

    mlat_range = [60,70]
    xcharsize = 0.8

    asi_var = 'thg_mlonimg'
    get_data, asi_var, times, mlon_image, limits=lim
    mlat_bins = lim.mlat_bins
    mlon_bins = lim.mlon_bins

    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nmlat_bin)
    ntime = n_elements(times)
    nmlon_bin = n_elements(mlon_bins)
    ewo = fltarr(ntime,nmlon_bin)
    foreach time, times, ii do begin
        foreach mlon, mlon_bins, jj do begin
            tmp = reform(mlon_image[ii,jj,mlat_index])
            ewo[ii,jj] = mean(tmp)
        endforeach
    endforeach
    ewo_var = 'thg_mlonimg_ewo'
    store_data, ewo_var, times, ewo, mlon_bins, limits={$
        spec: 1, $
        color_table: 49, $
        no_interp: 1, $
        ytitle: 'MLon (deg)', $
        ystyle: 1, $
        ztitle: 'Photon count (#)', $
        zlog: 0 , $
        yticklen: yticklen, $
        xticklen: xticklen }

    ; Convert to mlt_ewo.
    mlt_binsize = total(mlon_bins[0:1]*[-1,1])/15
    mlt_bins = make_bins(asi_mlt_range,mlt_binsize)
    nmlt_bin = n_elements(mlt_bins)

    mlt_ewo_var = 'thg_mlonimg_mlt_ewo'
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
    ystep = 1.0
    ytickv = make_bins(asi_mlt_range, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 2
    store_data, mlt_ewo_var, times, mlt_ewo, mlt_bins, limits={$
        spec: 1, $
        no_interp: 1, $
        color_table: 40, $
        ytitle: ' ', $
        ystyle: 1, $
        yrange: asi_mlt_range, $
        ytickv: ytickv, $
        yticks: yticks, $
        yminor: yminor, $
        ytickformat: '(A1)', $
        ztitle: 'Photon #', $
        zcharsize: 0.8, $
        zlog: 0 , $
        zrange: [0,160], $
        yticklen: yticklen, $
        xticklen: xticklen }


    tplot, mlt_ewo_var, trange=asi_time_range, position=tpos, noerase=1, nouttick=1

    xrange = asi_time_range
    xminor = 5
    xstep = xminor*60.
    xtickv = make_bins(xrange, xstep)
    xticks = n_elements(xtickv)-1
    xtickn = time_string(xtickv,tformat='hhmm')
    yrange = asi_mlt_range
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xtickname=xtickn, xticklen=xticklen, charsize=xcharsize, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yticklen=yticklen, ytickformat='(A1)', $
        position=tpos, /nodata, /noerase
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.9
    xyouts, tx,ty,/normal, 'e. Aurora', color=sgcolor('white')

    for ii=0, line_xs.length-1 do begin
        oplot, line_xs[ii], line_ys[ii], color=sgcolor('white')
    endfor

    ; Add lines.
    for ii=0, line_xx.length-1 do begin
        oplot, line_xx[ii], line_yy[ii], color=line_color, linestyle=line_linestyle
    endfor



    ; Zoom in of current.
    tpos = poss[*,0]
    j_ewo_var2 = j_ewo_var+'_tmp'
    copy_data, j_ewo_var, j_ewo_var2
    store_data, j_ewo_var2, limits={$
        spec: 1, $
        color_table: j_ewo_color_table, $
        no_color_scale: 1, $
        no_interp: 1, $
        ystyle: 1, $
        ytitle: 'MLT (hr)', $
        yrange: asi_mlt_range, $
        ytickv: ytickv, $
        yticks: yticks, $
        yminor: yminor, $
        yticklen: yticklen, $
        ytickformat: '(A1)', $
        xticklen: xticklen }
    tplot, j_ewo_var2, trange=asi_time_range, position=tpos, noerase=1, nouttick=1

    xtickn[-1] = ' '
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xtickname=xtickn, xticklen=xticklen, charsize=xcharsize, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yticklen=yticklen, $
        position=tpos, /nodata, /noerase
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.9
    xyouts, tx,ty,/normal, 'd. Zoom-in'

    for ii=0, line_xs.length-1 do begin
        oplot, line_xs[ii], line_ys[ii], color=sgcolor('white')
    endfor

    ; Add lines.
    for ii=0, line_xx.length-1 do begin
        oplot, line_xx[ii], line_yy[ii], color=line_color, linestyle=line_linestyle
    endfor

    if keyword_set(test) then stop
    sgclose
end
