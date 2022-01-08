;+
; Need the position and time of each spacecraft. Use any of the 3 or 4 to determine the direction of the plane
;-


pro test_plane_direction

    test = 0

;---Load data.
    store_data, '*', /delete
    id = '2014_0828_10'
    df_paper_load_data, event=event, id=id
    event_info = event[id]
    time_range = event_info.time_range
    probes = event_info.probes
    nprobe = n_elements(probes)


;---Add settings.
    df_info = dictionary()
    foreach probe, probes do df_info[probe] = dictionary()
    df_info['thd'].time_range = time_double(['2014-08-28/10:12:15', '2014-08-28/10:13:54'])
    df_info['g15'].time_range = time_double(['2014-08-28/10:16:04', '2014-08-28/10:17:40'])
    df_info['the'].time_range = time_double(['2014-08-28/10:20:00', '2014-08-28/10:22:30'])
    df_info['tha'].time_range = time_double(['2014-08-28/10:20:39', '2014-08-28/10:24:21'])
    df_info['rbspb'].time_range = time_double(['2014-08-28/10:45:50', '2014-08-28/10:52:20'])
    df_info['g13'].time_range = time_double(['2014-08-28/10:49:00', '2014-08-28/10:55:40'])
    foreach probe, probes do begin
        df_info[probe].time = mean(df_info[probe].time_range)
        df_info[probe].rgsm = get_var_data(probe+'_r_gsm', at=df_info[probe].time)
        df_info[probe].rsm = get_var_data(probe+'_r_sm', at=df_info[probe].time)
        df_info[probe].mlt = get_var_data(probe+'_mlt', at=df_info[probe].time)
        df_info[probe].mlon = get_var_data(probe+'_mlon', at=df_info[probe].time)
        df_info[probe].mlat = get_var_data(probe+'_mlat', at=df_info[probe].time)
    endforeach


;---Order the variables.
    if ~keyword_set(sort_by) then sort_by = 'mlt'
    sort_data = fltarr(nprobe)
    foreach probe, probes, ii do sort_data[ii] = (df_info[probe])[sort_by]
    sorted_probes = probes[sort(sort_data)]
    if n_elements(sort) eq nprobe then sorted_probes = sort


;---Recipe:
;   0. Assume the velocity is V = (Vx,Vy,Vz).
;   1. Find the earliest probe, use it as the origin (0).
;   2. For the i-th probe, R_(i,0) dot V/|V| = |V| T_(i,0), where
;       R_(i,0) = R_i - R_0, T_(i,0) = T_i - T_0
;   3. We can then fit R_(i,0) and T_(i,0), to get V/|V|^2.


;---Plot.
    fig_xsize = 4
    fig_ysize = 4
    lmargin = 8
    rmargin = 6
    tmargin = 2
    bmargin = 5
    psym = 6
    ticklen = -0.01
    full_ysz = 0.8
    half_ysz = 0.3
    label_size = 0.8
    lineskip = 0.3

    color_start = 50
    color_end = 250
    color_table = 40
    colors = round(smkarthm(color_start,color_end, nprobe, 'n'))
    for ii=0, nprobe-1 do colors[ii] = sgcolor(colors[ii], ct=color_table)

    sgopen, 0, xsize=fig_xsize, ysize=fig_xsize
    tpos = sgcalcpos(1,lmargin=lmargin,rmargin=rmargin,tmargin=tmargin,bmargin=bmargin, xchsz=xchsz, ychsz=ychsz)
    sgclose, /wdelete

    pos_panel = dictionary()
    xsm = fltarr(nprobe)
    ysm = fltarr(nprobe)
    foreach probe, sorted_probes, ii do begin
        rsm = df_info[probe].rsm
        xsm[ii] = rsm[0]
        ysm[ii] = rsm[1]
    endforeach
    xrange = minmax([xsm,-1,1])
    pos_panel.xrange = [ceil(xrange[1]),floor(xrange[0])]
    pos_panel.xrange = [1,-11]
    pos_panel.xtitle = 'SM X (Re)'
    pos_panel.xminor = 5
    pos_panel.xtickv = make_bins(pos_panel.xrange, pos_panel.xminor, /inner)
    pos_panel.xticks = n_elements(pos_panel.xtickv)-1
    yrange = minmax([ysm,-1,1])
    pos_panel.yrange = [ceil(yrange[1]),floor(yrange[0])]
    pos_panel.ytitle = 'SM Y (Re)'
    pos_panel.yminor = 5
    pos_panel.ytickv = make_bins(pos_panel.yrange, pos_panel.yminor, /inner)
    pos_panel.yticks = n_elements(pos_panel.ytickv)-1
    pos_panel.aspect_ratio = abs(double(pos_panel.yrange[1]-pos_panel.yrange[0])/(pos_panel.xrange[1]-pos_panel.xrange[0]))
    pos_panel.xticklen = ticklen/pos_panel.aspect_ratio
    pos_panel.yticklen = ticklen

    pos_panel.xsize = tpos[2]-tpos[0]
    pos_panel.ysize = pos_panel.xsize*pos_panel.aspect_ratio
    fig_ysize = (tmargin*ychsz+pos_panel.ysize+bmargin*ychsz)*fig_xsize

    if keyword_set(test) then begin
        file = test+1
        magnify = 2
    endif else begin
        file = join_path([srootdir(),'test_plane_direction.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify

    tpos = sgcalcpos(1,lmargin=lmargin,rmargin=rmargin,tmargin=tmargin,bmargin=bmargin)
    plot, pos_panel.xrange, pos_panel.yrange, $
        xstyle=5, ystyle=5, $
        /nodata, /noerase, position=tpos, /isotropic, _extra=pos_panel.tostruct()
    foreach probe, sorted_probes, ii do begin
        plots, xsm[ii], ysm[ii], color=colors[ii], psym=psym, symsize=label_size

        ; Add spacecraft short name.
        mission = resolve_probe(probe)
        label = strupcase(mission.short_name+mission.probe)

        tmp = convert_coord(xsm[ii], ysm[ii], /data, /to_normal)
        tx = tmp[0]+xchsz*1*label_size
        ty = tmp[1]
        xyouts, tx,ty,/normal, label, color=colors[ii], charsize=label_size
    endforeach

    ; Add earth and lines.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45
    plots, xs, ys
    foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

    plots, pos_panel.xrange, [0,0], linestyle=1
    plots, [0,0], pos_panel.yrange, linestyle=1

    plot, pos_panel.xrange, pos_panel.yrange, $
        xstyle=1, ystyle=1, /isotropic, $
        /nodata, /noerase, position=tpos, _extra=pos_panel.tostruct()

    r_mag0 = 1.
    v_mag0 = 20  ; scale vmag to a proper length on the plot.
    scale = r_mag0/v_mag0
    hsize = (size(file,/type) eq 7)? 80:10

;---3 probes.
    test_list = list()
    test_list.add, ['thd','g15','the']
    test_list.add, ['tha','rbspb','g13']
    test_list.add, ['g15','thd','tha']
    test_list.add, ['rbspb','the','tha']
    foreach test_probes, test_list do begin
        ;test_probes = ['thd','g15','the']
        ;test_probes = ['tha','rbspb','g13']
        ntest_probe = n_elements(test_probes)
        times = dblarr(ntest_probe)
        foreach probe, test_probes, ii do times[ii] = df_info[probe].time
        index = sort(times)
    
        dim_indexs = [0,1]
        ndim_index = n_elements(dim_indexs)
    
        the_probes = test_probes[index]
        test_info = df_info[the_probes]
        rr = dblarr(ndim_index,ntest_probe-1)
        rr_0 = (test_info[the_probes[0]].rsm)[dim_indexs]
        for ii=1, ntest_probe-1 do begin
            rr_i = (test_info[the_probes[ii]].rsm)[dim_indexs]
            rr[*,ii-1] = rr_i-rr_0
        endfor
        tt = dblarr(ntest_probe-1)
        tt_0 = test_info[the_probes[0]].time
        for ii=1, ntest_probe-1 do begin
            tt_i = test_info[the_probes[ii]].time
            tt[ii-1] = tt_i-tt_0
        endfor
        vv = la_linear_equation(rr,tt)
        v_hat = sunitvec(vv)
    
        rr_normal = dblarr(ntest_probe-1)
        for ii=0, ntest_probe-2 do rr_normal[ii] = sdot(rr[*,ii],v_hat)
    
        ;plot, tt, rr_normal, /ynozero, psym=-1
        fit_result = linfit(tt, rr_normal)
        v_mag = fit_result[1]*6378d
        
    
        x0 = 0
        foreach probe, test_probes do x0 += test_info[probe].rsm[0]
        x0 = x0/ntest_probe
        y0 = 0
        foreach probe, test_probes do y0 += test_info[probe].rsm[1]
        y0 = y0/ntest_probe
        x1 = x0+v_hat[0]*v_mag*scale
        y1 = y0+v_hat[1]*v_mag*scale
        
        arrow, x0,y0,x1,y1, /data, /solid, hsize=hsize
        tmp = convert_coord([x0,y0], /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]
        xyouts, tx+xchsz*0.5,ty-ychsz*half_ysz,/normal, $
            strjoin(strupcase(test_probes),','), charsize=label_size*0.6
    endforeach
    
    txs = pos_panel.xrange[1]+0.5+[0,r_mag0]
    tys = [0,0]
    for ii=0,1 do begin
        tmp = convert_coord([txs[ii],tys[ii]], /data, /to_normal)
        txs[ii] = tmp[0]
        tys[ii] = tmp[1]
    endfor
    scale = txs[1]-txs[0]
    
    ;txs = tpos[2]-xchsz*1-[0,scale]
    tys = tpos[3]-ychsz*0.5
    plots, txs,tys,/normal
    tx = mean(txs)
    ty = tys[0]-ychsz*1
    xyouts, tx,ty,/normal, sgnum2str(v_mag0)+' km/s', alignment=0.5, charsize=label_size


    if keyword_set(test) then stop
    sgclose
end

test_plane_direction
end
