;+
; Check the distribution of all DFs in the x-y plane, for
;   1. before
;   2. after physical contraints are applied.
;-

function bin_df_on_xy_plane, df_list, xbins=xbins, ybins=ybins

    ; Prepare bins.
    xrange = [-1,1]*9
    xstep = 0.5
    xbins = make_bins(xrange,xstep)
    nxbin = n_elements(xbins)-1
    yrange = [4,30]
    ystep = 2
    ybins = make_bins(yrange,ystep,/inner)
    nybin = n_elements(ybins)-1

    ; Extract related info.
    ndf = df_list.length
    mlts = fltarr(ndf)
    rxys = fltarr(ndf)
    rxxs = fltarr(ndf)
    ryys = fltarr(ndf)
    foreach df, df_list, ii do begin
        mlts[ii] = df.obs_mlt
        rxys[ii] = df.obs_rxy
        rxxs[ii] = df.obs_r_sm[0]
        ryys[ii] = df.obs_r_sm[1]
    endforeach

    ; Bin data.
    counts = fltarr(nxbin,nybin)+!values.f_nan
    xxs = mlts
    yys = rxys
    for ii=0,nxbin-1 do begin
        for jj=0,nybin-1 do begin
            index = where(xxs ge xbins[ii] and xxs lt xbins[ii+1] and yys ge ybins[jj] and yys lt ybins[jj+1], count)
            if count ne 0 then counts[ii,jj] = count
        endfor
    endfor

    return, counts

end



pro fig_check_midnight_df_distribution

    project = azim_df_load_project()
    data_dir = project.data_dir
    regions = ['around','pre','post']+'_midn'
    bin_list = list()

test = 0

;---Setp 1, all DF in ROI.

    ; Load events.
    files = list()
    foreach region, regions do begin
        files.add, join_path([data_dir,'azim_df_search_'+region+'_search_vertex.txt'])
    endforeach
    events = list()
    foreach file, files do begin
        events.add, azim_df_subgroup_read_file(file), /extract
    endforeach

    ; Extract unique DFs.
    df_list = list()
    foreach event, events do begin
        df_list.add, event.df_list, /extract
    endforeach
    ndf = df_list.length
    df_str_list = strarr(ndf)
    foreach df, df_list, ii do df_str_list[ii] = df.probe+'_'+time_string(df.obs_time,tformat='YYYY_MMDD_hhmm_ss')
    index = uniq(df_str_list, sort(df_str_list))
    uniq_df_list = df_list[index]
    bin_list.add, bin_df_on_xy_plane(uniq_df_list, xbins=xbins, ybins=ybins)


;---Step 2, the 10% filter.

    ; Load events.
    files = list()
    foreach region, regions do begin
        files.add, join_path([data_dir,'azim_df_search_'+region+'_search_subgroup.txt'])
    endforeach
    events = list()
    foreach file, files do begin
        events.add, azim_df_subgroup_read_file(file), /extract
    endforeach


    ; Extract unique DFs.
    df_list = list()
    foreach event, events do begin
        df_list.add, event.df_list, /extract
    endforeach
    ndf = df_list.length
    df_str_list = strarr(ndf)
    foreach df, df_list, ii do df_str_list[ii] = df.probe+'_'+time_string(df.obs_time,tformat='YYYY_MMDD_hhmm_ss')
    index = uniq(df_str_list, sort(df_str_list))
    uniq_df_list = df_list[index]
    bin_list.add, bin_df_on_xy_plane(uniq_df_list)



;---Step 3, require timing can be robustly done.

    ; Load evetns.
    events = azim_df_find_subgroup(project=project)

    ; Extract unique DFs.
    df_list = list()
    foreach event, events do begin
        df_list.add, event.df_list, /extract
    endforeach
    ndf = df_list.length
    df_str_list = strarr(ndf)
    foreach df, df_list, ii do df_str_list[ii] = df.probe+'_'+time_string(df.obs_time,tformat='YYYY_MMDD_hhmm_ss')
    index = uniq(df_str_list, sort(df_str_list))
    uniq_df_list = df_list[index]
    bin_list.add, bin_df_on_xy_plane(uniq_df_list)



;---Step 4, requires consistent propagation.

    ; Load evetns.
    events = azim_df_find_dfgroup(project=project)

    ; Extract unique DFs.
    df_list = list()
    foreach event, events do begin
        df_list.add, event.df_list, /extract
    endforeach
    ndf = df_list.length
    df_str_list = strarr(ndf)
    foreach df, df_list, ii do df_str_list[ii] = df.probe+'_'+time_string(df.obs_time,tformat='YYYY_MMDD_hhmm_ss')
    index = uniq(df_str_list, sort(df_str_list))
    uniq_df_list = df_list[index]
    bin_list.add, bin_df_on_xy_plane(uniq_df_list)



;---Make a plot.
    plot_dir = project.plot_dir
    plot_file = join_path([plot_dir,'fig_check_midnight_df_distribution.pdf'])
    if keyword_set(test) then plot_file = 0
    margins = [8,5,8,2]
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    xpan_size = 2
    ypan_size = xpan_size
    nxpanel = bin_list.length
    nypanel = 1
    xpad = [8,intarr(nxpanel)+2]
    ypad = 0.2
    xfig_size = xpan_size*nxpanel+total(xpad)*abs_xchsz+total(margins[[0,2]])*abs_xchsz
    yfig_size = ypan_size*nypanel+total(ypad)*abs_ychsz+total(margins[[1,3]])*abs_ychsz
    sgopen, plot_file, xsize=xfig_size, ysize=yfig_size, /inch
    ;sgopen, plot_file, xsize=6, ysize=4, /inch
    poss = sgcalcpos(nypanel,nxpanel, xchsz=xchsz, ychsz=ychsz, xpad=xpad, margins=margins)




    xrange = [10,-30]
    yrange = [1,-1]*20
    rxy_range = [4,30]
    mlt_range = [-1,1]*9
    xticklen = -0.02
    yticklen = -0.02
    fig_labels = letters(nxpanel)+'. Step '+['1','2/1','3/2','4/3']


;---Plot the abs # of DF.
    tpos = poss[*,0]
    cbpos = tpos[[2,1,2,3]]+xchsz*[1,0,2,0]
    top_color = 200
    bottom_color = 50
    count_ct = 55
    data_colors = smkarthm(bottom_color, top_color, 1, 'dx')
    ztitle = 'Log10 # of dipolarizations'
    data_range = [0.6,6e2]
    data_range = [1,3]
    sgcolorbar, data_colors, horizontal=0, position=cbpos, $
        zrange=data_range, ztitle=ztitle, ct=count_ct

    xpanel_id = 0

    xtitle = 'SM X (Re)'
    ytitle = (xpanel_id eq 0)? 'SM Y (Re)': ''
    ytickformat = (xpanel_id eq 0)? '': '(A1)'
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
        /nodata, /noerase, position=tpos, /iso

    data = alog10(bin_list[xpanel_id])
    zzs = bytscl(data, min=data_range[0], max=data_range[1], top=top_color-bottom_color)+bottom_color

    nxbin = n_elements(xbins)-1
    nybin = n_elements(ybins)-1
    rrs = ybins
    tts = (180+xbins)*15*constant('rad')
    for xid=0,nxbin-1 do begin
        for yid=0,nybin-1 do begin
            if finite(data[xid,yid],/nan) then continue
            the_tt = tts[xid:xid+1]
            the_rr = rrs[yid:yid+1]
            txxs = the_rr[[0,1,1,0,0]]*cos(the_tt[[0,0,1,1,0]])
            tyys = the_rr[[0,1,1,0,0]]*sin(the_tt[[0,0,1,1,0]])
            if product(minmax(txxs)-xrange[0]) lt 0 then continue
            if product(minmax(tyys)-yrange[0]) lt 0 then continue
            if product(minmax(txxs)-xrange[1]) lt 0 then continue
            if product(minmax(tyys)-yrange[1]) lt 0 then continue
            the_color = sgcolor(zzs[xid,yid], ct=count_ct)
            polyfill, txxs, tyys, color=the_color
        endfor
    endfor

    ; Circle for earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)

    ; Add magnetopause.
    magn_test = fltarr(nangle,3)
    magn_test[*,0] = circle_x*30
    magn_test[*,1] = circle_y*30
    pdyn = 10
    tmp = check_if_in_magn(magn_test, magn_pos=magn_pos, dynamic_pressure=pdyn)
    magn_pos = magn_pos[*,[0,2,1]]
    magn_pos = [magn_pos,magn_pos]
    magn_pos[nangle:nangle*2-1,1] *= -1


    ; Add more circles.
    foreach tmp, [5,10,20] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
    foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp

    ; Add mlt_range.
    foreach tmp, mlt_range do begin
        tt = tmp*15*constant('rad')+!dpi
        oplot, rxy_range*cos(tt), rxy_range*sin(tt)
    endforeach


    ; Add earth.
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y


    ; Add magnetopause.
    oplot, magn_pos[*,0], magn_pos[*,1]

    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    xyouts, tx,ty,/normal, fig_labels[xpanel_id]



;---Plot the ratios.
    data_range = [0.0,0.5]*100
    ;cbpos = [poss[0,0],poss[3,0],poss[2,-1],poss[3,0]]+ychsz*[0,0.5,0,1]
    cbpos = poss[[2,1,2,3],-1]+xchsz*[1,0,2,0]
    top_color = 200
    bottom_color = 50
    count_ct = 49
    data_colors = smkarthm(bottom_color, top_color, 1, 'dx')
    zero_color = sgcolor(bottom_color, ct=53)
    ztitle = 'Selection rate per bin (%)'
    sgcolorbar, data_colors, horizontal=0, position=cbpos, $
        zrange=data_range, ztitle=ztitle, ct=count_ct



    for xpanel_id=1,nxpanel-1 do begin
        tpos = poss[*,xpanel_id]
        xtitle = 'SM X (Re)'
        ytitle = (xpanel_id eq 0)? 'SM Y (Re)': ''
        ytickformat = (xpanel_id eq 0)? '': '(A1)'
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
            ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, ytickformat=ytickformat, $
            /nodata, /noerase, position=tpos, /iso

        data = bin_list[xpanel_id]/bin_list[xpanel_id-1]*100
        index = where(finite(bin_list[0]) and finite(data,/nan), count)
        if count ne 0 then data[index] = 0
        zzs = bytscl(data, min=data_range[0], max=data_range[1], top=top_color-bottom_color)+bottom_color

        nxbin = n_elements(xbins)-1
        nybin = n_elements(ybins)-1
        rrs = ybins
        tts = (180+xbins)*15*constant('rad')
        for xid=0,nxbin-1 do begin
            for yid=0,nybin-1 do begin
                if finite(data[xid,yid],/nan) then continue
                the_tt = tts[xid:xid+1]
                the_rr = rrs[yid:yid+1]
                txxs = the_rr[[0,1,1,0,0]]*cos(the_tt[[0,0,1,1,0]])
                tyys = the_rr[[0,1,1,0,0]]*sin(the_tt[[0,0,1,1,0]])
                if product(minmax(txxs)-xrange[0]) lt 0 then continue
                if product(minmax(tyys)-yrange[0]) lt 0 then continue
                if product(minmax(txxs)-xrange[1]) lt 0 then continue
                if product(minmax(tyys)-yrange[1]) lt 0 then continue
                the_color = (data[xid,yid] eq 0)? zero_color: sgcolor(zzs[xid,yid], ct=count_ct)
                polyfill, txxs, tyys, color=the_color
            endfor
        endfor

        ; Circle for earth.
        nangle = 50
        angles = smkarthm(0,2*!dpi,nangle,'n')
        circle_x = cos(angles)
        circle_y = sin(angles)

        ; Add magnetopause.
        magn_test = fltarr(nangle,3)
        magn_test[*,0] = circle_x*30
        magn_test[*,1] = circle_y*30
        pdyn = 10
        tmp = check_if_in_magn(magn_test, magn_pos=magn_pos, dynamic_pressure=pdyn)
        magn_pos = magn_pos[*,[0,2,1]]
        magn_pos = [magn_pos,magn_pos]
        magn_pos[nangle:nangle*2-1,1] *= -1


        ; Add more circles.
        foreach tmp, [5,10,20] do oplot, circle_x*tmp, circle_y*tmp, linestyle=1
        foreach tmp, rxy_range do oplot, circle_x*tmp, circle_y*tmp

        ; Add mlt_range.
        foreach tmp, mlt_range do begin
            tt = tmp*15*constant('rad')+!dpi
            oplot, rxy_range*cos(tt), rxy_range*sin(tt)
        endforeach


        ; Add earth.
        polyfill, circle_x<0, circle_y, color=sgcolor('silver')
        plots, circle_x, circle_y


        ; Add magnetopause.
        oplot, magn_pos[*,0], magn_pos[*,1]

        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, fig_labels[xpanel_id]

    endfor

    if keyword_set(test) then stop
    sgclose


end


fig_check_midnight_df_distribution
end
