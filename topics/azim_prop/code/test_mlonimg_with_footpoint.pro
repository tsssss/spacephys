;+
; Plot the aurora and draw spacecraft and magnetometers on it.
;-


;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data
    
    test = 1


;---Settings.
    einfo = get_var_data(info_var)
    time_range = einfo.aurora['time_range']
    the_time = time_double('2014-08-28/10:09:42')
    ;the_time = time_double('2014-08-28/10:13:00')

    mlat_range = [61,71]
    mlon_range = [-100d,-60]
    
    ; For tracing spacecraft.
    trace_dir = -1
    the_probes = ['thd','the']
    sc_colors = sgcolor(['gold','orange'])
    models = ['t89','t96','t01','t04s']
    model_syms = [6,5,4,1]
    height = 100

    ; For drawing magnetometers.
    mag_sites = ['whit','cmo','pokr','cigo','eagl','fykn','fsim']


;---Load mosaic image.
    mlonimg_var = 'thg_mlonimg'
    get_data, mlonimg_var, times, mlonimgs
    mlon_bins = get_setting(mlonimg_var, 'mlon_bins')
    mlat_bins = get_setting(mlonimg_var, 'mlat_bins')


;---Prepare for plotting.
    aspect_ratio = abs((mlat_range[1]-mlat_range[0])/(mlon_range[1]-mlon_range[0])*2.2)
    figure_size = 3d/[aspect_ratio,1]
    
    magnify = 2
    ticklen = -0.02
    top_color = 254
    colors = findgen(top_color+1)
    count_range = [10,500]
    color_table = 49
    count_range = [0,300]
    label_size = 0.8
    label_color = sgcolor('black')

    
    xtitle = 'MLon (deg)'
    xrange = mlon_range
    xticklen = ticklen
    xminor = 5
    dxmajor = 10
    xtickv = smkarthm(xrange[0],xrange[1],dxmajor,'dx')
    xticks = n_elements(xtickv)-1
    xtickn = string(xtickv,format='(I0)') & xtickn[0] = xtitle+'          '
    
    ytitle = 'MLat (deg)'
    yrange = mlat_range
    yticklen = ticklen*aspect_ratio
    yminor = 5
    dymajor = 5
    ytickv = make_bins(yrange, dymajor, /inner)
    yticks = n_elements(ytickv)-1
    ytickn = string(ytickv,format='(I0)')
    
    ztitle = 'Photon Count'
    zrange = count_range

    wid = 0
    sgopen, wid, xsize=figure_size[0], ysize=figure_size[1], /inch, magnify=magnify
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    figure_pos = sgcalcpos(1, lmargin=8, rmargin=8, tmargin=1, bmargin=4)
    colorbar_pos = [figure_pos[2]+xchsz*1, figure_pos[1], figure_pos[2]+xchsz*2, figure_pos[3]]
    sgclose
    wdelete, wid

    ii = where(times eq the_time, count)
    if count eq 0 then stop
    figdir = sparentdir(srootdir())
    ofn = figdir+'/test_mlonimg_with_footpoint.pdf'
    if keyword_set(test) then ofn = 0

    sgopen, ofn, xsize=figure_size[0], ysize=figure_size[1], /inch, magnify=magnify
    polyfill, [0,1,1,0,0], [0,0,1,1,0], /normal, color=sgcolor('white')

    current_time = times[ii]
    timg = reform(mlonimgs[ii,*,*])
    timg = bytscl(timg, min=count_range[0], max=count_range[1], top=top_color)
    if mlon_range[0] lt min(mlon_bins) then mlon_range[0] = min(mlon_bins)
    if mlon_range[1] gt max(mlon_bins) then mlon_range[1] = max(mlon_bins)
    index = where(mlon_bins ge mlon_range[0] and mlon_bins le mlon_range[1])
    timg = timg[index,*]
    if mlat_range[0] lt min(mlat_bins) then mlat_range[0] = min(mlat_bins)
    if mlat_range[1] gt max(mlat_bins) then mlat_range[1] = max(mlat_bins)
    index = where(mlat_bins ge mlat_range[0] and mlat_bins le mlat_range[1])
    timg = timg[*,index]

    ; Draw figure.
    tpos = figure_pos
    sgtv, timg, position=tpos, /resize, ct=color_table

    plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
    if keyword_set(add_grid) then begin
        plot, xrange, yrange, /nodata, position=tpos, $
            /noerase, $
            xstyle=1, xrange=xrange, xticks=xticks, xminor=0, xtickv=xtickv, xtitle='', xgridstyle=1, xtickformat='(A1)', xticklen=1, $
            ystyle=1, yrange=yrange, yticks=yticks, yminor=0, ytickv=ytickv, ytitle='', ygridstyle=1, ytickformat='(A1)', yticklen=1, $
            color=sgcolor('white')
    endif

    ; Add axes.
    axis, yaxis=0, ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, yticklen=yticklen, ytitle='', ytickname=ytickn
    axis, yaxis=1, ystyle=1, yrange=yrange, yticks=yticks, yminor=yminor, ytickv=ytickv, yticklen=yticklen, ytitle='', ytickformat='(A1)'
    axis, xaxis=0, xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xticklen=xticklen, xtitle='', xtickname=xtickn
    axis, xaxis=1, xstyle=1, xrange=xrange, xticks=xticks, xminor=xminor, xtickv=xtickv, xticklen=xticklen, xtitle='', xtickformat='(A1)'
    ; Add title.
    ;xyouts, tpos[0]-xchsz*6, tpos[1]-ychsz*1.8, /normal, xtitle, alignment=0
    xyouts, tpos[0]-xchsz*5, (tpos[1]+tpos[3])*0.5, /normal, ytitle, alignment=0.5, orientation=90
    ; Add label.
    xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, $
        time_string(current_time,tformat='YYYY-MM-DD/hh:mm:ss UT'), color=label_color


    ; Draw color bar.
    tpos = colorbar_pos
    sgcolorbar, colors, zrange=zrange, ztitle=ztitle, position=tpos, ct=color_table
    
    
;---Add magnetometers.
    psym = 6
    mag_color = sgcolor('red')
    tpos = figure_pos
    plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
    site_infos = themis_read_mag_metadata(sites=mag_sites)
    foreach site, mag_sites, ii do begin
        mlon = site_infos[ii].mlon
        mlat = site_infos[ii].mlat
        plots, mlon, mlat, psym=psym, color=mag_color, symsize=0.5
    endforeach
    
    
;---Add footpoint.
    
    tpos = figure_pos
    plot, xrange, yrange, /nodata, /noerase, position=tpos, xstyle=5, ystyle=5
    foreach sc, the_probes, ii do begin
        foreach model, models, jj do begin
            fmlon_var = sc+'_fmlon_'+model
            fmlat_var = sc+'_fmlat_'+model
            if tnames(fmlon_var) eq '' or tnames(fmlat_var) eq '' then begin
                r_var = sc+'_r_gsm'
                read_geopack_info, r_var, model=model, direction=trace_dir, h0=height
            endif
            fmlon = get_var_data(fmlon_var, at=the_time)
            fmlat = get_var_data(fmlat_var, at=the_time)
            ;print, fmlon, fmlat
            plots, fmlon, fmlat, psym=model_syms[jj], color=sc_colors[ii]
        endforeach
    endforeach
    
    foreach model, models, jj do begin
        tx = tpos[0]+xchsz*1+5*xchsz*jj
        ty = tpos[3]-1*ychsz
        plots, tx,ty+ychsz*0.3*label_size, /normal, psym=model_syms[jj], symsize=0.5, color=label_color
        xyouts, tx+xchsz*1,ty, strupcase(model), /normal, charsize=label_size, color=label_color
    endforeach
    
    foreach sc, the_probes, jj do begin
        tx = tpos[0]+xchsz*1+5*xchsz*jj
        ty = tpos[3]-2*ychsz
        xyouts, tx,ty, strupcase(sc), /normal, charsize=label_size, color=sc_colors[jj]
    endforeach
    tx = tpos[0]+xchsz*1+5*xchsz*n_elements(the_probes)
    xyouts, tx,ty, 'GMAG', /normal, charsize=label_size, color=mag_color

    if test then stop
    sgclose


end
