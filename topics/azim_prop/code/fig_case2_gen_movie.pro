;+
; Generate movie for case 2.
;-


;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range = einfo.time_range_asi

    figdir = shomedir()+'/'+id+'_event'
    movie_file = shomedir()+'/fig_'+id+'_asf_movie.mp4'

    mlat_range = [61,71]
    mlon_range = [-100d,-55]


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
    color_table = 40
    top_color = 254
    colors = findgen(top_color+1)
    count_range = [10,500]

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
    
    
    ntime = n_elements(times)
    figfiles = strarr(ntime)
    for ii=0, ntime-1 do begin
        ofn = figdir+'/thg_asf_mlonimg_'+time_string(times[ii],tformat='YYYY_MMDD_hhmm_ss.png')        
        if keyword_set(test) then ofn = 0
        figfiles[ii] = ofn

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
            time_string(current_time,tformat='YYYY-MM-DD/hh:mm:ss UT'), color=sgcolor('white')


        ; Draw color bar.
        tpos = colorbar_pos
        sgcolorbar, colors, zrange=zrange, ztitle=ztitle, position=tpos, ct=color_table
        
        if test then stop
        sgclose
    endfor


    if ~keyword_set(test) then begin
        spic2movie, figfiles, movie_file, 'png', 10
        foreach file, figfiles do file_delete, file
        file_delete, figdir
    endif
end
