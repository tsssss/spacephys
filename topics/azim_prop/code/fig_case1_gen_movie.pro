;+
; Generate movie for case 1.
;-
;

    test = 0
    reload = 0

;---Settings.
    time = time_double(['2014-08-28/04:55','2014-08-28/05:07'])
    ;time = time_double(['2014-08-28/05:05','2014-08-28/05:06'])
    sites = ['pina','kapu','snkq']
    min_elevs = [5,10,10]
    mlon_range = [-50,10]
    mlat_range = [59,69]
    
    id = '2014_0828_05'
    figdir = shomedir()+'/'+id+'_event'
    if keyword_set(add_footpoint) then begin
        if keyword_set(all_models) then begin
            movie_file = shomedir()+'/fig_'+id+'_asf_movie_all_models.mp4'
        endif else begin
            movie_file = shomedir()+'/fig_'+id+'_asf_movie_chosen_model.mp4'
        endelse
    endif else begin
        movie_file = shomedir()+'/fig_'+id+'_asf_movie.mp4'
    endelse
    
    
;---Load data.
    load = 0
    if tnames('thg_mlonimg') eq '' then load = 1
    if keyword_set(reload) then load = 1
    if load then themis_read_mlonimg, time, sites=sites, $
        min_elevs=min_elevs, mlon_range=mlon_range, mlat_range=mlat_range, renew_file=renew_file

    mos_info = get_var_data('thg_mlonimg_info')
    mlon_bins = mos_info.mlon_bins
    mlat_bins = mos_info.mlat_bins
    
    get_data, 'thg_mlonimg', times, mlonimgs
    

;---Prepare for plotting.
    ticklen = -0.02
    tplot_options, 'xticklen', ticklen
    tplot_options, 'yticklen', ticklen*0.5
    
    figure_size = [8d,3]
    sgopen, 0, xsize=figure_size[0], ysize=figure_size[1], /inch, magnify=2
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    wdelete, 0
    
    figure_pos = [0.10,0.15,0.90,0.95]
    colorbar_pos = [figure_pos[2]+xchsz*1, figure_pos[1], figure_pos[2]+xchsz*2, figure_pos[3]]
    
    top_color = 254
    colors = findgen(top_color+1)

    ntime = n_elements(times)
    figfiles = []
    for ii=0, ntime-1 do begin
        ofn = figdir+'/thg_asf_mlonimg_'+site+'_'+time_string(times[ii],tformat='YYYY_MMDD_hhmm_ss.png')
        if test eq 1 then ofn = 0 else figfiles = [figfiles, ofn]
        sgopen, ofn, xsize=figure_size[0], ysize=figure_size[1], /inch, magnify=2
        polyfill, [0,1,1,0,0], [0,0,1,1,0], /normal, color=sgcolor('white')

        current_time = times[ii]

        ; mlon is the x-axis.
        timg = reform(mlonimgs[ii,*,*])
        timg = bytscl(timg, min=count_range[0], max=count_range[1], top=top_color)
        index = where(mlon_bins ge mlon_range[0] and mlon_bins le mlon_range[1])
        timg = timg[index,*]
        index = where(mlat_bins ge mlat_range[0] and mlat_bins le mlat_range[1])
        timg = timg[*,index]

        xtitle = 'MLon (deg)'
        ytitle = 'MLat (deg)'
        xrange = mlon_range
        yrange = mlat_range
        xminor = 5
        xtickv = smkarthm(xrange[0],xrange[1],xminor,'dx')
        xticks = n_elements(xtickv)-1
        xtickn = string(xtickv,format='(I0)') & xtickn[0] = ' '; & xtickn[xticks] = ' '
        xticklen = -0.02
        yminor = 2
        ytickv = smkarthm(yrange[0],yrange[1],yminor,'dx')
        yticks = n_elements(ytickv)-1
        ytickn = string(ytickv,format='(I0)'); & ytickn[0] = ' ' & ytickn[yticks] = ' '
        yticklen = -0.01

        tpos = figure_pos
        sgtv, timg, position=tpos, /resize, ct=ct
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
        xyouts, tpos[0]-xchsz*6, tpos[1]-ychsz*1.8, /normal, xtitle, alignment=0
        xyouts, tpos[0]-xchsz*5, (tpos[1]+tpos[3])*0.5, /normal, ytitle, alignment=0.5, orientation=90
        ; Add label.
        xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.5, /normal, $
            strupcase(site)+' '+time_string(current_time,tformat='YYYY-MM-DD/hh:mm:ss UT'), color=sgcolor('white')

        ; Add s/c footpoint.
        if keyword_set(add_footpoint) then begin
            probes = einfo.probes
            probe_colors = sgcolor(['magenta','red'])
            pres = einfo.pres
            models = einfo.models
            model_index = where(models eq einfo.model_chosen)
            if keyword_set(all_models) then model_index = indgen(n_elements(models))
            model_symbols = [7,4,1,6]

            foreach pre0, pres, jj do begin
                probe = probes[jj]
                get_data, pre0+'fpt_mlon', uts, data
                fptmlon = sinterpol(data, uts, current_time)
                get_data, pre0+'fpt_mlat', uts, data
                fptmlat = sinterpol(data, uts, current_time)
                nmodel = n_elements(model_index)
                for kk=0, nmodel-1 do begin
                    tk = model_index[kk]
                    plots, fptmlon[tk], fptmlat[tk], psym=model_symbols[tk], color=probe_colors[jj], symsize=0.8

                    tx = figpos[0]+xchsz*25+xchsz*kk*6
                    ty = figpos[3]+ychsz*0.5
                    xyouts, tx, ty, /normal, color=sgcolor('black'), strupcase(models[tk])
                    plots, tx-xchsz*1, ty+ychsz*0.35, /normal, psym=model_symbols[tk], symsize=0.8

                    ; Add legend for sc.
                    xyouts, tpos[0]+xchsz*8*jj, tpos[3]+ychsz*0.5, /normal, color=probe_colors[jj], 'RBSP-'+strupcase(probe)
                endfor
            endforeach
        endif


        tpos = colorbar_pos
        sgcolorbar, colors, zrange=count_range, ztitle='Photon Count', position=tpos, ct=ct
        if test then stop
        sgclose
    endfor

    if n_elements(figfiles) ne 0 then begin
        spic2movie, figdir, movie_file, 'png', 10
        foreach file, figfiles do file_delete, file
        file_delete, figdir
    endif
end