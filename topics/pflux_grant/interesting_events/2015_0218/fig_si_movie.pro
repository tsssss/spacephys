;+
; Generate the movie for supplemental information.
;
; Adopted from fig_asf_movie and fig_overview.
;-


;---Load data and settings.
    event_info = _2015_0218_02_load_data()
    test = 0

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = time_double(['2015-02-18/02:05','2015-02-18/02:15'])
    models = ['t89']
    model_colors = sgcolor(['red','green','blue','purple'])

    root_dir = srootdir()
    plot_dir = join_path([root_dir,'asf'])
    movie_file = join_path([root_dir,'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD')+'.mp4'])

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [8,4,2,1]

;---ASI.
    sites = event_info['site']
    mlt_range = [-2.5,-0.5]
    mlat_range = [59,69]
    asi_ct = 1
    top_color = 254

    themis_read_mltimg, time_range, sites=sites
    mltimg_var = 'thg_mltimg'
    get_data, mltimg_var, times, mltimgs
    index = where_pro(times, '[]', time_range)
    times = times[index]
    mltimgs = mltimgs[index,*,*]

    plot_files = []
    foreach time, times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mltimg_'+time_string(time,tformat='YYYYYY_MMDD_hhmm_ss')+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        sgopen, plot_file, xsize=5, ysize=4
        tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, margins=margins)

        mltimg = reform(mltimgs[time_id,*,*])
        mlt_bins = get_setting(mltimg_var, 'mlt_bins')
        mlt_index = where_pro(mlt_bins, '[]', mlt_range)
        mlt_bins = mlt_bins[mlt_index]
        mltimg = mltimg[mlt_index,*]
        mlat_bins = get_setting(mltimg_var, 'mlat_bins')
        mlat_index = where_pro(mlat_bins, '[]', mlat_range)
        mlat_bins = mlat_bins[mlat_index]
        mltimg = mltimg[*,mlat_index]

        tpos[3] = tpos[3]-ychsz*2.5
        cbpos = tpos
        cbpos[1] = tpos[3]+ychsz*0.2
        cbpos[3] = cbpos[1]+ychsz*0.5

        zrange = [0,250]
        ztitle = 'Photon count (#)'
        zzs = bytscl(mltimg, min=zrange[0],max=zrange[1], top=top_color)
        sgtv, zzs, ct=asi_ct, position=tpos, resize=1
        sgcolorbar, findgen(top_color), horizontal=1, ztitle=ztitle, zrange=zrange, ct=asi_ct, position=cbpos

        xtitle = 'MLT (h)'
        xrange = mlt_range+24
        xstep = 0.5
        xtickv = make_bins(xrange,xstep,/inner)
        xticks = n_elements(xtickv)-1
        xminor = 5

        ytitle = 'MLat (deg)'
        yrange = mlat_range
        ystep = 5
        ytickv = make_bins(yrange,ystep,/inner)
        yticks = n_elements(ytickv)-1
        yminor = 5

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Add grid.
        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xtitle='', xtickformat='(A1)', xtickv=xtickv, xticks=xticks, xminor=1, xticklen=1, xgridstyle=1, $
            ystyle=1, ylog=0, yrange=yrange, ytitle='', ytickformat='(A1)', ytickv=ytickv, yticks=yticks, yminor=1, yticklen=1, ygridstyle=1, $
            position=tpos, nodata=1, noerase=1, ynozero=1, color=sgcolor('gray')

        ; Draw axes.
        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, $
            ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
            position=tpos, nodata=1, noerase=1, ynozero=1

    ;---Add footpoint.
        foreach model, models, model_id do begin
            color = model_colors[model_id]

            fmlt = get_var_data(prefix+'fmlt_'+model, at=time)+24
            fmlat = get_var_data(prefix+'fmlat_'+model, at=time)
            tmp = convert_coord(fmlt, fmlat, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            plots, tx,ty,normal=1, psym=6, symsize=0.5, color=color

            ; Add model name.
            tx = tpos[0]+xchsz*1+xchsz*5*model_id
            ty = tpos[3]-ychsz*1.5
            xyouts, tx,ty,normal=1, strupcase(model), color=color
        endforeach

    ;---Add label.
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = strupcase(sites[0])+' '+time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('white')

        if keyword_set(test) then stop
        sgclose
    endforeach

    spic2movie, movie_file, plot_files=plot_files

end
