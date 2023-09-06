;+
; Generate a movie for this event, with footpoint.
;-


;---Load data and settings.
    event_info = _2015_0218_02_load_data()
    test = 1

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = time_double(['2015-02-18/02:00','2015-02-18/03:00'])
    models = ['t89','t96','t01','t04s']

    time_range = time_double(['2015-02-18/02:00','2015-02-18/02:20'])
    models = ['t89']

    root_dir = srootdir()
    plot_dir = join_path([root_dir,'asf'])
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

    foreach time, times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mltimg_'+time_string(time,tformat='YYYYYY_MMDD_hhmm_ss')+'.png'])
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
        dir = -1
        r_gsm = get_var_data(prefix+'r_gsm', at=time)
        rx = r_gsm[0]
        ry = r_gsm[1]
        rz = r_gsm[2]
        re = constant('re')
        r0 = 100d/re+1
        model_colors = sgcolor(['red','green','blue','purple'])
        igrf_types = [0]
        psyms = [1,6]

        foreach model, models, model_id do begin
            par_var = model+'_par'
            if check_if_update(par_var, time_range) then sgeopack_par, time_range, model
            par = get_var_data(par_var, at=time)

            t89 = (model eq 't89')? 1: 0
            t96 = (model eq 't96')? 1: 0
            t01 = (model eq 't01')? 1: 0
            t04s = (model eq 't04s')? 1: 0
            storm = (model eq 't04s')? 1: 0

            color = model_colors[model_id]

            foreach igrf_type, igrf_types, igrf_id do begin
                psym = psyms[igrf_id]

                geopack_trace, rx,ry,rz, dir, par, xf,yf,zf, r0=r0, igrf=igrf_type, t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm

                f_gsm = [xf,yf,zf]
                f_sm = cotran(f_gsm, time, 'gsm2sm')
                fmlat = asin(f_sm[2]/r0)*deg
                fmlt = atan(f_sm[1],f_sm[0])*deg/15+12
                if fmlt lt 0 then fmlt += 24
                plots, fmlt, fmlat, psym=psym, color=color
                
                tx = tpos[0]+xchsz*1+xchsz*5*model_id
                ty = tpos[3]-ychsz*1.5
                xyouts, tx,ty,normal=1, strupcase(model), color=color
            endforeach
        endforeach
        
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = strupcase(sites[0])+' '+time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('white')

        if keyword_set(test) then stop
        sgclose
    endforeach
    
    movie_file = join_path([root_dir,'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD')+'.mp4'])
    spic2movie, plot_dir, movie_file

end
