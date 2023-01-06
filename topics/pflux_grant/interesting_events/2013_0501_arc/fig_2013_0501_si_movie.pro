;+
; Generate the movie for supplemental information.
;
; Adopted from fig_si_movie
;-

function fig_2013_0501_si_movie, movie_file, event_info=event_info


;---Load data and settings.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
test = 0

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = event_info['asi_time_range']
    model_setting = event_info['model_setting']
    external_model = model_setting['external_model']
    internal_model = model_setting['internal_model']
    models = model_setting['models']
    asi_setting = event_info['asi_setting']
    site = (asi_setting['sites'])[0]

    if n_elements(movie_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        movie_file = join_path([plot_dir,'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD')+'.mp4'])
    endif

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [7,4,3,1]

;---ASI.
    mlt_range = [-2,0.5]
    mlat_range = [58,68]
    asi_ct = 49
    top_color = 254

    xtitle = 'MLT (h)'
    xrange = mlt_range+24
    xstep = 0.5
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtickn = string(xtickv,format='(F4.1)')
    foreach tx, xtickv, ii do begin
        if tx lt 24 then continue
        xtickn[ii] = strtrim(string(xtickv[ii]-24,format='(F4.1)'),2)
    endforeach

    ytitle = 'MLat (deg)'
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    yminor = 5


    zrange = [0,1e4]
    ztitle = 'Brightness Count (#)'

    mltimg_var = 'thg_asf_mlt_image_rect'
    get_data, mltimg_var, times, mltimgs
    index = lazy_where(times, '[]', time_range)
    times = times[index]
    mltimgs = mltimgs[index,*,*]

    plot_files = []
    foreach time, times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mlt_image_rect',$
            'thg_asf_mlt_image_rect_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, xsize=4, ysize=2.4, magnify=magn
        tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz, margins=margins)

        mltimg = reform(mltimgs[time_id,*,*])
        mlt_bins = get_setting(mltimg_var, 'mlt_bins')
        mlt_index = lazy_where(mlt_bins, '[]', mlt_range)
        mlt_bins = mlt_bins[mlt_index]
        mltimg = mltimg[mlt_index,*]
        mlat_bins = get_setting(mltimg_var, 'mlat_bins')
        mlat_index = lazy_where(mlat_bins, '[]', mlat_range)
        mlat_bins = mlat_bins[mlat_index]
        mltimg = mltimg[*,mlat_index]

        tpos[3] = tpos[3]-ychsz*2.5
        cbpos = tpos
        cbpos[1] = tpos[3]+ychsz*0.2
        cbpos[3] = cbpos[1]+ychsz*0.5

        zzs = bytscl(mltimg, min=zrange[0],max=zrange[1], top=top_color)
        sgtv, zzs, ct=asi_ct, position=tpos, resize=1
        sgcolorbar, findgen(top_color), horizontal=1, ztitle=ztitle, zrange=zrange, ct=asi_ct, position=cbpos

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Add grid.
        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xtitle='', xtickformat='(A1)', xtickv=xtickv, xticks=xticks, xminor=1, xticklen=1, xgridstyle=1, xtickname=xtickn, $
            ystyle=1, ylog=0, yrange=yrange, ytitle='', ytickformat='(A1)', ytickv=ytickv, yticks=yticks, yminor=1, yticklen=1, ygridstyle=1, $
            position=tpos, nodata=1, noerase=1, ynozero=1, color=sgcolor('gray')

        ; Draw axes.
        plot, xrange, yrange, $
            xstyle=1, xlog=0, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, xtickname=xtickn, $
            ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
            position=tpos, nodata=1, noerase=1, ynozero=1

    ;---Add footpoint.
        fmlt = get_var_data(prefix+'fmlt_'+internal_model, at=time)+24
        fmlat = get_var_data(prefix+'fmlat_'+internal_model, at=time)
        model_index = where(models eq external_model)
        color = sgcolor('red')
        fmlt = fmlt[model_index]
        fmlat = fmlat[model_index]

        tmp = convert_coord(fmlt, fmlat, /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]
        plots, tx,ty,normal=1, psym=6, symsize=0.5, color=color
        msg = 'RBSP-B'
        xyouts, tx,ty+ychsz*0.5,normal=1, msg, color=color, alignment=0.5

;        ; Add model name.
;        tx = tpos[0]+xchsz*1
;        ty = tpos[3]-ychsz*1.5
;        xyouts, tx,ty,normal=1, strupcase(external_model), color=color


    ;---Add label.
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = strupcase(site)+' '+time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('black')

        if keyword_set(test) then stop
        sgclose
    endforeach

    fig2movie, movie_file, fig_files=plot_files
    return, movie_file

end


print, fig_2013_0501_si_movie(event_info=event_info)
end