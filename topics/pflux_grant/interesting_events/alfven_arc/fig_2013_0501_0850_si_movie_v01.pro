;+
; Generate the movie for supplemental information.
;
; Adopted from fig_si_movie
;-

function fig_2013_0501_0850_si_movie_v01, movie_file, event_info=event_info


;---Load data and settings.
    if n_elements(event_info) eq 0 then event_info = _2013_0501_0850_load_data()
test = 1

    probe = event_info['probe']
    prefix = event_info['prefix']

    time_range = event_info['asi_time_range']
    model_setting = event_info['model_setting']
    internal_model = model_setting['internal_model']
    models = model_setting['models']
    models = ['t89']
    asi_setting = event_info['asi_setting']
    site = (asi_setting['sites'])[0]
    time_step = 3d

    if n_elements(movie_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        movie_file = join_path([plot_dir,$
            'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD_hhmm')+'_v01.mp4'])
    endif

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [1,1,0.5,4]

;---ASI.
    mlt_range = [-1,1]*6
    mlat_range = [55,70]
    asi_ct = 49
    top_color = 254
    poss = panel_pos(pansize=[6,3], margins=margins, fig_size=fig_size)
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    zrange = [2e2,1e4]
    log_zrange = alog10(zrange)

    tpos = reform(poss)
    cbpos = tpos
    cbpos[0] = tpos[0]+xchsz*1
    cbpos[2] = tpos[2]-xchsz*1
    cbpos[1] = tpos[3]+ychsz*1
    cbpos[3] = cbpos[1]+ychsz*0.8
    ztitle = 'N-hem ASI (#)'
    log_ztickv = make_bins(log_zrange,1, inner=1)
    ztickv = 10^log_ztickv
    zticks = n_elements(ztickv)-1
    zminor = 10
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    zlog = 1
    asi_linestyle = 1
    label_size = 0.8

    min_mlat = 50
    yrange = [min_mlat,90]
    ystep = 10.
    ytickv0 = make_bins(yrange, ystep, /inner)
    ytickv = (ytickv0-min_mlat)/(90-min_mlat)
    ytickn = string(ytickv0,format='(I0)')
    ;ytickn[1:*:2] = ' '

    xtickv0 = [0.,6,12,18]
    xtickv = xtickv0/24*2*!dpi
    xtickn = string(xtickv0,format='(I02)')

    nangle = 50
    dangle = 2*!dpi/nangle
    angles = make_bins([0,2*!dpi], dangle)

    mlt_image_var = 'thg_asf_mlt_image'
    get_data, mlt_image_var, times, mlt_images
    index = where_pro(times, '[]', time_range)
    times = times[index]
    npx = n_elements(mlt_images[0,0,*])
    mlt_images = mlt_images[index,*,0:npx*0.5-1]

    plot_files = []
    foreach time, times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mlt_image',$
            'thg_asf_mlt_image_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, size=fig_size, magnify=magn, test=test

        mlt_image = reform(mlt_images[time_id,*,*])
        zzs = bytscl(alog10(mlt_image), min=log_zrange[0],max=log_zrange[1], top=top_color)
        sgtv, zzs, ct=asi_ct, position=tpos

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ;---Draw axes.
        plot, [-1,1], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        ;plots, [-1,1], [0,0], linestyle=asi_linestyle
        plots, [0,0], [-1,1], linestyle=asi_linestyle

        foreach val, ytickv, val_id do begin
            txs = val*cos(angles)
            tys = val*sin(angles)
            plots, txs,tys, /data, linestyle=asi_linestyle
            ;if float(ytickn[val_id]) eq min_mlat then continue
            tt = 1-val
            rr = -0.25*!dpi
            tx = tt*cos(rr)
            ty = tt*sin(rr)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size
        endforeach

        plots, [0,0], [-1,1], /data, linestyle=asi_linestyle
        plots, [-1,1], [0,0], /data, linestyle=asi_linestyle
        foreach val, xtickv, val_id do begin
            tmp = val-!dpi*0.5
            tx = cos(tmp)
            ty = sin(tmp)
            msg = xtickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            case msg of
                '00': xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
                '12': xyouts, tx,ty-ychsz*label_size*0.7,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
                '06': xyouts, tx,ty-ychsz*label_size*0.3,/normal, alignment=1, msg, charsize=label_size, color=mltimg_linecolor
                '18': xyouts, tx,ty-ychsz*label_size*0.3,/normal, alignment=0, msg, charsize=label_size, color=mltimg_linecolor
            endcase
        endforeach




    ;---Add footpoint.
        suffix = '_'+internal_model+'_north'
        fmlts = get_var_data(prefix+'fmlt'+suffix, at=time)+24
        fmlats = get_var_data(prefix+'fmlat'+suffix, at=time)
        model_colors = sgcolor(['red','green','blue','black'])
        foreach external_model, models, model_index do begin
            sc_color = event_info['sc_color']
            model_color = model_colors[model_index]
            fmlt = fmlts[model_index]
            fmlat = fmlats[model_index]

            tr = (90-fmlat)/(90-min_mlat)
            tt = (fmlt*15-90)*!dtor
            tx = tr*cos(tt)
            ty = tr*sin(tt)
            tmp = convert_coord(tx, ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            plots, tx,ty,normal=1, psym=6, symsize=0.5, color=model_color

            ; Add model name.
            tx = tpos[0]+xchsz*(1+model_index*5)
            ty = tpos[3]-ychsz*1.5
            xyouts, tx,ty,normal=1, strupcase(external_model), color=model_color
            
        endforeach



    ;---Colorbar
        sgcolorbar, findgen(top_color), horizontal=1, $
            ztitle=ztitle, zrange=log_zrange, ct=asi_ct, position=cbpos, $
            ztickv=log_ztickv, zticks=zticks, ztickname=ztickn, zminor=zminor


    ;---Add label.
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('black')



        if keyword_set(test) then stop
        sgclose
    endforeach

    fig2movie, movie_file, fig_files=plot_files, fps=20
    return, movie_file

end


print, fig_2013_0501_0850_si_movie_v01(event_info=event_info)
end