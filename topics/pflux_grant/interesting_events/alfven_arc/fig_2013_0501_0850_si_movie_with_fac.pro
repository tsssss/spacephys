;+
; Generate the movie for supplemental information.
; To print the current, glat/glon.
; Adopted from fig_si_movie
;-

function fig_2013_0501_0850_si_movie_with_fac, movie_file, event_info=event_info


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
        movie_file = join_path([plot_dir,$
            'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD_hhmm')+'_with_fac_v01.mp4'])
    endif

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [1,1,1,3.5]

;---ASI.
    zrange = [0,8e3]
    mlt_range = [-1,1]*6
    mlat_range = [55,70]
    asi_ct = 49
    top_color = 254
    poss = panel_pos(pansize=[6,3], margins=margins, fig_size=fig_size)
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]


    tpos = reform(poss)
    cbpos = tpos
    cbpos[0] = tpos[0]+xchsz*1
    cbpos[2] = tpos[2]-xchsz*1
    cbpos[1] = tpos[3]+ychsz*1
    cbpos[3] = cbpos[1]+ychsz*0.8
    ztitle = 'N-hem ASI (#)'
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
        plot_file = join_path([plot_dir,'thg_asf_mlt_image_with_fac',$
            'thg_asf_mlt_image_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, size=fig_size, magnify=magn, test=test

        mlt_image = reform(mlt_images[time_id,*,*])
        zzs = bytscl(mlt_image, min=zrange[0],max=zrange[1], top=top_color)
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
        sgcolorbar, findgen(top_color), horizontal=1, ztitle=ztitle, zrange=zrange, ct=asi_ct, position=cbpos


    ;---Add label.
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('black')



    ;---Add FAC.
        j_ver = get_var_data('thg_j_ver', at=time, limits=lim)
        label_size = 0.7
        pixel_glats = lim.glat_grids
        pixel_glons = lim.glon_grids
        pixel_mlons = lim.mlon_grids
        pixel_mlats = lim.mlat_grids
        pixel_mlts = mlon2mlt(pixel_mlons, time)
        index = where(pixel_mlts gt mlt_range[0] and pixel_mlts lt mlt_range[1] and $
            pixel_mlats gt min_mlat and pixel_mlats lt 90, npixel)
        rrs = (90-pixel_mlats[index])/(90-min_mlat)
        tts = (pixel_mlts[index]*15-90)*constant('rad')
        xxs = rrs*cos(tts)
        yys = rrs*sin(tts)

        zzs = j_ver[index]
        the_glats = pixel_glats[index]
        the_glons = pixel_glons[index]
        fac_zrange = [-1,1]*1e5
        ccs = bytscl(zzs, min=fac_zrange[0], max=fac_zrange[1])
        dc = 80
        ccs[where(zzs ge 0)] = 128+dc
        ccs[where(zzs lt 0)] = 128-dc
        ct = 70
        tmp = smkarthm(0,2*!dpi,30,'n')
        txs = cos(tmp)
        tys = sin(tmp)
        usersym, txs, tys
        
        symszs = (abs(zzs/fac_zrange[1]))^0.25*0.5
        glat_list = list()
        glon_list = list()
        plot, [-1,1], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
            
        for ii=0,npixel-1 do begin
            cc = sgcolor(255-ccs[ii], ct=ct)
            symsz = symszs[ii]
            ;symsz = 0.5
            ;cc = (zzs[ii] ge 0)? sgcolor('blue'): sgcolor('red')
            plots, xxs[ii], yys[ii], color=cc, psym=8, symsize=symsz
            tmp = convert_coord(xxs[ii],yys[ii], data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]
;            msg = ''
;            if glat_list.where(the_glats[ii]) eq !null then begin
;                msg += strtrim(string(the_glats[ii],format='(F4.1)'),2)+'!C'
;                glat_list.add, the_glats[ii]
;            endif
;            if glon_list.where(the_glons[ii]) eq !null then begin
;                msg += strtrim(string(the_glons[ii],format='(F6.1)'),2)+'!C'
;                glon_list.add, the_glons[ii]
;            endif
;            if msg eq '' then begin
;                msg = strtrim(string(zzs[ii]/1e3,format='(F6.1)'),2)
;                xyouts, tx,ty, msg, normal=1, charsize=label_size, color=cc
;            endif else begin
;                xyouts, tx,ty, msg, normal=1, charsize=label_size
;            endelse
        endfor


        if keyword_set(test) then stop
        sgclose
    endforeach

    fig2movie, movie_file, fig_files=plot_files
    return, movie_file

end


print, fig_2013_0501_0850_si_movie_with_fac(event_info=event_info)
end