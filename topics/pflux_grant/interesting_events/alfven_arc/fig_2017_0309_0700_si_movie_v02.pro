;+
; Generate the movie for supplemental information.
;
; Adopted from fig_si_movie
;-

function fig_2017_0309_0700_si_movie_v02, movie_file, event_info=event_info


;---Load data and settings.
test = 0
    id = '2017_0309_0700'
    version = 'v02'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)

;---ASI setting.
    time_range = event_info['time_range']
    asi_setting = event_info.ground['asi_setting']
    asi_sites = sort_uniq(asi_setting.sites)
    time_step = 3d

;---SC setting.
    sc_list = list()
    sc_list.add, event_info.themis.thd
    sc_list.add, event_info.themis.the
    nprobe = n_elements(sc_list)




;---Prepare MLT image for ASI.
    mlt_image_var = 'thg_asf_mlt_image'
    get_data, mlt_image_var, times, asi_mlt_images
    index = where_pro(times, '[]', time_range)
    common_times = times[index]
    ntime = n_elements(common_times)
    npx = n_elements(asi_mlt_images[0,0,*])
    asi_mlt_images = asi_mlt_images[index,0:npx*0.5-1,0:npx*0.5-1]
    

;---Settings.
    if n_elements(movie_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        movie_file = join_path([plot_dir,$
            'thg_asf_movie_'+event_info.id+'_'+version+'.mp4'])
    endif

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [1,0.5,1,4]
    poss = panel_pos(pansize=[3,3], margins=margins, fig_size=fig_size)
    top_color = 254
    label_size = 1

    ; Settings.
    mlt_range = [-1d,0]*6
    min_mlat = 50d
    mlat_range = [min_mlat,90]
    ct_ssusi = 70
    ;ct_ssusi = 49
    ssusi_zrange = [-1,1]*20
    ;ssusi_zrange = [0,1]*2
    ct_asi = 49
    asi_zlog = 1
    asi_zrange = [5e2,2e4]
    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins([-!dpi,0], dangle)

    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0]+1,xrange[1]-1, 2, 'dx')
    xtick_minor = make_bins(xrange, 1, inner=1)
    xtickn = xtickv
    index = where(xtickv lt 0, count)
    if count ne 0 then xtickn[index] += 24
    xtickn = string(xtickn,format='(I02)')
    
    ; ylabels.
    ytick_pos = (-1+1d/12)*!dpi
    ytick_pos = (-1)*!dpi
    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange, 10, inner=1)
    ytick_minor = make_bins(yrange, ystep, inner=1)
    ytickn = string(ytickv,format='(I0)')




;---ASI.
    sgopen, 0, size=fig_size, xchsz=xchsz, ychsz=ychsz
    asi_ct = 49
    asi_zrange = [5e2,2e4]
    asi_log_zrange = alog10(asi_zrange)
    asi_tpos = reform(poss[*,0])
    asi_cbpos = asi_tpos
    asi_cbpos[0] = asi_tpos[0]+xchsz*1
    asi_cbpos[2] = asi_tpos[2]-xchsz*1
    asi_cbpos[1] = asi_tpos[3]+ychsz*0.8
    asi_cbpos[3] = asi_cbpos[1]+ychsz*0.5
    asi_ztitle = 'N-hem ASI (#)'
    asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
    asi_ztickv = 10^asi_log_ztickv
    asi_zticks = n_elements(asi_ztickv)-1
    asi_zminor = 9  ; 10 doesn't work??
    asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')
    asi_zlog = 1
    asi_linestyle = 1
    asi_zticklen = xticklen_chsz*ychsz/(asi_cbpos[3]-asi_cbpos[1])


;---Loop through the times.
    plot_files = []
    foreach time, common_times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mlt_image_'+version,$
            'thg_asf_mlt_image_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'_'+version+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, size=fig_size, magnify=magn, test=test, xchsz=xchsz, ychsz=ychsz

    ;---ASI.
        ; Colorbar
        sgcolorbar, findgen(top_color), horizontal=1, $
            ztitle=asi_ztitle, zrange=asi_zrange, ct=asi_ct, position=asi_cbpos, $
            ztickv=asi_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor, log=asi_zlog, zticklen=asi_zticklen
        
        tpos = asi_tpos
        asi_mlt_image = reform(asi_mlt_images[time_id,*,*])
        zzs = bytscl(alog10(asi_mlt_image), min=asi_log_zrange[0],max=asi_log_zrange[1], top=top_color)
        sgtv, zzs, ct=asi_ct, position=tpos

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Draw axes.
        plot, [-1,0], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        oplot, [0,0], [-1,1], color=sgcolor('silver')
        oplot, [-1,1], [0,0], color=sgcolor('silver')

        ; circles for ytickv.
        foreach yminor, ytick_minor, val_id do begin
            rr = (yminor-min_mlat)/(90-min_mlat)
            txs = rr*cos(angles)
            tys = rr*sin(angles)
            linestyle = 1
            index = where(ytickv eq yminor, count)
            if count ne 0 then linestyle = 0
            oplot, txs,tys, linestyle=linestyle, color=sgcolor('silver')
        endforeach

        ; lines for xickv.
        foreach xminor, xtick_minor, val_id do begin
            linestyle = 1
            index = where(xtickv eq xminor, count)
            if count ne 0 then linestyle = 0
            
            tt = (xminor*15-90)*constant('rad')
            txs = [0,1]*cos(tt)
            tys = [0,1]*sin(tt)

            plots, txs,tys, data=1, linestyle=linestyle, color=sgcolor('silver')
        endforeach

        ; add yticknames.
        foreach yminor, ytickv, val_id do begin
            rr = 1-(yminor-min_mlat)/(90-min_mlat)
            tt = ytick_pos
            tx = rr*cos(tt)
            ty = rr*sin(tt)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,normal=1, alignment=0.5, msg, charsize=label_size
        endforeach

        ; add xticknames.
        foreach xminor, xtickv, val_id do begin
            tmp = (xminor*15-90)*constant('rad')
            rr = xtickn_pos
            tx = rr*cos(tmp)
            ty = rr*sin(tmp)
            msg = xtickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            xyouts, tx,ty-ychsz*0.3,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
        endforeach
        


    ;---Add footpoint.
        foreach sc_info, sc_list, sc_id do begin            
            model_setting = sc_info['model_setting']
            internal_model = model_setting['internal_model']
            external_model = model_setting['external_model']
            models = model_setting.models
            model_index = where(models eq external_model)

            probe = sc_info['probe']
            prefix = sc_info['prefix']
            sc_name = sc_info['sc_name']
            sc_color = sc_info['sc_color']

            suffix = '_'+internal_model+'_north'
            fmlts = get_var_data(prefix+'fmlt'+suffix, at=time)
            fmlats = get_var_data(prefix+'fmlat'+suffix, at=time)

            fmlt = fmlts[model_index]
            fmlat = fmlats[model_index]

            tr = (90-fmlat)/(90-min_mlat)
            tt = (fmlt*15-90)*!dtor
            tx = tr*cos(tt)
            ty = tr*sin(tt)
            tmp = convert_coord(tx, ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            plots, tx,ty,normal=1, psym=6, symsize=1, color=sc_color
            msg = strupcase(sc_name+'-'+probe)
;            if probe eq 'd' then begin
;                xyouts, tx-xchsz*1,ty-ychsz*0.35,normal=1,alignment=1, msg, color=sc_color, charsize=label_size
;            endif else begin
;                xyouts, tx,ty-ychsz*1,normal=1,alignment=0.5, msg, color=sc_color, charsize=label_size
;            endelse
            if probe eq 'd' then tx = tx-xchsz*2
            xyouts, tx-xchsz*0.5,ty-ychsz*1.0, alignment=0.5,normal=1, $
                strupcase(sc_name)+'-'+strupcase(probe), color=sc_color, charsize=sc_label_size
        endforeach


        ; Add panel label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        msg = strjoin(strupcase(asi_sites),' ')+' | '+time_string(time,tformat='hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, charsize=label_size


        if keyword_set(test) then stop
        sgclose
    endforeach

    fig2movie, movie_file, fig_files=plot_files, fps=20
    return, movie_file

end


print, fig_2017_0309_0700_si_movie_v02(event_info=event_info)
end