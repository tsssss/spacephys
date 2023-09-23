;+
; Generate the movie for supplemental information.
;
; Adopted from fig_si_movie
;-

function alfven_arc_gen_asi_movie, input_time_range, sites=sites, min_elevs=min_elevs, event_info=event_info, $
    test=test, calibration_method=calibration_method, merge_method=merge_method


;---Load data and settings.
    probe = event_info['probe']
    prefix = event_info['prefix']

    internal_model = 'dipole'
    models = 't89'
    time_step = 3d
    version = 'v01'
    time_range = time_double(input_time_range)


;---Prepare MLT image for ASI.
    mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
    get_data, mlt_image_var, times, asi_mlt_images
    index = where_pro(times, '[]', time_range)
    common_times = times[index]
    ntime = n_elements(common_times)
    npx = n_elements(asi_mlt_images[0,0,*])
    asi_mlt_images = asi_mlt_images[index,*,0:npx*0.5-1]


;---Prepare MLT image for vertical J.
    j_ver_var = 'thg_j_ver_mlt_image'
    get_data, j_ver_var, times, jvers, limits=lim
    j_ver_npx = lim.image_size[0]
    j_ver_mltimgs = fltarr(ntime,j_ver_npx,j_ver_npx)
    for ii=0,j_ver_npx-1 do for jj=0,j_ver_npx-1 do j_ver_mltimgs[*,ii,jj] = interpol(jvers[*,ii,jj],times, common_times)
    jver_mlt_images = j_ver_mltimgs[*,*,0:j_ver_npx*0.5-1]
    

;---Settings.
    if n_elements(movie_file) eq 0 then begin
        plot_dir = event_info['plot_dir']
        movie_file = join_path([plot_dir,$
            'thg_asf_movie_'+time_string(time_range[0],tformat='YYYY_MMDD_hhmm')+'_'+version+'.mp4'])
    endif

    xticklen_chsz = -0.2
    yticklen_chsz = -0.5
    margins = [1,1,0.5,4]
    poss = panel_pos(pansize=[6,3], margins=margins, fig_size=fig_size, xpans=[1,1], xpad=2)
    top_color = 254
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    label_size = 1

    nangle = 50
    dangle = !dpi/nangle
    angles = make_bins([-!dpi,0], dangle)

    ; xlabels.
    xrange = [-1d,1]*6
    xstep = 6d
    xtickv = make_bins(xrange, xstep, inner=1)
    xtickv_tt = make_bins(xrange, 1, inner=1)
    xtickn = ['06','00','18']


    ; ylabels.
    min_mlat = 50d
    yrange = [min_mlat,90]
    ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    ytickv_rr = (ytickv-min_mlat)/(90-min_mlat)
    ytickn = string(ytickv,format='(I0)')
    ytickn[1:*:2] = ' '




;---ASI.
    asi_ct = 49
    asi_zrange = [2e2,2e4]
    asi_log_zrange = alog10(asi_zrange)
    asi_tpos = reform(poss[*,0])
    asi_cbpos = asi_tpos
    asi_cbpos[0] = asi_tpos[0]+xchsz*1
    asi_cbpos[2] = asi_tpos[2]-xchsz*1
    asi_cbpos[1] = asi_tpos[3]+ychsz*1
    asi_cbpos[3] = asi_cbpos[1]+ychsz*0.8
    asi_ztitle = 'N-hem ASI (#)'
    asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
    asi_ztickv = 10^asi_log_ztickv
    asi_zticks = n_elements(asi_ztickv)-1
    asi_zminor = 10
    asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')
    asi_zlog = 1
    asi_linestyle = 1


;---Jver.
    jver_zrange = [-1,1]*80
    jver_ct = 70
    jver_tpos = reform(poss[*,1])
    jver_cbpos = jver_tpos
    jver_cbpos[0] = jver_tpos[0]+xchsz*1
    jver_cbpos[2] = jver_tpos[2]-xchsz*1
    jver_cbpos[1] = jver_tpos[3]+ychsz*1
    jver_cbpos[3] = jver_cbpos[1]+ychsz*0.8
    jver_ztitle = 'N-hem J!Nverticl!N (kA)'
    jver_ztickv = smkarthm(jver_zrange[0],jver_zrange[1], 5,'n')
    jver_zticks = n_elements(jver_ztickv)-1
    jver_zminor = 10
    jver_ztickn = string(jver_ztickv,format='(I0)')
    jver_zlog = 0


;---Loop through the times.
    plot_files = []
    foreach time, common_times, time_id do begin
        plot_file = join_path([plot_dir,'thg_asf_mlt_image_v02',$
            'thg_asf_mlt_image_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'_'+version+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, size=fig_size, magnify=magn, test=test

    ;---ASI.
        ; Colorbar
        sgcolorbar, findgen(top_color), horizontal=1, $
            ztitle=asi_ztitle, zrange=asi_log_zrange, ct=asi_ct, position=asi_cbpos, $
            ztickv=asi_log_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor
        
        tpos = asi_tpos
        asi_mlt_image = reform(asi_mlt_images[time_id,*,*])
        zzs = bytscl(alog10(asi_mlt_image), min=asi_log_zrange[0],max=asi_log_zrange[1], top=top_color)
        sgtv, zzs, ct=asi_ct, position=tpos

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Draw axes.
        plot, [-1,1], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        plots, [0,0], [-1,1], linestyle=asi_linestyle
        plots, [-1,1], [0,0], /data, linestyle=asi_linestyle

        ; circles for ytickv.
        foreach rr, ytickv_rr, val_id do begin
            txs = rr*cos(angles)
            tys = rr*sin(angles)
            plots, txs,tys, data=1, linestyle=asi_linestyle
        endforeach

        ; add yticknames.
        foreach rr, ytickv_rr, val_id do begin
            rr = 1-rr
            tt = -0.25*!dpi
            tx = rr*cos(tt)
            ty = rr*sin(tt)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,normal=1, alignment=0.5, msg, charsize=label_size
        endforeach

        ; lines for xickv.
        foreach tt, xtickv_tt, val_id do begin
            tt = (tt*15-90)*constant('rad')
            txs = [0,1]*cos(tt)
            tys = [0,1]*sin(tt)
            plots, txs,tys, data=1, linestyle=asi_linestyle
        endforeach

        ; add xticknames.
        foreach tt, xtickv, val_id do begin
            tmp = (tt*15-90)*constant('rad')
            tx = cos(tmp)
            ty = sin(tmp)
            msg = xtickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            case msg of
                '00': xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
                '12': xyouts, tx,ty-ychsz*label_size*0.7,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
                '06': xyouts, tx,ty-ychsz*label_size*0.3,/normal, alignment=0, msg, charsize=label_size, color=mltimg_linecolor
                '18': xyouts, tx,ty-ychsz*label_size*0.3,/normal, alignment=1, msg, charsize=label_size, color=mltimg_linecolor
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

;            ; Add model name.
;            tx = tpos[0]+xchsz*(1+model_index*5)
;            ty = tpos[3]-ychsz*1.5
;            xyouts, tx,ty,normal=1, strupcase(external_model), color=model_color
            
        endforeach


        ; Add label.
        tx = tpos[0]+xchsz*1
        ty = tpos[1]+ychsz*0.3
        msg = time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,normal=1, msg, color=sgcolor('black')



    ;---Jver.

        ; Colorbar
        sgcolorbar, findgen(top_color), horizontal=1, $
            ztitle=jver_ztitle, zrange=jver_zrange, ct=jver_ct, position=jver_cbpos, $
            ztickv=jver_ztickv, zticks=jver_zticks, ztickname=jver_ztickn, zminor=jver_zminor

        tpos = jver_tpos
        jver_mlt_image = reform(jver_mlt_images[time_id,*,*])
        zzs = bytscl((jver_mlt_image), min=jver_zrange[0],max=jver_zrange[1], top=top_color)
        sgtv, zzs, ct=jver_ct, position=tpos

        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Draw axes.
        plot, [-1,1], [-1,0], /nodata, /noerase, $
            xstyle=5, ystyle=5, position=tpos
        ;plots, [-1,1], [0,0], linestyle=asi_linestyle
        plots, [0,0], [-1,1], linestyle=asi_linestyle
        plots, [-1,1], [0,0], /data, linestyle=asi_linestyle

        ; circles for ytickv.
        foreach rr, ytickv_rr, val_id do begin
            txs = rr*cos(angles)
            tys = rr*sin(angles)
            plots, txs,tys, data=1, linestyle=asi_linestyle
        endforeach

        ; add yticknames.
        foreach rr, ytickv_rr, val_id do begin
            rr = 1-rr
            tt = -0.25*!dpi
            tx = rr*cos(tt)
            ty = rr*sin(tt)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,normal=1, alignment=0.5, msg, charsize=label_size
        endforeach

        ; lines for xickv.
        foreach tt, xtickv_tt, val_id do begin
            tt = (tt*15-90)*constant('rad')
            txs = [0,1]*cos(tt)
            tys = [0,1]*sin(tt)
            plots, txs,tys, data=1, linestyle=asi_linestyle
        endforeach

        ; add xticknames.
        foreach tt, xtickv, val_id do begin
            tmp = (tt*15-90)*constant('rad')
            tx = cos(tmp)
            ty = sin(tmp)
            msg = xtickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]
            case msg of
                '00': xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
                '12': xyouts, tx,ty-ychsz*label_size*0.7,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
                '06': xyouts, tx,ty-ychsz*label_size*0.3,/normal, alignment=0, msg, charsize=label_size, color=mltimg_linecolor
                '18': xyouts, tx,ty-ychsz*label_size*0.3,/normal, alignment=1, msg, charsize=label_size, color=mltimg_linecolor
            endcase
        endforeach


        ; Add footpoint.
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

;            ; Add model name.
;            tx = tpos[0]+xchsz*(1+model_index*5)
;            ty = tpos[3]-ychsz*1.5
;            xyouts, tx,ty,normal=1, strupcase(external_model), color=model_color

        endforeach


        ; Add label.
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

test = 0

event_info = orderedhash()
event_info['2015_0312'] = dictionary($
    'time_range', ['2015-03-12/08:00','2015-03-12/11:00'], $
    'probes', ['b'], $
    'bad_e_time_ranges', list(['2015-03-12/09:04:51','2015-03-12/09:04:54']), $
    'asi_setting', dictionary($
        'sites', ['fykn','mcgr'], $
        'best_site', 'fykn', $
        'min_elevs', float([5,10]), $
        'merge_method', 'max_elev', $
        'calibration_method', 'simple', $
        'mlt_range', [-6d,0], $
        'mlat_range', [55d,70] ) )
    
;event_info['2015_0302_11'] = dictionary($
;    'time_range', ['2015-03-02/10:30','2015-03-02/11:30'], $
;    'probes', ['a'], $
;    ;    'asi_sites', ['nrsq','snkq','gill','rank','chbg'], $
;    'asi_sites', ['nrsq','snkq'], $
;    'asi_min_elevs', !null, $
;    'asi_has_moon', 0 )
;event_info['2015_0105_00'] = dictionary($
;    'time_range', ['2015-01-05/00:10','2015-01-05/02:00'], $
;    'probes', ['a'], $
;;    'asi_sites', ['nrsq','snkq','gill','rank','chbg'], $
;    'asi_sites', ['nrsq','snkq'], $
;    'asi_min_elevs', [1,1]*2.5, $
;    'asi_has_moon', 1 )

foreach info, event_info do begin
    time_range = info.time_range
    
    asi_setting = info.asi_setting
    tmp = alfven_arc_load_ground_data(time_range, asi_setting=asi_setting)
    stop
    
    probes = info.probes
    foreach probe, probes do begin
        tmp = alfven_arc_load_rbsp_data(time_range, probe=probe, _extra=info.tostruct())
    endforeach

    print, alfven_arc_gen_asi_movie(time_range, sites=sites, min_elevs=min_elevs, test=test, merge_method='max_elev', event_info=pinfo)
    stop
endforeach
end