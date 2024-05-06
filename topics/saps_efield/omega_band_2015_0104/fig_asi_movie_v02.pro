
function fig_asi_movie_v02, test=test, event_info=event_info

    version = 'v02'
    id = '2015_0104'


    if n_elements(event_info) eq 0 then event_info = saps_efield_load_data(id)

    time_range = event_info.time_range
    mlt_image_var = 'thg_asf_mlt_image'


;---RBSP fpt.
    sc_list = list()
    probes = ['b']
    internal_model = 'igrf'
    external_model = 't89'
    sc_colors = sgcolor(['purple','blue'])
    foreach probe, probes, probe_id do begin
        prefix = 'rbsp'+probe+'_'
        sc_name = strupcase('rbsp')
        sc_color = sc_colors[probe_id]
        probe_info = dictionary($
            'prefix', prefix, $
            'probe', probe, $
            'sc_name', sc_name, $
            'sc_color', sc_color, $
            'internal_model', internal_model, $
            'external_model', external_model )
        sc_list.add, probe_info
    endforeach
    

;---Settings.
    add_jver = 0
    mlt_type = 'night'
    asi_zlog = 0
    if asi_zlog eq 0 then asi_zrange = [-1,1]*5e3 else asi_zrange = [5e1,2e4]
    asi_ct = 70

    time_step = 3d
    common_times = make_bins(time_range, time_step)
    ntime = n_elements(common_times)
    asi_sites = (event_info['asi_setting'])['sites']
    
    plot_dir = event_info['plot_dir']
    if n_elements(movie_file) eq 0 then begin
        movie_file = join_path([plot_dir,$
            'thg_asf_movie_'+id+'_'+version+'.mp4'])
    endif

    margins = [1,0.5,1,4]
    if mlt_type eq 'pre_midn' then begin
        pansize = [3,3]
    endif else if mlt_type eq 'post_midn' then begin
        pansize = [3,3]
    endif else if mlt_type eq 'night' then begin
        pansize = [6,3]
    endif else if mlt_type eq 'all_mlt' then begin
        pansize = [6,6]
    endif

    tmp = get_abs_chsz()
    abs_xchsz = tmp[0]
    abs_ychsz = tmp[1]

    uniform_ticklen = -abs_ychsz*0.15
    if add_jver then begin
        poss = panel_pos(nxpan=2, $
            pansize=pansize, margins=margins, fig_size=fig_size, xpad=2)
        asi_tpos = reform(poss[*,0])
        jver_tpos = reform(poss[*,1])
    endif else begin
        poss = panel_pos(pansize=pansize, margins=margins, fig_size=fig_size)
        asi_tpos = reform(poss[*,0])
    endelse
    top_color = 254
    label_size = 1
    xchsz = abs_xchsz/fig_size[0]
    ychsz = abs_ychsz/fig_size[1]

    asi_cbpos = asi_tpos
    asi_cbpos[0] = asi_tpos[0]+xchsz*1
    asi_cbpos[2] = asi_tpos[2]-xchsz*1
    asi_cbpos[1] = asi_tpos[3]+ychsz*0.8
    asi_cbpos[3] = asi_cbpos[1]+ychsz*0.5
    asi_zticklen = uniform_ticklen/xchsz*1/fig_size[0]


    if add_jver then begin
        jver_cbpos = jver_tpos
        jver_cbpos[0] = jver_tpos[0]+xchsz*1
        jver_cbpos[2] = jver_tpos[2]-xchsz*1
        jver_cbpos[1] = jver_tpos[3]+ychsz*0.8
        jver_cbpos[3] = jver_cbpos[1]+ychsz*0.5
        jver_zticklen = uniform_ticklen/xchsz*1/fig_size[0]
    endif


    ; Settings.
    if mlt_type eq 'pre_midn' then begin
        mlt_range = [-1d,0]*6
    endif else if mlt_type eq 'post_midn' then begin
        mlt_range = [0d,1]*6
    endif else if mlt_type eq 'night' then begin
        mlt_range = [-1,1]*6
    endif else if mlt_type eq 'all_mlt' then begin
        mlt_range = [-1,1]*12
    endif
    min_mlat = 50d
    mlat_range = [min_mlat,90]
    color_top = 254
    nangle = 50
    dangle = !dpi/nangle
    if mlt_type eq 'pre_midn' then begin
        angles = make_bins([-!dpi*0.5,0], dangle)
    endif else if mlt_type eq 'post_midn' then begin
        angles = make_bins([-!dpi,-!dpi*0.5], dangle)
    endif else if mlt_type eq 'night' then begin
        angles = make_bins([-!dpi,0], dangle)
    endif else if mlt_type eq 'all_mlt' then begin
        angles = make_bins([-!dpi,!dpi], dangle)
    endif
    ; xlabels.
    xtickn_pos = (90.-min_mlat)/(90.-min_mlat+6)
    xrange = mlt_range
    xstep = 6d
    xtickv = smkarthm(xrange[0]+2,xrange[1]-2, 2, 'dx')
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

;---asi.
    asi_var = 'thg_asf_mlt_image'
    get_data, asi_var, times, mlt_images
    npx = n_elements(mlt_images[0,0,*])
    if mlt_type eq 'pre_midn' then begin
        mlt_images = mlt_images[*,0:npx*0.5-1,0:npx*0.5-1]
    endif else if mlt_type eq 'post_midn' then begin
        mlt_images = mlt_images[*,npx*0.5-1:npx-1,0:npx*0.5-1]
    endif else if mlt_type eq 'night' then begin
        mlt_images = mlt_images[*,*,0:npx*0.5-1]
    endif else if mlt_type eq 'all_mlt' then begin
        ; do nothing. no need to crop.
    endif
    
    image_size = size(reform(mlt_images[0,*,*]), dimensions=1)
    asi_mlt_images = fltarr([ntime,image_size])
    for ii=0,image_size[0]-1 do for jj=0,image_size[1]-1 do asi_mlt_images[*,ii,jj] = interpol(mlt_images[*,ii,jj],times, common_times)
    index = where(finite(asi_mlt_images,nan=1), count)
    if count ne 0 then asi_mlt_images[index] = 0

    if asi_zlog eq 1 then begin
        asi_log_zrange = alog10(asi_zrange)
        asi_log_ztickv = make_bins(asi_log_zrange,1, inner=1)
        asi_ztickv = 10^asi_log_ztickv
        asi_zticks = n_elements(asi_ztickv)-1
        asi_zminor = 9
        asi_ztickn = '10!U'+string(asi_log_ztickv,format='(I0)')

        asi_zzs = bytscl(alog10(asi_mlt_images), min=asi_log_zrange[0],max=asi_log_zrange[1], top=color_top)
    endif else begin
        asi_ztickv = sort_uniq([asi_zrange,0])
        asi_zticks = n_elements(asi_ztickv)-1
        asi_zminor = 10
        asi_ztickn = string(asi_ztickv,format='(I0)')

        asi_zzs = bytscl((asi_mlt_images), min=asi_zrange[0],max=asi_zrange[1], top=color_top)
    endelse
    asi_ztitle = 'N-hem ASI (#)'
    asi_zticks = n_elements(asi_ztickv)-1
    asi_linestyle = 1


;---Jver.
    jver_var = 'thg_j_ver_mlt_image'
    if add_jver then begin
        get_data, jver_var, times, mlt_images, limits=lim
        npx = n_elements(mlt_images[0,0,*])
        if mlt_type eq 'pre_midn' then begin
            mlt_images = mlt_images[*,0:npx*0.5-1,0:npx*0.5-1]
        endif else if mlt_type eq 'post_midn' then begin
            mlt_images = mlt_images[*,npx*0.5-1:npx-1,0:npx*0.5-1]
        endif else if mlt_type eq 'night' then begin
            mlt_images = mlt_images[*,*,0:npx*0.5-1]
        endif else if mlt_type eq 'all_mlt' then begin
            ; do nothing. no need to crop.
        endif
        image_size = size(reform(mlt_images[0,*,*]), dimensions=1)
        jver_mlt_images = fltarr([ntime,image_size])
        for ii=0,image_size[0]-1 do for jj=0,image_size[1]-1 do jver_mlt_images[*,ii,jj] = interpol(mlt_images[*,ii,jj],times, common_times)

        jver_zrange = [-1,1]*80
        jver_ct = 70
        jver_zlog = 0
        jver_ztitle = 'N-hem J vertical (kA)'
        jver_ztickv = smkarthm(jver_zrange[0],jver_zrange[1], 5,'n')
        jver_zticks = n_elements(jver_ztickv)-1
        jver_zminor = 10
        jver_ztickn = string(jver_ztickv,format='(I0)')


        jver_zzs = bytscl(jver_mlt_images, min=jver_zrange[0], max=jver_zrange[1], top=color_top)
    endif


;---Loop through the times.
    plot_files = []
    foreach time, common_times, time_id do begin
        if time_id mod 10 ne 0 then continue
    
        plot_file = join_path([plot_dir,'thg_asf_mlt_image_'+version,$
            'thg_asf_mlt_image_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'_'+version+'.png'])
        plot_files = [plot_files, plot_file]
        if keyword_set(test) then plot_file = 0
        if keyword_set(test) then magn = 1 else magn = 2
        sgopen, plot_file, size=fig_size, magnify=magn, test=test, xchsz=xchsz, ychsz=ychsz

    ;---ASI.
        ; colobar.
        sgcolorbar, findgen(top_color), horizontal=1, $
            ztitle=asi_ztitle, zrange=asi_zrange, ct=asi_ct, position=asi_cbpos, $
            ztickv=asi_ztickv, zticks=asi_zticks, ztickname=asi_ztickn, zminor=asi_zminor, log=asi_zlog, zticklen=asi_zticklen
        
        tpos = asi_tpos
        zzs = reform(asi_zzs[time_id,*,*])
        sgtv, zzs, ct=asi_ct, position=tpos

    ;---Jver.
        if keyword_set(add_jver) then begin
            ; colobar.
            sgcolorbar, findgen(top_color), horizontal=1, $
                ztitle=jver_ztitle, zrange=jver_zrange, ct=jver_ct, position=jver_cbpos, $
                ztickv=jver_ztickv, zticks=jver_zticks, ztickname=jver_ztickn, zminor=jver_zminor, log=jver_zlog, zticklen=jver_zticklen

            tpos = jver_tpos
            zzs = reform(jver_zzs[time_id,*,*])
            sgtv, zzs, ct=jver_ct, position=tpos
        endif

    

    ;---Add axis, labels, etc.
        info_list = list()
        info_list.add, dictionary($
            'msg', ['a) North | white light',strjoin(strupcase(asi_sites),' ')+' | '+$
            time_string(time,tformat='hh:mm:ss')+' UT'], $
            'hemisphere', 'north', $
            'position', asi_tpos, $
            'ct', asi_ct )
        if keyword_set(add_jver) then begin
            info_list.add, dictionary($
                'msg', ['a-2) North | J vertical', ''], $
                'hemisphere', 'north', $
                'position', jver_tpos, $
                'ct', jver_ct )
        endif

        foreach the_info, info_list do begin
            tpos = the_info.position

            ; Add labels, etc.
            ; Draw axes.
            if mlt_type eq 'pre_midn' then begin
                xrange = [-1,0]
                yrange = [-1,0]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,-1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
            endif else if mlt_type eq 'post_midn' then begin
                xrange = [0,1]
                yrange = [-1,0]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
            endif else if mlt_type eq 'night' then begin
                xrange = [-1,1]
                yrange = [-1,0]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,-1], [0,0], color=sgcolor('silver')
                plots, [0,1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
            endif else if mlt_type eq 'all_mlt' then begin
                xrange = [-1,1]
                yrange = [-1,1]
                plot, xrange, yrange, nodata=1, noerase=1, $
                    xstyle=5, ystyle=5, position=tpos
                plots, [0,-1], [0,0], color=sgcolor('silver')
                plots, [0,1], [0,0], color=sgcolor('silver')
                plots, [0,0], [0,-1], color=sgcolor('silver')
                plots, [0,0], [0,1], color=sgcolor('silver')
            endif

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

            ; add sc fpt.
            if n_elements(sc_list) ne 0 then begin
                foreach sc_info, sc_list do begin
                    prefix = sc_info['prefix']
                    probe = sc_info['probe']
                    sc_name = sc_info['sc_name']
                    sc_color = sc_info['sc_color']
                    internal_model = (sc_info.haskey('internal_model'))? sc_info['internal_model']: 'dipole'
                    external_model = (sc_info.haskey('external_model'))? sc_info['external_model']: 'dipole'

                    suffix = '_'+internal_model+'_'+external_model+'_north'
                    fmlt = get_var_data(prefix+'fmlt'+suffix, at=time)+24
                    fmlat = abs(get_var_data(prefix+'fmlat'+suffix, at=time))

                    tr = (90-fmlat)/(90-min_mlat)
                    tt = (fmlt*15-90)*!dtor
                    tx = tr*cos(tt)
                    ty = tr*sin(tt)
                    tmp = convert_coord(tx, ty, /data, /to_normal)
                    tx = tmp[0]
                    ty = tmp[1]
                    plots, tx,ty,normal=1, psym=6, symsize=label_size, color=sc_color
                    msg = strupcase(sc_name)+'-'+strupcase(probe)
                    ;msg = strupcase(probe)
                    xyouts, tx-xchsz*1.2,ty-ychsz*0.5, alignment=1,normal=1, $
                        msg, color=sc_color, charsize=sc_label_size
                endforeach
            endif
            
            
            
            ; Add panel label.
            tx = tpos[0]+xchsz*0.5
            ty = tpos[1]+ychsz*0.2
            msgs = the_info.msg
            xyouts, tx,ty,normal=1, msgs[1], charsize=label_size
            xyouts, tx,ty+ychsz*1,normal=1, msgs[0]
        endforeach
    
        if keyword_set(test) then stop
        sgclose
    endforeach

    fig2movie, movie_file, fig_files=plot_files, fps=25
    return, movie_file

end


print, fig_asi_movie_v02(test=0, event_info=event_info)
end
