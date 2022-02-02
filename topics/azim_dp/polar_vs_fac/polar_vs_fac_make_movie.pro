
test = 0

event_list = list()
event_list.add, time_double(['2008-01-19/06:00','2008-01-19/09:00'])
;event_list.add, time_double(['2008-01-21/07:00','2008-01-21/10:00'])

;event_list.add, time_double(['2008-01-29/02:30','2008-01-29/04:30'])
;event_list.add, time_double(['2008-02-10/10:00','2008-02-10/13:00'])
;event_list.add, time_double(['2007-03-16/02:30','2007-03-16/03:30'])
;event_list.add, time_double(['2007-05-07/07:30','2007-05-07/08:00'])
;event_list.add, time_double(['2007-09-23/04:30','2007-09-23/06:00'])
;event_list.add, time_double(['2007-10-03/05:00','2007-10-03/06:00'])
;event_list.add, time_double(['2007-10-30/01:50','2007-10-30/02:50'])
;event_list.add, time_double(['2007-11-05/10:20','2007-11-05/10:50'])
;event_list.add, time_double(['2007-11-08/05:30','2007-11-08/06:30'])
;event_list.add, time_double(['2008-01-05/08:00','2008-01-05/10:00'])

root_dir = join_path([homedir(),'polar_vs_fac'])
foreach time_range, event_list do begin
    event_id = time_string(time_range[0],tformat='YYYY_MMDD_hh')
    plot_dir = join_path([root_dir,event_id])
    movie_file = join_path([root_dir,'polar_vs_fac_movie_'+event_id+'.mp4'])

    themis_read_current_mltimg, time_range
    polar_read_mltimg, time_range

    uvi_var = 'po_mltimg'
    get_data, uvi_var, times, uvis, limits=lim
    uvi_npx = n_elements(reform(uvis[0,*,0]))
    fac_var = 'thg_j_ver_mltimg'
    get_data, fac_var, common_times, facs
    ntime = n_elements(common_times)


;---Interp UVI to common time.
    uvi_mltimgs = fltarr(ntime,uvi_npx,uvi_npx)
    for ii=0,uvi_npx-1 do begin
        for jj=0,uvi_npx-1 do begin
            uvi_mltimgs[*,ii,jj] = interpol(uvis[*,ii,jj],times,common_times)
        endfor
    endfor


;---Get FAC circle MLT image.
    min_mlat = 50d
    mlt_bins = get_setting(fac_var, 'mlt_bins')
    mlat_bins = get_setting(fac_var, 'mlat_bins')

    nmlat_bin = n_elements(mlat_bins)
    npixel = nmlat_bin*2+1
    nmlt_bin = n_elements(mlt_bins)
    nmlat_bin = n_elements(mlat_bins)

    old_image_size = [nmlt_bin,nmlat_bin]
    mlt_2d = mlt_bins # (fltarr(nmlat_bin)+1)
    mlat_2d = (fltarr(nmlt_bin)+1) # mlat_bins
    sphere = 1

    fac_npx = npixel
    fac_mltimgs = fltarr(ntime,fac_npx,fac_npx)
    foreach time, common_times, time_id do begin
        old_image = reform(facs[time_id,*,*])
        get_mlt_image, old_image, mlat_2d, mlt_2d, min_mlat, sphere, mcell=npixel, new_image
        fac_mltimgs[time_id,*,*] = congrid(new_image,fac_npx,fac_npx, interp=0)
    endforeach


;---Plot UVI and FAC side by side.
    fig_xsize = 5
    fig_ysize = 2.5
    magnify = 2
    uvi_zrange = [40d,200]
    fac_zrange = [-1d,1]*8e1
    uvi_ct = 49
    fac_ct = 70
    plot_files = strarr(ntime)
    foreach time, common_times, time_id do begin
        base = 'polar_vs_fac_movie_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'.png'
        plot_file = join_path([plot_dir,base])
        plot_files[time_id] = plot_file
        if keyword_set(test) then plot_file = 0
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, magnify=2
        poss = sgcalcpos(1,2, xpad=0, xchsz=xchsz, ychsz=ychsz, position=[0,0,1,1])

    ;---Settings.
        xtickv0 = [0.,6,12,18]
        xtickv = xtickv0/24*2*!dpi
        xtickn = string(xtickv0,format='(I02)')

        yrange = [min_mlat,90]
        ystep = 10.
        ytickv0 = make_bins(yrange, ystep, /inner)
        ytickv = (ytickv0-min_mlat)/(90-min_mlat)
        ytickn = string(ytickv0,format='(I0)')
        ytickn[1:*:2] = ' '

        label_size = 0.8

        nangle = 50
        dangle = 2*!dpi/nangle
        angles = make_bins([0,2*!dpi], dangle)
        linestyle = 1
        mltimg_linestyle = 0
        mltimg_linecolor = sgcolor('silver')
        top_color = 254
        bottom_color = 1
        reverse_ct = 1

    ;---UVI.
        tpos = reform(poss[*,0])
        timg = reform(uvi_mltimgs[time_id,*,*])
        zz = bytscl(timg, uvi_zrange[0], uvi_zrange[1], top=top_color-bottom_color)+bottom_color
        sgtv, zz, position=tpos, ct=uvi_ct

    ;---Add axis.
        plot, [-1,1], [-1,1], /nodata, /noerase, $
            xstyle=5, ystyle=5, $
            position=tpos

        foreach val, ytickv, val_id do begin
            txs = val*cos(angles)
            tys = val*sin(angles)
            plots, txs,tys, /data, linestyle=mltimg_linestyle, color=mltimg_linecolor
            ;if float(ytickn[val_id]) eq min_mlat then continue
            tt = 1-val
            rr = 0.25*!dpi
            tx = tt*cos(rr)
            ty = tt*sin(rr)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
        endforeach

        plots, [0,0], [-1,1], /data, linestyle=mltimg_linestyle, color=mltimg_linecolor
        plots, [-1,1], [0,0], /data, linestyle=mltimg_linestyle, color=mltimg_linecolor
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

    ;---Add colorbar.
        cbpos = tpos
        cbpos[0] = tpos[0]+xchsz*1
        cbpos[2] = (tpos[0]+tpos[2])*0.5-xchsz*1
        cbpos[3] = tpos[3]-ychsz*2.5
        cbpos[1] = cbpos[3]-ychsz*0.5
        uvi_ztitle = 'Aurora (#)'

        sgcolorbar, reverse(make_bins([bottom_color,top_color],1)), ct=uvi_ct, reverse_ct=reverse_ct, $
            zrange=uvi_zrange, ztitle=' ', position=cbpos, horizontal=1
        tx = (cbpos[0]+cbpos[2])*0.5
        ty = cbpos[3]+ychsz*1.2
        msg = uvi_ztitle
        xyouts, tx,ty,msg, alignment=0.5, normal=1, charsize=label_size


    ;---FAC.
        tpos = reform(poss[*,1])
        timg = reform(fac_mltimgs[time_id,*,*])
        zz = bytscl(timg*1e-3, fac_zrange[0], fac_zrange[1], top=top_color-bottom_color)+bottom_color
        sgtv, zz, position=tpos, ct=fac_ct, reverse_ct=reverse_ct



    ;---Add axis.
        plot, [-1,1], [-1,1], /nodata, /noerase, $
            xstyle=5, ystyle=5, $
            position=tpos

        foreach val, ytickv, val_id do begin
            txs = val*cos(angles)
            tys = val*sin(angles)
            plots, txs,tys, /data, linestyle=mltimg_linestyle, color=mltimg_linecolor
            ;if float(ytickn[val_id]) eq min_mlat then continue
            tt = 1-val
            rr = 0.25*!dpi
            tx = tt*cos(rr)
            ty = tt*sin(rr)
            msg = ytickn[val_id]
            tmp = convert_coord(tx,ty, /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*label_size*0.35
            xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size, color=mltimg_linecolor
        endforeach

        plots, [0,0], [-1,1], /data, linestyle=mltimg_linestyle, color=mltimg_linecolor
        plots, [-1,1], [0,0], /data, linestyle=mltimg_linestyle, color=mltimg_linecolor
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

    ;---Add horizontal currents.
        j_scale = 10000
        j_hor_start_color = sgcolor('silver')
        j_hor_color = sgcolor('dark_slate_grey')
        j_hor_color = sgcolor('black')
        j_hor_thick = 1

        tmp = smkarthm(0,2*!dpi,10, 'n')
        txs = cos(tmp)
        tys = sin(tmp)
        usersym, txs,tys, /fill


        tx = (tpos[0]+tpos[2])*0.5+xchsz*5
        ty = tpos[3]-ychsz*2
        xyouts, tx,ty+ychsz*0.3, 'Horizontal current', /normal, alignment=0.5, charsize=label_size, color=j_hor_color
        len = 1e3/j_scale
        tmp = convert_coord(len,0, /data, /to_normal)
        x_len = tmp[0]-mean(tpos[[0,2]])
        txs = tx+[-1,0]*x_len-xchsz*0.4
        tys = ty+[0,0]-ychsz*0.45+ychsz*0.2
        plots, txs,tys, /normal, color=j_hor_color;, thick=j_hor_thick
        xyouts, tx+xchsz*0.0, ty-ychsz*label_size*0.8+ychsz*0.2, /normal, '1 A/m', color=j_hor_color, charsize=label_size

        var = 'thg_j_hor'
        get_data, var, uts, j_hor
        time_index = where(uts eq time)
        j_hor = reform(j_hor[time_index,*,*,*])
        glat_bins = get_setting(var, 'glatbins')
        glon_bins = get_setting(var, 'glonbins')
        nglat_bin = n_elements(glat_bins)
        nglon_bin = n_elements(glon_bins)
        glat_bins = (fltarr(nglon_bin)+1) # (get_setting(var, 'glatbins'))
        glon_bins = (get_setting(var, 'glonbins')) # (fltarr(nglat_bin)+1)
        ; Convert to mlon/mlt bins.
        geotoapex, glat_bins, glon_bins, '', mlat_bins, mlon_bins
        mlt_bins = mlon2mlt(mlon_bins, time)
        r_bins = (90-mlat_bins)/(90-min_mlat)
        t_bins = (mlt_bins-6)*15*constant('rad')
        x_bins = r_bins*cos(t_bins)
        y_bins = r_bins*sin(t_bins)
        for x_id=0,nglon_bin-1 do begin
            for y_id=0,nglat_bin-1 do begin
                if mlat_bins[x_id,y_id] lt min_mlat then continue
                xx = x_bins[x_id,y_id]
                yy = y_bins[x_id,y_id]
                if product(xx-[-1,1]) ge 0 then continue
                if product(yy-[-1,1]) ge 0 then continue
                jj = reform(j_hor[x_id,y_id,*])
                angle = atan(yy,xx)+!dpi
                jx = jj[0]*cos(angle) + jj[1]*sin(angle)
                jy = jj[0]*sin(angle) - jj[1]*cos(angle)
                x1 = xx+jx/j_scale
                y1 = yy+jy/j_scale
                plots, xx,yy, psym=8, symsize=0.2, color=j_hor_start_color
                oplot, [xx,x1],[yy,y1], color=j_hor_color, thick=j_hor_thick
            endfor
        endfor




    ;---Add colorbar.
        cbpos = tpos
        cbpos[0] = tpos[0]+xchsz*1
        cbpos[2] = (tpos[0]+tpos[2])*0.5-xchsz*1
        cbpos[3] = tpos[3]-ychsz*2.5
        cbpos[1] = cbpos[3]-ychsz*0.5
        fac_ztitle = 'FAC (kA)'

        sgcolorbar, reverse(make_bins([bottom_color,top_color],1)), ct=fac_ct, reverse_ct=reverse_ct, $
            zrange=fac_zrange, ztitle=' ', position=cbpos, horizontal=1
        tx = (cbpos[0]+cbpos[2])*0.5
        ty = cbpos[3]+ychsz*1.2
        msg = fac_ztitle
        xyouts, tx,ty,msg, alignment=0.5, normal=1, charsize=label_size


    ;---Add label.
        tpos = poss[*,0]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        msg = time_string(time, tformat='YYYY-MM-DD/hh:mm:ss')+' UT'
        xyouts, tx,ty,msg, normal=1, charsize=label_size

        if keyword_set(test) then stop
        wait, 0.1
        sgclose
    endforeach

stop
    if keyword_set(test) then stop
    spic2movie, movie_file, plot_files=plot_files

endforeach

end
