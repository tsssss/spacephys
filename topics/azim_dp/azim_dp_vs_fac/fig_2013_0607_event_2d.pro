;+
; Plot 2D FAC with westward currents at representative times.
;-

times = time_double('2013-06-07/'+['04:47:20','05:21:20'])
times = time_double('2013-06-07/'+['04:47','04:59','05:21','05:42'])

time_range = time_double(['2013-06-07/02:30','2013-06-07/07:00'])
themis_read_current_mltimg, time_range

test = 0
plot_file = join_path([homedir(),'Dropbox','mypapers','dp_vs_fac','plot',$
    'fig_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_2_v01.pdf'])
if keyword_set(test) then plot_file = 0


;---Settings.
    ct = 70
    reverse_ct = 1
    zrange = [-1,1]*300
    ztitle = 'Vertical current (kA), negative value for upward current'

    xtickv0 = [0.,6,12,18]
    xtickv = xtickv0/24*2*!dpi
    xtickn = string(xtickv0,format='(I02)')

    min_mlat = 40
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
    top_color = 255
    bottom_color = 0

    j_hor_color = sgcolor('dark_slate_grey')
    j_hor_color = sgcolor('black')
;    j_hor_color = sgcolor('forest_green')
;    j_hor_color = sgcolor('dark_magenta')
;    j_hor_color = sgcolor('medium_slate_blue')
    j_hor_thick = 1

    fig_label = 'f-'
    fig_label = ''

    tmp = smkarthm(0,2*!dpi,10, 'n')
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs,tys, /fill


;---Get the size of the figure, and position of the panels.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    ; Get the sizes of the left panels.
    margins = [10.5,4,7,1]
    ; xsize of the main panel.
    panel_x_ratio = 1. ; inch per hour.
    panel_xsize = abs(total(time_range*[-1,1]))/3600*panel_x_ratio
    fig_xsize = panel_xsize+(total(margins[[0,2]]))*abs_xchsz
    fig_xsize = 5
    nypan = 2
    nxpan = n_elements(times)/nypan

    margins = [2,1,8,1]
    xpad = 0.8
    ypad = 0.4
    panel_xsize = (fig_xsize-total(margins[[0,2]])*abs_xchsz-(nxpan-1)*xpad*abs_xchsz)/nxpan
    panel_ysize = panel_xsize
    fig_ysize = panel_ysize*nypan+ypad*(nypan-1)*abs_ychsz+(total(margins[[1,3]]))*abs_ychsz

    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
    poss = sgcalcpos(nypan,nxpan, xpad=xpad,ypad=ypad, margins=margins)

    foreach time, times, time_id do begin
        tmp = array_indices([nxpan,nypan], time_id, /dimension)
        tpos = poss[*,tmp[0],tmp[1]]

        ; Read FAC and horizontal currents.
        var = 'thg_j_ver_mltimg'
        get_data, var, uts, mltimgs
        time_index = where(uts eq time)
        mltimg = mltimgs[time_index,*,*]
        mlt_bins = get_setting(var, 'mlt_bins')
        mlat_bins = get_setting(var, 'mlat_bins')

        nmlat_bin = n_elements(mlat_bins)
        npixel = nmlat_bin*2+1
        nmlt_bin = n_elements(mlt_bins)
        nmlat_bin = n_elements(mlat_bins)

        old_image_size = [nmlt_bin,nmlat_bin]
        mlt_2d = mlt_bins # (fltarr(nmlat_bin)+1)
        mlat_2d = (fltarr(nmlt_bin)+1) # mlat_bins
        sphere = 1
        old_image = reform(mltimg)
        get_mlt_image, old_image, mlat_2d, mlt_2d, min_mlat, sphere, mcell=npixel, new_image
        mltimg_circ = new_image

        zz = bytscl(mltimg_circ*1e-3, min=zrange[0], max=zrange[1], top=top_color)
        sgtv, zz, position=tpos, ct=ct, reverse_ct=reverse_ct



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
        tx = tpos[0]+xchsz*5
        ty = tpos[3]-ychsz*2.3
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
        apexfile = join_path([homedir(),'Projects','idl','spacephys','aurora','image','support','mlatlon.1997a.xdr'])
        geo2apex, glat_bins, glon_bins, mlat_bins, mlon_bins
        mlt_bins = mlon2mlt(mlon_bins, time)
        r_bins = (90-mlat_bins)/(90-min_mlat)
        t_bins = (mlt_bins-6)*15*constant('rad')
        x_bins = r_bins*cos(t_bins)
        y_bins = r_bins*sin(t_bins)
        for x_id=0,nglon_bin-1 do begin
            for y_id=0,nglat_bin-1 do begin
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
                plots, xx,yy, psym=8, symsize=0.2, color=j_hor_color
                plots, [xx,x1],[yy,y1], /data, color=j_hor_color, thick=j_hor_thick
            endfor
        endfor



;        tx = tpos[0]+xchsz*0.5
;        ty = tpos[1]+ychsz*0.3
;        xyouts, tx,ty, /normal, time_string(time)

        tx = tpos[0]+xchsz*0.25
        ty = tpos[3]-ychsz*label_size
        msg = fig_label+string(time_id+1,format='(I0)')+'. '+time_string(time,tformat='hh:mm')
        xyouts, tx,ty,/normal, msg

    endforeach

    cbpos = poss[*,nxpan-1,nypan-1]
    cbpos[1] = poss[1,0,nypan-1]
    cbpos[3] = poss[3,0,0]
    cbpos[0] = cbpos[2]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.8
    sgcolorbar, reverse(make_bins([bottom_color,top_color],1)), ct=ct, reverse_ct=reverse_ct, $
        zrange=zrange, ztitle=ztitle, position=cbpos



;---Done.
    if keyword_set(test) then stop
    sgclose

end
