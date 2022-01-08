;+
; Plot magnetic field and SC position in SM.
;-


;---Settings.
    time_range = time_double(['2015-02-18/02:05','2015-02-18/02:15'])
    plot_time = time_double('2015-02-18/02:10')
    probe = 'a'
    best_model = 't89'
    mlat_range = [61.5,64.5]
    plot_mlats = make_bins([60,70], 2, /inner)
    xrange = [0,8]
    yrange = [-0.5,4]
    fig_file = join_path([homedir(),'Dropbox','mypapers','pflux_vs_beads','fig_2015_0218_orbit.pdf'])
    test = 0

    psym = 1
    xticklen_chsz = -0.20
    yticklen_chsz = -0.30


    prefix = 'rbsp'+probe+'_'
    r_var = prefix+'r_gsm'
    if check_if_update(r_var, time_range) then begin
        rbsp_read_orbit, time_range, probe=probe, coord='gsm'
    endif

    fmlt_var = prefix+'fmlt_'+best_model
    fmlon_var = prefix+'fmlon_'+best_model
    fmlat_var = prefix+'fmlat_'+best_model
    if check_if_update(fmlt_var, time_range) then begin
        read_geopack_info, r_var, model=best_model, direction=-1, refine=1
    endif
    plot_mlt = get_var_data(fmlt_var, at=plot_time)
    plot_mlon = get_var_data(fmlon_var, at=plot_time)


    r_sm = cotran(get_var_data(r_var, at=plot_time), plot_time, 'gsm2sm')


    sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    pan_xsize = 2   ; inch.
    pan_ysize = abs(pan_xsize/total(xrange*[-1,1])*total(yrange*[-1,1]))
    margins = [8,4,2,1]
    fig_xsize = pan_xsize+total(margins[[0,2]])*abs_xchsz
    fig_ysize = pan_ysize+total(margins[[1,3]])*abs_ychsz
    tpos = double(margins)
    tpos[[0,2]] *= abs_xchsz/fig_xsize
    tpos[[1,3]] *= abs_ychsz/fig_ysize
    tpos[[2,3]] = 1-tpos[[2,3]]

    if keyword_set(test) then fig_file = 0
    sgopen, fig_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz


    ; Set up coord.
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos


    ; Add B field.
    deg = constant('deg')
    rad = constant('rad')
    tilt = geopack_recalc(plot_time)
    dir = 1
    par_var = best_model+'_par'
    par = get_var_data(par_var, at=plot_time)
    fline_color = sgcolor('silver')
    foreach plot_mlat, plot_mlats do begin
        lat = plot_mlat*rad
        lon = plot_mlon*rad
        r_mag = [cos(lat)*[cos(lon),sin(lon)],sin(lat)]
        r_gsm = cotran(r_mag, plot_time, 'mag2gsm')

        geopack_trace, r_gsm[0],r_gsm[1],r_gsm[2], $
            dir, par, xf,yf,zf, fline=fline, t89=1
        npoint = n_elements(fline[*,0])
        fline = cotran(fline, plot_time+dblarr(npoint), 'gsm2sm')
        oplot, snorm(fline[*,0:1]), fline[*,2], color=fline_color
    endforeach

;    plot_mlat = get_var_data(fmlat_var, at=plot_time)
;    lat = plot_mlat*rad
;    lon = plot_mlon*rad
;    r_mag = [cos(lat)*[cos(lon),sin(lon)],sin(lat)]
;    r_gsm = cotran(r_mag, plot_time, 'mag2gsm')
;    geopack_trace, r_gsm[0],r_gsm[1],r_gsm[2], $
;        dir, par, xf,yf,zf, fline=fline, t89=1
;    npoint = n_elements(fline[*,0])
;    fline = cotran(fline, plot_time+dblarr(npoint), 'gsm2sm')
;    tx = snorm(fline[*,0:1])
;    ty = fline[*,2]
;    oplot, tx, ty



    ; Add SC position.
    tx = snorm(r_sm[0:1])
    ty = r_sm[2]
    plots, tx,ty,/data, psym=psym
    tmp = convert_coord(tx,ty, /data, /to_normal)
    tx = tmp[0]
    ty = tmp[1]
    xyouts, tx+xchsz*0.5,ty,/normal, 'RBSP-'+strupcase(probe)


    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    ;polyfill, circle_x>0, circle_y, color=sgcolor('silver')
    oplot, circle_x, circle_y


    ; Add grids.
    grid_linestyle = 1
    foreach tx, make_bins(xrange, 2, /inner) do begin
        oplot, tx+[0,0], yrange, linestyle=grid_linestyle
    endforeach
    foreach ty, make_bins(yrange, 2, /inner) do begin
        oplot, xrange, ty+[0,0], linestyle=grid_linestyle
    endforeach


    xtitle = 'SM Rxy (Re), at MLT = '+$
        strtrim(string(plot_mlt,format='(F5.1)'),2)+' h'
    xstep = 2
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)
    xminor = 4
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])

    ytitle = 'SM Z (Re)'
    ystep = 2
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)
    yminor = 4
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])


    ; Draw box.
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticks=xticks, xtickv=xtickv, xminor=xminor, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticks=yticks, ytickv=ytickv, yminor=yminor, yticklen=yticklen, $
        nodata=1, noerase=1, position=tpos


    ; Add label.
    msg = 'd. '+strupcase(best_model)+' model'
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.9
    xyouts, tx,ty,/normal, msg

    if keyword_set(test) then stop
    sgclose

end
