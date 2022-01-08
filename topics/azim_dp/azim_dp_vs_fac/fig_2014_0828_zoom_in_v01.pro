;+
; Show local M-I coupling.
;-


;---Settings.
    test = 0
    test_panels = 0
    plot_file = join_path([homedir(),'Dropbox','mypapers','dp_vs_fac','fig_2014_0828_zoom_in.pdf'])
    if keyword_set(test) then plot_file = 0

    ; For line plots of Pflux, BBFs, etc.
    mid_time_range = time_double(['2014-08-28/10:00','2014-08-28/10:40'])
    mid_probes = 'th'+['a','d','e']
    mid_colors = constant('rgb')
    injection_probes = ['g15','1991-080']
    injection_energy_range = [50,500]   ; keV.
    data_time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])

    ; For add DP.
    dp_info = dictionary($
        'probes', ['thd','the','tha','g15','rbspb','g13'], $
        'dp_times', time_double('2014-08-28/'+['10:10:40','10:13:20','10:15:15','10:20:25','10:44:55','10:47:25']), $
        'omega', 2.1, $
        'angular_width', 10., $
        'plot_time_range', time_double(['2014-08-28/10:07','2014-08-28/10:55']), $
        'plot_mlt_range', [0,7] )
    _2014_0828_10_load_data
    nprobe = n_elements(dp_info.probes)
    dp_mlts = fltarr(nprobe)
    dp_rxys = fltarr(nprobe)
    ndim = 3
    dp_r_sms = fltarr(nprobe,ndim)
    foreach probe, dp_info.probes, probe_id do begin
        time = dp_info.dp_times[probe_id]
        r_gsm = get_var_data(probe+'_r_gsm', at=time)
        r_sm = cotran(r_gsm, time, 'gsm2sm')
        dp_r_sms[probe_id,*] = r_sm
        dp_mlts[probe_id] = pseudo_mlt(r_sm)
        dp_rxys[probe_id] = snorm(r_sm[0:1])
    endforeach
    dp_info['dp_mlts'] = dp_mlts
    dp_info['dp_rxys'] = dp_rxys
    dp_info['dp_r_sms'] = dp_r_sms


    ; For 2D map.
    probes = ['thd','the','g15','tha','1991-080','g13','rbspb']
    short_names = strupcase(['thd','the','g15','tha','Lanl','g13','rbb'])
    nprobe = n_elements(probes)
    pos_times = time_double('2014-08-28/'+['10:10','10:20'])
    coord = 'sm'
    pos_xrange = [2,-15]
    pos_yrange = [2,-12]
    pos_zrange = [0,6]

    ; For mapping.
    model = 't01'
    the_time = time_double('2014-08-28/10:20')
    tilt = geopack_recalc(the_time)
    line_info = list()
    up_current = dictionary('lat_range', [60,67], 'color', sgcolor('salmon'))
    ;down_current = dictionary('lat_range', [68,72], 'color', sgcolor('sky_blue'))
    down_current = dictionary('lat_range', [68,72], 'color', sgcolor('dodger_blue'))
    rad = constant('rad')
    foreach lat, make_bins(up_current.lat_range,1)*rad do line_info.add, dictionary('lat',lat, 'color', up_current.color)
    foreach lat, make_bins(down_current.lat_range,1)*rad do line_info.add, dictionary('lat',lat, 'color', down_current.color)
    foreach lat, [50,60,70,80,90,-80,-90,100,110,120,130]*rad do line_info.add, dictionary('lat',lat, 'color', sgcolor('black'))
    nline = line_info.length


    ; Probe colors.
    top_color = 240
    bottom_color = 90
    probe_colors = smkarthm(bottom_color,top_color,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)


;---Calc panel positions.
    sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    xpads = [9.,11]
    ypad = 0.4
    margins = [7,4,6,1]

    ; left panels.
    left_pan_xsize = 2.
    nleft_pan = 2
    left_ypans = [$
        abs(left_pan_xsize/total(pos_xrange*[-1,1])*total(pos_zrange*[-1,1])),$
        abs(left_pan_xsize/total(pos_xrange*[-1,1])*total(pos_yrange*[-1,1]))]
    left_pan_ysize = total(left_ypans)+abs_ychsz*ypad

    ; middle panels.
    mid_pan_xsize = 2.
    mid_pan_ysize = left_pan_ysize
    nmid_pan = 5
    second_color = sgcolor('red')
    test_times = time_double(['2014-08-28/10:10','2014-08-28/10:20'])

    ; right panels.
    right_pan_ysize = left_pan_ysize
    nright_pan = 4
    right_ypans = [1,1,1,1.5]
    right_ypads = ypad*[1,1,1]
    pan_ysize = (right_pan_ysize-total(right_ypads)*abs_ychsz)/total(right_ypans)*right_ypans[0]
    right_pan_xsize = pan_ysize*2

    ; Figure.
    fig_ysize = left_pan_ysize+total(margins[[1,3]])*abs_ychsz
    fig_xsize = left_pan_xsize+mid_pan_xsize+right_pan_xsize+total(xpads)*abs_xchsz+total(margins[[0,2]])*abs_xchsz

    left_pan_pos = [$
        margins[0]*abs_xchsz,$
        margins[1]*abs_ychsz,$
        margins[0]*abs_xchsz+left_pan_xsize,$
        margins[1]*abs_ychsz+left_pan_ysize ]

    mid_pan_pos = left_pan_pos
    mid_pan_pos[0] = left_pan_pos[2]+xpads[0]*abs_xchsz
    mid_pan_pos[2] = mid_pan_pos[0]+mid_pan_xsize

    right_pan_pos = mid_pan_pos
    right_pan_pos[0] = mid_pan_pos[2]+xpads[1]*abs_xchsz
    right_pan_pos[2] = right_pan_pos[0]+right_pan_xsize

    ; Normalized.
    left_pan_pos[[0,2]] /= fig_xsize
    left_pan_pos[[1,3]] /= fig_ysize
    mid_pan_pos[[0,2]] /= fig_xsize
    mid_pan_pos[[1,3]] /= fig_ysize
    right_pan_pos[[0,2]] /= fig_xsize
    right_pan_pos[[1,3]] /= fig_ysize


    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz

    left_poss = sgcalcpos(nleft_pan, position=left_pan_pos, ypad=ypad, ypans=left_ypans)
    mid_poss = sgcalcpos(nmid_pan, position=mid_pan_pos, ypad=ypad)
    right_poss = sgcalcpos(nright_pan, position=right_pan_pos, ypad=right_ypads, ypans=right_ypans)

    if keyword_set(test_panels) then begin
        for ii=0,nleft_pan-1 do begin
            tpos = left_poss[*,ii]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, xtitle='X', ytitle='Y', $
                position=tpos, nodata=1, noerase=1
        endfor

        for ii=0,nmid_pan-1 do begin
            tpos = mid_poss[*,ii]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, xtitle='X', ytitle='Y', $
                position=tpos, nodata=1, noerase=1
        endfor

        for ii=0,nright_pan-1 do begin
            tpos = right_poss[*,ii]
            plot, [0,1],[0,1], $
                xstyle=1, ystyle=1, xtitle='X', ytitle='Y', $
                position=tpos, nodata=1, noerase=1
        endfor
        cbpos = right_poss[*,0]
        cbpos[1] = right_poss[1,nright_pan-2]
        cbpos[0] = right_poss[2,0]+xchsz*0.8
        cbpos[2] = cbpos[0]+xchsz*0.5
        sgcolorbar, position=cbpos, zrange=[0,1], ztitle='Z'

        if keyword_set(test) then stop
        sgclose
        stop
    endif

    xticklen_chsz = -0.20
    yticklen_chsz = -0.40


;---Left panels.

    ; Panel a-1.
    tpos = left_poss[*,0]
    label_size = 0.8

    times = the_time+dblarr(nline)
    ndim = 3
    r_sm = dblarr(nline,ndim)
    foreach line, line_info, ii do begin
        lat = line.lat
        r_sm[ii,*] = [-cos(lat),0,sin(lat)]
    endforeach

    h0 = 100.   ; km.
    re = constant('re')
    r0 = 1+h0/re
    r_sm *= r0
    r_gsm = cotran(r_sm, times, 'sm2gsm')
    flines = list()


    xrange = pos_xrange
    yrange = pos_zrange
    xstep = 5
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = ystep

    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    ytitle = 'SM Z (Re)'

    ; Set up coord.
    plot, xrange, yrange, /nodata, $
        xstyle=1, xrange=xrange, xticklen=1, xgridstyle=1, xtickformat='(A1)', $
        xticks=xticks, xtickv=xtickv, xminor=1, $
        ystyle=1, yrange=yrange, yticklen=1, ygridstyle=1, ytickformat='(A1)', $
        yticks=yticks, ytickv=ytickv, yminor=1, $
        position=tpos, /iso, color=sgcolor('silver')


    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)>min(yrange)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    plots, circle_x, circle_y


    ; Prepare model info.
    sgeopack_par, the_time+[-1,1]*600, model
    par = get_var_data(model+'_par', at=the_time)
    routine = (model eq 't04s')? 'ts04': model
    routine = 'geopack_'+routine
    t89 = (model eq 't89')? 1: 0
    t96 = (model eq 't96')? 1: 0
    t01 = (model eq 't01')? 1: 0
    t04s = (model eq 't04s')? 1: 0
    storm = (model eq 't04s')? 1: 0


    min_bmag = 20.
    min_bmag_color = sgcolor('dark_gray')
    min_bmag_color = sgcolor('silver')
    foreach line, line_info, ii do begin
        dir = (r_gsm[ii,2]>0)? 1: -1
        geopack_trace, r_gsm[ii,0], r_gsm[ii,1], r_gsm[ii,2], dir, par, xf, yf, zf, $
            fline=fline, /igrf, r0=r0, t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm
        ndata = n_elements(fline[*,0])
        fline = cotran(fline, the_time+dblarr(ndata), 'gsm2sm')
        oplot, fline[*,0], fline[*,2], color=line.color

        ndata = n_elements(fline[*,0])
        b_gsm = fltarr(ndata,ndim)
        for jj=0,ndata-1 do begin
            geopack_igrf_gsm, fline[jj,0],fline[jj,1],fline[jj,2], bxp,byp,bzp
            call_procedure, routine, par, fline[jj,0],fline[jj,1],fline[jj,2], tbx,tby,tbz
            b_gsm[jj,*] = [bxp,byp,bzp]+[tbx,tby,tbz]
        endfor
        bmag = snorm(b_gsm)
        index = where(bmag le min_bmag, count)
        if count ne 0 then oplot, fline[index,0], fline[index,2], color=min_bmag_color
    endforeach


    plot, xrange, yrange, /nodata, /noerase, $
        xstyle=1, xrange=xrange, xticklen=xticklen, xtitle=' ', $
        xticks=xticks, xtickv=xtickv, xminor=xminor, xtickformat='(A1)', $
        ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
        yticks=yticks, ytickv=ytickv, yminor=yminor, $
        position=tpos, /iso


    yspace = 0.9
    msg = 'J Downward!C '+strjoin(string(down_current.lat_range,'(I0)'),'-')+' deg'
    tx = -8
    ty = 4
    xyouts, tx,ty,alignment=0.5, /data, msg, color=sgcolor('blue'), charsize=label_size
    msg = 'J Upward!C'+strjoin(string(up_current.lat_range,'(I0)'),'-')+' deg'
    tx = -5.7
    ty = 2.0
    xyouts, tx,ty,alignment=0.5, /data, msg, color=sgcolor('red'), charsize=label_size
    msg = 'Current!Cclosure?'
    tx = -12.5
    ty = 1.5
    xyouts, tx,ty,alignment=0.5, /data, msg, color=min_bmag_color, charsize=label_size


    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.9
    xyouts, tx,ty, /normal, 'a-1. XZ'


    ; Panel a-2.
    tpos = left_poss[*,1]

    pos_coord = list()
    foreach probe, probes do begin
        r_gsm = get_var_data(probe+'_r_gsm', at=pos_times)
        r_coord = cotran(r_gsm, pos_times, 'gsm2'+coord)
        if n_elements(r_coord) eq 3 then r_coord = reform(r_coord,1,3)
        pos_coord.add, r_coord
    endforeach
    pos_coord = pos_coord.toarray()

    xrange = pos_xrange
    yrange = pos_yrange
    zrange = pos_zrange

    xstep = 5
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = ystep

    zstep = 2
    ztickv = make_bins(zrange, zstep, /inner)
    zticks = n_elements(ztickv)-1
    zminor = zstep

    xtitle = strupcase(coord)+' X (Re)'
    ytitle = strupcase(coord)+' Y (Re)'
    ztitle = strupcase(coord)+' Z (Re)'


;---Setup coord.
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    plot, xrange, yrange, /iso, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        nodata=1, noerase=1, position=tpos

    plots, xrange, [0,0], linestyle=1
    plots, [0,0], yrange, linestyle=1

    ; Add earth.
    nangle = 50
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    polyfill, circle_x<0, circle_y, color=sgcolor('silver')
    polyfill, circle_x>0, circle_y, color=sgcolor('white')
    plots, circle_x, circle_y
    foreach r, [5,10] do oplot, circle_x*r, circle_y*r, linestyle=1


    ; Add DP.
    rad = constant('rad')
    deg = constant('deg')
    arrow_hsize = keyword_set(test)? 4:80
    arrow_solid = 1

    nangle = 50
    omega_2d = dp_info.omega

    rxy_extent = minmax(dp_info.dp_rxys)
    mlt_extent = minmax(dp_info.dp_mlts)
    angle_extent = (mlt_extent*15+180)*rad
    ww_angle = dp_info.angular_width*rad

    ; Shape.
    mlt_extent_color = sgcolor('blue')
    shape_color = sgcolor('light_cyan')
    rr = 0
    the_angle = angle_extent[0]*(1-rr)+angle_extent[1]*rr
    angles = smkarthm(the_angle-ww_angle*0.5, the_angle+ww_angle*0.5,nangle, 'n')
    trs = smkarthm(rxy_extent[1],rxy_extent[0],nangle,'n')
    txs = [rxy_extent[1]*cos(angles),trs*cos(angles[-1]),rxy_extent[0]*reverse(cos(angles)),reverse(trs)*cos(angles[0])]
    tys = [rxy_extent[1]*sin(angles),trs*sin(angles[-1]),rxy_extent[0]*reverse(sin(angles)),reverse(trs)*sin(angles[0])]
    ;txs = min(panel_xrange)>txs<max(panel_xrange)
    ;tys = min(panel_yrange)>tys<max(panel_yrange)
    polyfill, txs, tys, /data, color=shape_color
    plots, txs,tys, /data, color=mlt_extent_color

    ; MLT extent.
    rr = 0.8
    the_rxy = rxy_extent[0]*(1-rr)+rxy_extent[1]*rr
    angles = smkarthm(angle_extent[0], angle_extent[1], nangle, 'n')
    txs = the_rxy*cos(angles)
    tys = the_rxy*sin(angles)
    plots, txs, tys, /data, color=mlt_extent_color
    arrow, txs[-2],tys[-2], txs[-1],tys[-1], /data, $
        hsize=arrow_hsize*1.5, solid=arrow_solid, color=mlt_extent_color
    tmp = convert_coord(0,ychsz, /normal,/to_data)-convert_coord(0,0, /normal,/to_data)
    drxy = tmp[1]*0.3
    tmp = the_rxy+[-1,1]*drxy
    plots, tmp*cos(angles[0]), tmp*sin(angles[0]), /data, color=mlt_extent_color

    tx = the_rxy*cos(angles[0])
    ty = the_rxy*sin(angles[0])
    tmp = convert_coord(tx,ty, /data,/to_normal)
    tx = tmp[0]+xchsz*1
    ty = tmp[1]
    ww_center = ww_angle*mean(rxy_extent)
    ww_10re = ww_angle*10
    xyouts, tx,ty,/normal, $
        'W = '+string(ww_angle*deg,format='(I0)')+' deg', charsize=label_size


    ; Add SC locations.
    line_color = sgcolor('silver')
    psym = 8
    symsize = 0.5
    nangle = 20
    angles = smkarthm(0,2*!dpi,nangle,'n')
    circle_x = cos(angles)
    circle_y = sin(angles)
    usersym, circle_x, circle_y, /fill
    foreach probe, probes, probe_id do begin
        xx = reform(pos_coord[probe_id,1,0])
        yy = reform(pos_coord[probe_id,1,1])
        zz = reform(pos_coord[probe_id,1,2])
        color = probe_colors[probe_id]
        plots, xx, yy, psym=psym, symsize=symsize, color=color

        tmp = max(sqrt(xx^2+yy^2), index)
        tmp = min(xx, index)
        tmp = convert_coord(xx[index], yy[index], /data, /to_normal)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]+ychsz*0.0
        xyouts, tx,ty,/normal, strupcase(short_names[probe_id]), color=color, charsize=label_size
    endforeach


    ; Draw box.
    plot, xrange, yrange, /iso, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen, $
        nodata=1, noerase=1, position=tpos

    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*0.9
    xyouts, tx,ty, /normal, 'a-2. XY'





;---Middle panels.
    time_step = 3
    common_times = make_bins(mid_time_range, time_step)
    ntime = n_elements(common_times)
    nmid_probe = n_elements(mid_probes)

    ; Pflux and BBF.
    foreach probe, mid_probes, probe_id do begin
        prefix = probe+'_'

        ; Pflux.
        get_data, prefix+'pf_fac', times, pffac
        pffac = get_var_data(prefix+'pf_fac', at=common_times)
        cmap = get_var_data(prefix+'cmap', at=common_times)
        pf_para_norm = cmap*pffac[*,0]
        store_data, prefix+'pflux', common_times, pf_para_norm

        ; BBF.
        u_var = prefix+'u_gsm'
        if check_if_update(u_var, data_time_range) then begin
            the_probe = strmid(probe,2,1)
            themis_read_ion_vel, data_time_range, probe=the_probe
        endif
        u_gsm = get_var_data(u_var, at=common_times)
        r_gsm = get_var_data(prefix+'r_gsm', at=common_times)
        bbf = -vec_dot(u_gsm, sunitvec(r_gsm))
        store_data, prefix+'bbf', common_times, bbf
    endforeach


    ; injection.
    max_nebin = 4
    foreach probe, injection_probes do begin
        prefix = probe+'_'
        var = prefix+'kev_e_flux'
        get_data, var, times, fluxes, ebins, limits=lim
        index = lazy_where(ebins, '[]', injection_energy_range, count=nebin)
        fluxes = fluxes[*,index]
        ebins = ebins[index]
        if nebin gt max_nebin then begin
            fluxes = fluxes[*,-max_nebin:-1]
            ebins = ebins[-max_nebin:-1]
            ;fluxes = fluxes[*,0:*:2]
            ;ebins = ebins[0:*:2]
        endif
        index = lazy_where(times, '[]', mid_time_range)
        times = times[index]
        fluxes = fluxes[index,*]
        store_data, var+'_plot', times, fluxes, ebins
        yrange_log = alog10(minmax(fluxes))
        yrange = 10d^[floor(yrange_log[0]),ceil(yrange_log[1])]
        ;yrange = (var eq prefix+'kev_e_flux')? [1,1e7]: [1,10e5]
        ytickv_log = make_bins(alog10(yrange),1)
        ytickv = 10d^ytickv_log
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytickn = '10!U'+string(ytickv_log,format='(I0)')
        ytickn[0:*:2] = ' '
        add_setting, var+'_plot', /smart, {$
            display_type: 'list', $
            ylog:1, $
            yrange: yrange, $
            yticks: yticks, $
            ytickv: ytickv, $
            yminor: yminor, $
            ytickname: ytickn, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'e!U-!N flux'}
    endforeach


    ; Start to plot.
    xrange = mid_time_range
    xstep = 60*10
    xtickv = make_bins(xrange, xstep)
    xtickn = time_string(xtickv,tformat='hh:mm')
;    xtickn[0] += '!C'+time_string(xtickv[0],tformat='YYYY-MM-DD')
    xtickn[0] = time_string(xtickv[0],tformat='YYYY!CMTH DD')
    xticks = n_elements(xtickv)-1
    xminor = 6

    plot_vars = ['th'+['d','e','a'], injection_probes+'_kev_e_flux_plot']

    foreach probe, 'th'+['d','e','a'], probe_id do begin
        prefix = probe+'_'
        tpos = mid_poss[*,probe_id]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        ; Pflux.
        yys = get_var_data(prefix+'pflux', times=xxs)
        case probe of
            'thd': begin
                yrange = [-1,2]*50
                end
            'the': begin
                yrange = [-1,2]*30
                end
            'tha': begin
                yrange = [-1,2]*100
                end
        endcase
        yminor = 5
        yticks = 2
        ytitle = '(mW/m!U2!N)'
        xtickformat='(A1)'
        plot, xxs, yys, $
            xstyle=1, xrange=xrange, xlog=0, xtickv=xtickv, xtickname=xtickn, $
            xticks=xticks, xticklen=xticklen, xtickformat=xtickformat, xminor=xminor, $
            ystyle=9, yrange=yrange, ylog=0, $
            yticks=yticks, yticklen=yticklen, yminor=yminor, $
            position=tpos, noerase=1, nodata=0
        plots, xrange, [0,0], linestyle=1


        yys = get_var_data(prefix+'bbf', times=xxs)
        ytitle = '(km/s)'
        case probe of
            'thd': yrange = [-1,2]*200
            'the': yrange = [-1,2]*200
            'tha': yrange = [-1,2]*200
        endcase
        axis, yaxis=1, save=1, ystyle=1, yrange=yrange, ylog=0, $
            yticks=yticks, yticklen=yticklen, yminor=yminor, color=second_color
        plots, xxs, yys, color=second_color

        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*0.9
        xyouts, tx,ty, /normal, alignment=1, 'b-'+string(probe_id+1,format='(I0)')+'. '+strupcase(probe)
    endforeach

    ytitle = 'S!D||!N @100 km (mW/m!U2!N)'
    tx = tpos[0]-ychsz*1.2
    ty = (mid_poss[1,0]+mid_poss[3,2])*0.5
    xyouts, tx,ty,/normal, alignment=0.5, orientation=90, ytitle

    ytitle = 'Inward ion vel (km/s)'
    tx = tpos[2]+ychsz*1.2
    xyouts, tx,ty,/normal, alignment=0.5, orientation=90, ytitle, color=second_color



    label_size = 0.7
    foreach panel_id, [3,4] do begin
        tpos = mid_poss[*,panel_id]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        var = plot_vars[panel_id]
        get_data, var, xxs, yys, zzs, limits=lim

        yrange = lim.yrange
        ytickv = lim.ytickv
        yticks = lim.yticks
        ytickn = lim.ytickname
        yminor = 10

        if panel_id eq 3 then ytickn[-1] = ' '

        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, xlog=0, $
            ystyle=5, yrange=yrange, ylog=1, $
            position=tpos, noerase=1, nodata=1

        nbin = n_elements(zzs)
        dy = (tpos[3]-tpos[1]-label_size*ychsz)/(nbin-1)
        label_ys = tpos[3]-dy*(findgen(nbin))
        foreach color, lim.colors, color_id do begin
            plots, xxs, yys[*,color_id], color=color
            tx = tpos[2]+xchsz*1
            ty = label_ys[color_id]-ychsz*label_size*0.8
            xyouts, tx,ty,/normal, lim.labels[color_id], color=color, charsize=label_size
        endforeach

        xtickformat = (panel_id eq 4)? '': '(A1)'
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xlog=0, xtickv=xtickv, xtickname=xtickn, $
            xticks=xticks, xticklen=xticklen, xtickformat=xtickformat, xminor=xminor, $
            ystyle=1, yrange=yrange, ylog=1, ytickv=ytickv, ytickn=ytickn, $
            yticks=yticks, yticklen=yticklen, ytickformat=ytickformat, yminor=yminor, $
            position=tpos, noerase=1, nodata=1

        tx = tpos[2]-xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        the_probe = (panel_id eq 3)? 'G15': 'LANL'
        xyouts, tx,ty, /normal, alignment=1, 'b-'+string(panel_id+1,format='(I0)')+'. '+the_probe
    endforeach
    tx = tpos[0]-ychsz*1.5
    ty = (mid_poss[1,3]+mid_poss[3,4])*0.5
    xyouts, tx,ty,/normal, alignment=0.5, orientation=90, lim.ytitle


    tpos = mid_pan_pos
    yrange = [0,1]
    plot, xrange, yrange, $
        xstyle=5, ystyle=5, xrange=xrange, yrange=yrange, $
        position=tpos, noerase=1, nodata=1
    foreach tx, test_times do begin
        plots, tx+[0,0], yrange, linestyle=1
    endforeach



;---Right panel.
    asf_times = time_double('2014-08-28/'+$
        ['10:09','10:11','10:13'])
    ewo_time_range = time_double(['2014-08-28/10:09','2014-08-28/10:15'])
    mlat_range = [61,70]
    mlt_range = [-1.2,1]
    asf_color_table = 49
    zrange = [0,300]
    ewo_zrange = [0,150]
    ewo_color_table = 49
    ztitle = 'Brightness count (#)'
    ztickv = make_bins(zrange, 100)
    zticks = n_elements(ztickv)-1

    yrange = mlat_range
    ystep = 5
    ytickv = make_bins(yrange,ystep,/inner)
    yticks = n_elements(ytickv)-1
    yminor = ystep

    xrange = mlt_range
    xstep = 1
    xtickv = make_bins(xrange,xstep,/inner)
    xticks = n_elements(xtickv)-1
    xminor = 5
    xtitle = 'MLT (h)'


    themis_read_mltimg, ewo_time_range
    mltimg_var = 'thg_mltimg'
    get_data, mltimg_var, times, mltimg
    mlt_bins = get_setting(mltimg_var, 'mlt_bins')
    mlat_bins = get_setting(mltimg_var, 'mlat_bins')
    mlt_index = lazy_where(mlt_bins, '[]', mlt_range)
    mlat_index = lazy_where(mlat_bins, '[]', mlat_range)


    foreach panel_id, [0,1,2] do begin
        tpos = right_poss[*,panel_id]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

        index = where(times eq asf_times[panel_id])
        img = reform(mltimg[index[0],mlt_index,mlat_index])
        zzs = bytscl(img, min=zrange[0], max=zrange[1], top=254)

        sgtv, zzs, resize=1, position=tpos, ct=asf_color_table
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xticklen=xticklen, xminor=xminor, xtickformat='(A1)', $
            ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yticklen=yticklen, yminor=yminor, $
            position=tpos, nodata=1, noerase=1

        tx = tpos[2]-xchsz*0.5
        ty = tpos[1]+ychsz*0.2
        xyouts, tx,ty,/normal, alignment=1, time_string(asf_times[panel_id],tformat='hh:mm'), charsize=label_size
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.9
        xyouts, tx,ty,/normal, 'c-'+string(panel_id+1,format='(I0)')
    endforeach
    tx = tpos[0]-ychsz*0.8
    ty = (right_poss[1,0]+right_poss[3,2])*0.5
    xyouts, tx,ty,/normal, alignment=0.5, orientation=90, 'MLat (deg)'

    cbpos = right_poss[*,0]
    cbpos[1] = right_poss[1,nright_pan-2]
    cbpos[0] = right_poss[2,0]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.5
    sgcolorbar, zrange=zrange, ztitle=ztitle, ztickv=ztickv, zticks=zticks, ct=asf_color_table, position=cbpos


    ; EWOgram.
    ewo_var = 'thg_asf_ewo'
    if check_if_update(ewo_var, ewo_time_range) then $
        themis_read_mltimg_ewo, ewo_time_range, mlat_range=mlat_range
    get_data, ewo_var, times, ewo, mlt_bins
    index = lazy_where(mlt_bins, '[]', mlt_range)
    ewo = ewo[*,index]
    index = lazy_where(times, '[]', ewo_time_range)
    ewo = ewo[index,*]

    zzs = bytscl(transpose(ewo), min=ewo_zrange[0], max=ewo_zrange[1], top=254)
    zzs = reverse(zzs,2)

    tpos = right_poss[*,-1]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    sgtv, zzs, resize=1, position=tpos, ct=ewo_color_table


    ; Draw box.
    yrange = reverse(ewo_time_range)
    ystep = -60*3.
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    ytickn = time_string(ytickv,tformat='hh:mm')
    yminor = abs(ystep/60)
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xticks=xticks, xtickv=xtickv, xticklen=xticklen, xminor=xminor, xtitle=xtitle, $
        ystyle=1, yrange=yrange, yticks=yticks, ytickv=ytickv, yticklen=yticklen, yminor=yminor, ytickname=ytickn, $
        position=tpos, nodata=1, noerase=1
        
        
    ; Add lines.
    yys = dp_info.plot_time_range
    xxs = dp_info.plot_mlt_range[0]+[0,dp_info.omega*total(dp_info.plot_time_range*[-1,1])/15/60]
    dx = 5*dp_info.omega/15
    ;oplot, xxs-dx, yys, color=sgcolor('white')
    
    ; The auroral east-west lines.
    time0 = time_double('2014-08-28/10:12:30')
    mlt0 = 0.1
    mlts = [-1,1]
    line_color = sgcolor('black')
    times = time_double(['2014-08-28/10:13:30','2014-08-28/10:13:00'])
    foreach time, times, time_id do begin
        xxs = [mlt0,mlts[time_id]]
        yys = [time0,time]
        oplot, xxs, yys-30, color=line_color
        the_omega = -total(xxs*[-1,1])/total(yys*[-1,1])*15*60  ; eastward prop is +, westward prop is -.
        tx = mean(xxs)
        ty = time0-90
        alignment = 0.5
        xyouts, tx,ty,/data, alignment=alignment, $
            strtrim(string(the_omega,format='(F5.1)'),2)+'!Cdeg/min', color=line_color, charsize=label_size
    endforeach


    ewo_ztitle = 'Avg. brightness (#)'
    cbpos = tpos[*,0]
    cbpos[0] = tpos[2]+xchsz*0.8
    cbpos[2] = cbpos[0]+xchsz*0.5

    ztitle = 'Avg. count (#)'
    zrange = ewo_zrange
    zstep = 50
    ztickv = make_bins(zrange, zstep)
    zticks = n_elements(ztickv)-1
    sgcolorbar, zrange=zrange, ztitle=ztitle, ztickv=ztickv, zticks=zticks, ct=ewo_color_table, position=cbpos


    tx = tpos[2]-xchsz*0.5
    ty = tpos[1]+ychsz*0.2
    xyouts, tx,ty,/normal, alignment=1, 'c-4'

    if keyword_set(test) then stop
    sgclose

end
