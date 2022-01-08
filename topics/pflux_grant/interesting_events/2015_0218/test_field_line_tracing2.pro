;+
; Test the equatorial footpoint of certain field lines for different models.
;-


    event_info = _2015_0218_02_load_data()
    test = 0


    probe = event_info['probe']
    prefix = event_info['prefix']
    iono_mlats = [63.2d]
    iono_mlts = [22.1d]
    time_range = time_double('2015-02-18/'+['02:05','02:15'])
    time_step = 60d
    times = make_bins(time_range, time_step)
    ntime = n_elements(times)
    ndim = 3
    rad = constant('rad')
    deg = constant('deg')

    dir = 1
    mlat = iono_mlats[0]
    mlt = iono_mlts[0]
    lat = mlat*rad
    lon = (mlt+12)*rad
    r_sm = transpose([cos(lon)*cos(lat),sin(lon)*cos(lat),sin(lat)])
    r_gsm = cotran(r_sm, time, 'sm2gsm')
    rx = r_gsm[0]
    ry = r_gsm[1]
    rz = r_gsm[2]

    models = ['t89','t96','t01']
    colors = sgcolor(['red','green','blue'])

    xrange = [0,7]
    yrange = [-1,1]*2.5
    xtitle = 'SM Rxy (Re)'
    ytitle = 'SM Z (Re)'
    xticklen = -0.02
    yticklen = -0.02

    xsize = abs(total(xrange*[-1,1]))
    ysize = abs(total(yrange*[-1,1]))
    tpos = panel_pos(0, 'scale', pansize=5*[1,ysize/xsize], fig_size=figsz)

    sgopen, 0, xsize=figsz[0], ysize=figsz[1], xchsz=xchsz, ychsz=ychsz
    plot, xrange, yrange, $
        xstyle=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
        ystyle=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, $
        position=tpos, iso=1, nodata=1, noerase=1


    time = time_double('2015-02-18/02:10')
    tilt = geopack_recalc(time)
    foreach model, models, model_id do begin
        sgeopack_par, time_range, model
        par = get_var_data(model+'_par', at=time)

        t89 = 0
        t96 = 0
        t01 = 0
        case model of
            't89': t89 = 1
            't96': t96 = 1
            't01': t01 = 1
        endcase

        foreach igrf, [0,1], igrf_id do begin
            geopack_trace, rx,ry,rz, dir, par, tilt=tilt, fx,fy,fz, igrf=igrf, fline=fline, $
                t89=t89, t96=t96, t01=t01
            fline = cotran(fline, time, 'gsm2sm')

            txs = snorm(fline[*,0:1])
            tys = fline[*,2]
            oplot, txs,tys, color=colors[model_id], linestyle=igrf_id

            print, (['Dipole','IGRF'])[igrf_id]
            print, atan(tys[0],txs[0])*deg
            print, atan(tys[-1],txs[-1])*deg
        endforeach
    endforeach

    ; Pure dipole.
;    iono_mlats = [63.7d]
;    iono_mlts = [22.1d]
;    time_range = time_double('2015-02-18/'+['02:05','02:15'])
;    time_step = 60d
;    times = make_bins(time_range, time_step)
;    ntime = n_elements(times)
;    ndim = 3
;    rad = constant('rad')
;    deg = constant('deg')
;
;    dir = 1
;    mlat = iono_mlats[0]
;    mlt = iono_mlts[0]
;    lat = mlat*rad
;    lon = (mlt+12)*rad
;    r_sm = transpose([cos(lon)*cos(lat),sin(lon)*cos(lat),sin(lat)])
;    r_gsm = cotran(r_sm, time, 'sm2gsm')
;    rx = r_gsm[0]
;    ry = r_gsm[1]
;    rz = r_gsm[2]
    
    foreach igrf, [0,1], igrf_id do begin
        geopack_trace, rx,ry,rz, dir, 0, tilt=tilt, fx,fy,fz, igrf=igrf_id, fline=fline, $
            t89=0, t96=0, t01=0
        fline = cotran(fline, time, 'gsm2sm')

        txs = snorm(fline[*,0:1])
        tys = fline[*,2]
        oplot, txs,tys, color=sgcolor('black'), linestyle=igrf_id
    endforeach


    ; Print models.
    foreach model, models, model_id do begin
        tx = tpos[0]+xchsz*1+model_id*5*xchsz
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,/normal, strupcase(model), color=colors[model_id]
    endforeach

    foreach type, ['Dipole','IGRF'], type_id do begin
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*(type_id+2)
        plots, tx+[0,2]*xchsz, ty+[0,0]+ychsz*0.3, normal=1, linestyle=type_id
        xyouts, tx+xchsz*3,ty,/normal, type
    endforeach


    oplot, xrange, [0,0], linestyle=1
    oplot, [0,0], yrange, linestyle=1
    tmp = smkarthm(0,2*!dpi,40,'n')
    circ_x = cos(tmp)
    circ_y = sin(tmp)
    oplot, circ_x, circ_y

    sgclose

stop


end
