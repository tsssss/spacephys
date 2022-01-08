;+
; Test temporal evolution of the equatorial footpoint of certain field lines.
;-


    event_info = _2015_0218_02_load_data()
    test = 0


    probe = event_info['probe']
    prefix = event_info['prefix']
    iono_mlats = [63d,64]
    iono_mlts = [22.2d]
    time_range = time_double('2015-02-18/'+['02:05','02:15'])
    time_step = 60d
    times = make_bins(time_range, time_step)
    ntime = n_elements(times)
    ndim = 3
    rad = constant('rad')
    deg = constant('deg')

    dir = 1
    igrf = 1
    model = 't96'
    sgeopack_par, time_range, model
    
    xrange = [0,6]
    yrange = [-1,3]
    xtitle = 'SM Rxy (Re)'
    ytitle = 'SM Z (Re)'
    
    eq_rxy = fltarr(ntime)
    foreach time, times, time_id do begin
        tilt = geopack_recalc(time)
        par = get_var_data(model+'_par', at=time)
        
        sgopen, 0, xsize=5, ysize=4
        tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz)
        plot, xrange, yrange, $
            xstyle=1, xrange=xrange, xtitle=xtitle, $
            ystyle=1, yrange=yrange, ytitle=ytitle, $
            position=tpos, iso=1, nodata=1, noerase=1

        foreach mlat, iono_mlats, mlat_id do begin
            foreach mlt, iono_mlts, mlt_id do begin
                lat = mlat*rad
                lon = (mlt+12)*rad
                r_sm = transpose([cos(lon)*cos(lat),sin(lon)*cos(lat),sin(lat)])
                r_gsm = cotran(r_sm, time, 'sm2gsm')
                rx = r_gsm[0]
                ry = r_gsm[1]
                rz = r_gsm[2]
                geopack_trace, rx,ry,rz, dir, par, tilt=tilt, fx,fy,fz, igrf=igrf, t96=1, fline=fline
                fline = cotran(fline, time, 'gsm2sm')
                
                if mlat_id eq 0 and mlt_id eq 0 then begin
                    txs = snorm(fline[*,0:1])
                    tys = fline[*,2]
                    oplot, txs, tys
                    tmp = max(txs, index)
                    eq_rxy[time_id] = txs[index]
                endif
            endforeach
        endforeach

        oplot, xrange, [0,0], linestyle=1
        oplot, [0,0], yrange, linestyle=1
        tmp = smkarthm(0,2*!dpi,40,'n')
        circ_x = cos(tmp)
        circ_y = sin(tmp)
        oplot, circ_x, circ_y
        
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*2
        xyouts, tx,ty, /normal, time_string(time)

        sgclose
    endforeach

    plot, eq_rxy, ynozero=1
stop


end
