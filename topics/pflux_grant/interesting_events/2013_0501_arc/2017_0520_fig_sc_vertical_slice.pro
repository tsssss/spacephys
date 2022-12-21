
test = 0

the_time = time_double('2017-05-20/01:35')
mission_probe = 'arase'

models = ['t89','t96','t01','t04s']
foreach model, models do begin
    igrf = 1
    t89_par = 2
    dir = 1
    xrange = [1.2,-10]
    yrange = [-2,4]

    line_mlats = [60,65,70,75,80]
    nline = n_elements(line_mlats)
    lines = list()
    ps = geopack_recalc(the_time)
    label_size = 0.8
    if n_elements(h0) eq 0 then h0 = 100d
    r0 = h0/constant('re')+1
    rad = constant('rad')
    deg = constant('deg')


    plot_file = join_path([srootdir(),time_string(the_time,tformat='YYYY_MMDD_fig_sc_vertical_slice_'+mission_probe+'_'+model+'_v01.pdf')])
    if keyword_set(test) then plot_file = 0

    ;---Get the position.
    time_range = the_time+[-1,1]*600
    pinfo = resolve_probe(mission_probe)
    routine = pinfo['routine_name']+'_read_orbit'
    probe = pinfo['probe']
    prefix = pinfo['prefix']
    r_gsm_var = call_function(routine, time_range, probe=probe, coord='gsm')

    ;---Get MLT.
    r_sm = transpose(get_var_data(r_gsm_var, at=the_time))
    the_mlt = pseudo_mlt(r_sm)
    line_mlts = the_mlt+dblarr(nline)

    par_var = geopack_read_par(time_range, model=model, t89_par=t89_par)
    par = get_var_data(par_var, at=the_time)

    t89 = (model eq 't89')? 1: 0
    t96 = (model eq 't96')? 1: 0
    t01 = (model eq 't01')? 1: 0
    t04s = (model eq 't04s')? 1: 0
    index = strpos(model, 's')
    storm = (index[0] ge 0)? 1: 0

    for i=0,nline-1 do begin
        tmp1 = (24-line_mlts[i])*15*rad  ; angle between mlt and midnight.
        tmp2 = line_mlats[i]*rad
        v_sm = [-cos(tmp2)*cos(tmp1),cos(tmp2)*sin(tmp1),sin(tmp2)]
        v0 = cotran(v_sm, the_time, 'sm2gsm')
        geopack_trace, v0[0],v0[1],v0[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+the_time, 'gsm2sm')
        srotate, fline, (24-the_mlt)*15*rad, 2
        lines.add, fline
    endfor

    step = 5
    if n_elements(xrange) eq 0 then begin
        xr = []
        foreach line, lines do xr = [xr, minmax(lines[*,0])]
        xr = minmax(make_bins(xr,step))
    endif else xr = xrange
    if n_elements(yrange) eq 0 then begin
        yr = []
        foreach line, lines do yr = [yr, minmax(lines[*,2])]
        yr = minmax(make_bins(yr,step))
    endif else yr = yrange

    pansize = abs([total(xr*[-1,1]),total(yr*[-1,1])])
    pansize = pansize/pansize[0]*4
    tpos = panel_pos(plot_file, fig_size=fig_size, pansize=pansize)
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    xticklen_chsz = -0.15   ; in ychsz.
    yticklen_chsz = -0.30   ; in xchsz.
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xminor = 2
    xtickv = make_bins(xr, xminor, inner=1)
    xticks = n_elements(xtickv)
    yminor = 2
    ytickv = make_bins(yr, yminor, inner=1)
    yticks = n_elements(ytickv)

    ; Add axis.
    plot, xr, yr, $
        xstyle=1, xrange=xr, xtitle='SM X (Re)', xminor=xminor, xticks=xticks, xtickv=xtickv, $
        ystyle=1, yrange=yr, ytitle='SM Y (Re)', yminor=yminor, yticks=yticks, ytickv=ytickv, $
        nodata=1, noerase=1, position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen

    ; Add SC.
    the_r_gsm = get_var_data(r_gsm_var, at=the_time)
    the_r_sm = transpose(cotran(the_r_gsm, the_time, 'gsm2sm'))
    srotate, the_r_sm, (24-the_mlt)*15*rad, 2
    plots, the_r_sm[0], the_r_sm[2], psym=1, color=sgcolor('red')
    tmp = convert_coord(the_r_sm[0],the_r_sm[2], data=1, to_normal=1)
    tx = tmp[0]+xchsz*0.5
    ty = tmp[1]
    msg = strupcase(mission_probe)
    xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sgcolor('red')

    foreach dir, [-1,1] do begin
        geopack_trace, the_r_gsm[0],the_r_gsm[1],the_r_gsm[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+the_time, 'gsm2sm')
        srotate, fline, (24-the_mlt)*15*rad, 2
        oplot, fline[*,0], fline[*,2], linestyle=1, color=sgcolor('red')
        ; Add Fpt.
        if dir eq -1 then begin
            f_gsm = [xf,yf,zf]
            f_mag = cotran(f_gsm, the_time, 'gsm2mag')
            fmlat = asin(f_mag[2]/r0)*deg
            tmp = convert_coord(fline[0,0],fline[0,2], data=1, to_normal=1)
            tx = tmp[0]+xchsz*0.5
            ty = tmp[1]-ychsz*2
            msg = 'F/MLat: '+strtrim(string(fmlat,format='(F4.1)'))+' deg'
            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sgcolor('red')
            ty = tmp[1]-ychsz*1
            msg = string(the_mlt,format='(F4.1)')
            if the_mlt le 0 then msg = string(the_mlt+24,format='(F4.1)')
            msg = 'MLT: '+strtrim(msg,2)+' h'
            xyouts, tx,ty,normal=1, msg, charsize=label_size, color=sgcolor('red')
        endif
    endforeach


    foreach line, lines do begin
        oplot, line[*,0], line[*,2]
    endforeach

    ; Add earth.
    tmp = smkarthm(0,2*!dpi, 40, 'n')
    txs = cos(tmp)
    tys = sin(tmp)
    polyfill, txs>0, tys, color=sgcolor('white')
    polyfill, txs<0, tys, color=sgcolor('grey')
    plots, txs, tys

    ; Add labeling.
    tx = tpos[0]+xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = time_string(the_time,tformat='YYYY-MM-DD/hh:mm')+' UT'
    xyouts, tx,ty,normal=1, msg, charsize=label_size
    ty = tpos[3]-ychsz*2
    msg = 'Model: '+strupcase(model)
    xyouts, tx,ty,normal=1, msg, charsize=label_size



    if keyword_set(test) then stop
    sgclose
endforeach

end