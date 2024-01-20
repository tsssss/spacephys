;+
; Figure to show the vertical B field.
;-


function fig_vertical_config_v01, event_info=event_info


    sgopen, 0, size=[6,6]
    my_pos = sgcalcpos(1)

;---Input settings.
    pos_xrange = [7,-8]   ; actually Rxy.
    pos_yrange = [-1,1]*4   ; actually z.
    config_time = time_double('2015-03-17/12:00')
    config_time_range = time_double(['2015-03-17','2015-03-19'])
    tpos = my_pos
    xrange = pos_xrange
    yrange = pos_yrange



;---Model settings.
    external_model = 't04s'
    ;external_model = 't89'
    internal_model = 'dipole'
    dir = 1

    par_var = geopack_read_par(config_time_range, model=external_model, t89_par=t89_par)
    if external_model eq 't89' then par = 2 else par = get_var_data(par_var, at=config_time)

    model_info = geopack_resolve_model(external_model)
    t89 = model_info.t89
    t96 = model_info.t96
    t01 = model_info.t01
    t04s = model_info.ts04
    storm = model_info.storm


;---Constants.
    rad = constant('rad')
    deg = constant('deg')


;---The field lines.
    fline_mlt = 0d
    fline_mlats = smkarthm(40,72,2,'dx')
    fline_mlats = [45,50,55,60,65]
    nfline = n_elements(fline_mlats)*2
    ndim = 3
    f_gsms = fltarr(nfline,ndim)
    for ii=0,nfline*0.5-1 do begin
        fline_pp = fline_mlt*15*rad     ; phi
        fline_tt = fline_mlats[ii]*rad  ; theta
        the_f_sm = transpose([cos(fline_tt)*cos(fline_pp),cos(fline_tt)*sin(fline_pp),sin(fline_tt)])
        f_gsms[ii,*] = cotran_pro(the_f_sm, config_time, 'sm2gsm')
        fline_pp += !dpi
        the_f_sm = transpose([cos(fline_tt)*cos(fline_pp),cos(fline_tt)*sin(fline_pp),sin(fline_tt)])
        f_gsms[ii+nfline*0.5,*] = cotran_pro(the_f_sm, config_time, 'sm2gsm')
    endfor

    flines = list()
    for ii=0,nfline-1 do begin
        v0 = f_gsms[ii,*]
        geopack_trace, v0[0],v0[1],v0[2], dir, par, xf,yf,zf, fline=fline, $
            t89=t89, t96=t96, t01=t01, ts04=t04s, storm=storm, igrf=igrf, r0=r0
        ntime = n_elements(fline)/3
        fline = cotran(fline, fltarr(ntime)+config_time, 'gsm2sm')
        srotate, fline, (24-fline_mlt)*15*rad, 2
        flines.add, fline
    endfor


    
;---Plot.
    step = 5
    if n_elements(xrange) eq 0 then begin
        xr = []
        foreach line, flines do xr = [xr, minmax(line[*,0])]
        xr = minmax(make_bins(xr,step))
    endif else xr = xrange
    if n_elements(yrange) eq 0 then begin
        yr = []
        foreach line, flines do yr = [yr, minmax(line[*,2])]
        yr = minmax(make_bins(yr,step))
    endif else yr = yrange

    fig_size = double([!d.x_size,!d.y_size])
    tmp = get_charsize()
    xchsz = tmp[0]
    ychsz = tmp[1]
    uniform_ticklen = -ychsz*0.15*fig_size[1]
    
    xticklen = uniform_ticklen/(tpos[3]-tpos[1])/fig_size[1]
    yticklen = uniform_ticklen/(tpos[2]-tpos[0])/fig_size[0]

    xminor = 2
    xtickv = make_bins(xr, xminor, inner=1)
    xticks = n_elements(xtickv)
    xtickn = string(xtickv,format='(I0)')
    ;    xtitle = 'SM R!DXY!N (Re)'
    xtickn[-1] = 'SM R!DXY!N (Re)        '
    xtitle = ' '
    yminor = 2
    ytickv = make_bins(yr, yminor, inner=1)
    yticks = n_elements(ytickv)
    ytitle = 'SM Z (Re)'
    
    xstep = 1
    xgrids = make_bins(xr, xstep, inner=1)
    ystep = 1
    ygrids = make_bins(yr, ystep, inner=1)

    ; Add axis.
    plot, xr, yr, $
        xstyle=1, xrange=xr, xtitle=xtitle, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, $
        ystyle=1, yrange=yr, ytitle=ytitle, yminor=yminor, yticks=yticks, ytickv=ytickv, $
        nodata=1, noerase=1, position=tpos, iso=1, $
        xticklen=xticklen, yticklen=yticklen
    foreach val, xgrids do begin
        oplot, val+[0,0], yr, linestyle=1
    endforeach
    foreach val, ygrids do begin
        oplot, xr, val+[0,0], linestyle=1
    endforeach

    foreach fline, flines do begin
        oplot, fline[*,0], fline[*,2]
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
;    msg = time_string(the_time,tformat='YYYY-MM-DD/hh:mm')+' UT'
;    xyouts, tx,ty,normal=1, msg, charsize=label_size
;    ty = tpos[3]-ychsz*2
    msg = 'Model: '+strupcase(external_model)
    xyouts, tx,ty,normal=1, msg, charsize=label_size

stop

end

print, fig_vertical_config_v01()
end