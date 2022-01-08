;+
; Plot traced field lines in a vertical cut.
;-

test = 0

    rad = constant('rad')

    model = 't01'
    the_time = time_double('2014-08-28/10:20')
    tilt = geopack_recalc(the_time)

    line_info = list()
    up_current = dictionary('lat_range', [60,67], 'color', sgcolor('salmon'))
    down_current = dictionary('lat_range', [68,72], 'color', sgcolor('sky_blue'))
    foreach lat, make_bins(up_current.lat_range,1)*rad do line_info.add, dictionary('lat',lat, 'color', up_current.color)
    foreach lat, make_bins(down_current.lat_range,1)*rad do line_info.add, dictionary('lat',lat, 'color', down_current.color)
    foreach lat, [50,60,70,80,90,-80,-90,100,110,120,130]*rad do line_info.add, dictionary('lat',lat, 'color', sgcolor('black'))
    ;foreach lat, [60,66,72]*rad do line_info.add, dictionary('lat',lat, 'color', sgcolor('red'))
    nline = line_info.length

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


    xrange = [5,-20]
    yrange = [0,1]*7
    xstep = 5
    xtickv = make_bins(xrange, xstep, /inner)
    xticks = n_elements(xtickv)-1
    xminor = xstep
    ystep = 5
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = ystep

    panel_xsize = 5
    panel_ysize = abs(panel_xsize*total(yrange*[-1,1])/total(xrange*[-1,1]))
    margins = [7,5,3,1]
    sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete
    fig_xsize = panel_xsize+total(margins[[0,2]])*abs_xchsz
    fig_ysize = panel_ysize+total(margins[[1,3]])*abs_ychsz

    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    ofn = join_path([root_dir,'case2_fig_field_line_trace_'+model+'_v2.pdf'])
    if keyword_set(test) then ofn = 0
    sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize
    tpos = sgcalcpos(1, margins=margins, xchsz=xchsz, ychsz=ychsz)
    xticklen_chsz = -0.30
    yticklen_chsz = -0.40
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])

    xtitle = 'SM X (Re)'
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
        xstyle=1, xrange=xrange, xticklen=xticklen, xtitle=xtitle, $
        xticks=xticks, xtickv=xtickv, xminor=xminor, $
        ystyle=1, yrange=yrange, yticklen=yticklen, ytitle=ytitle, $
        yticks=yticks, ytickv=ytickv, yminor=yminor, $
        position=tpos, /iso

    tx = tpos[2]-xchsz*0.5
    ty = tpos[3]-ychsz*1
    msg = 'Model ('+strupcase(model)+') magnetic field lines at local midnight'
    xyouts, tx,ty,alignment=1, /normal, msg, charsize=label_size

    yspace = 0.9
    label_size = 0.9
    msg = 'Downward Current!C '+strjoin(string(down_current.lat_range,'(I0)'),'-')+' deg'
    tx = -8
    ty = 3.5
    xyouts, tx,ty,alignment=0.5, /data, msg, color=down_current.color, charsize=label_size
    msg = 'Upward Current!C'+strjoin(string(up_current.lat_range,'(I0)'),'-')+' deg'
    tx = -5.5
    ty = 1.5
    xyouts, tx,ty,alignment=0.5, /data, msg, color=up_current.color, charsize=label_size
    msg = 'Proxy for high-beta plasma!C|B|<'+string(min_bmag,format='(I0)')+' nT, current closure?'
    tx = -15
    ty = 1.0
    xyouts, tx,ty,alignment=0.5, /data, msg, color=min_bmag_color, charsize=label_size

    ;msg = 'Model (T89) runtime: '+time_string(the_time)+' UT'
    ;ty = ty-ychsz*4*yspace
    ;xyouts, tx,ty,alignment=1, /normal, msg, charsize=label_size

    if keyword_set(test) then stop
    sgclose

end
