;+
; Load B data and calculate tilt angle and plot the position of spacecraft.
;-
;

function test_bield_tilt_angle_resolve_probe_type, probe0, the_probe=the_probe

    probe = strlowcase(probe0)
    missions = ['rbsp','th','g']
    full_name = dictionary()
    full_name.rbsp = 'rbsp'
    full_name.th = 'themis'
    full_name.g = 'goes'

    foreach mission, missions do begin
        length = strlen(mission)
        if strmid(probe,0,length) eq mission then begin
            the_probe = strmid(probe,length)
            return, full_name[mission]
        endif
    endforeach

end


pro test_bfield_tilt_angle, probes=probes, time_range=time_range, plot_bmag=plot_bmag, event=event

    deg = 180d/!dpi
    trace_dir = -1
    model = 't89'
    probes = event['probes']
    time_range = event['time_range']
    id = event['id']

    test = 0


    foreach probe, probes do begin
        mission = test_bield_tilt_angle_resolve_probe_type(probe, the_probe=the_probe)

        ; Load xxx_r_gsm, xxx_b_gsm.
        call_procedure, mission+'_read_orbit', time_range, probe=the_probe
        call_procedure, mission+'_read_bfield', time_range, probe=the_probe

        ; Convert xxx_b_gsm to xxx_b_sm and get tilt angle.
        bgsm_var = probe+'_b_gsm'
        if tnames(bgsm_var) eq '' then continue
        get_data, bgsm_var, times, bgsm
        bsm = cotran(bgsm, times, 'gsm2sm')
        tilt = atan(bsm[*,2],sqrt(bsm[*,0]^2+bsm[*,1]^2))*deg
        store_data, probe+'_b_tilt', times, tilt, limits={ytitle:'(deg)',labels:strupcase(probe)}


        ; Convert xxx_r_gsm to xxx_r_mag and get MLT.
        rgsm_var = probe+'_r_gsm'
        if tnames(rgsm_var) eq '' then continue
        get_data, rgsm_var, times, rgsm
        rmag = cotran(rgsm, times, 'gsm2mag')
        mlat = asin(rmag[*,2]/snorm(rmag))*deg
        mlon = atan(rmag[*,1],rmag[*,0])*deg
        mlt = mlon2mlt(mlon, times)

        pre0 = probe+'_'
        mlat_var = pre0+'mlat'
        store_data, mlat_var, times, mlat
        add_setting, mlat_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLat'}

        mlon_var = pre0+'mlon'
        store_data, mlon_var, times, mlon
        add_setting, mlon_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLon'}

        mlt_var = pre0+'mlt'
        store_data, mlt_var, times, mlt
        add_setting, mlt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'MLT'}


        ; Calculates |B|-|B|_model.
        ; Get model field.
        read_geopack_info, probe+'_r_gsm', model=model, direction=trace_dir

        bgsm_var = probe+'_b_gsm'
        get_data, bgsm_var, times, bgsm
        bmag = snorm(bgsm)

        bmodmag = snorm(get_var_data(probe+'_bmod_gsm_'+model, at=times))
        dbmag = bmag-bmodmag
        store_data, probe+'_db_mag', times, dbmag, limits={ytitle:'(nT)',labels:strupcase(probe)}
    endforeach


;---Filter out the probes that do not have data.
    nprobe = n_elements(probes)
    flags = bytarr(nprobe)
    for ii=0, nprobe-1 do begin
        if tnames(probes[ii]+'_b_gsm') eq '' then flags[ii] = 1
        if tnames(probes[ii]+'_r_gsm') eq '' then flags[ii] = 1
    endfor
    index = where(flags eq 0, nprobe)
    if nprobe eq 0 then return
    probes = probes[index]


;---Prepare to plot.
    fig_xsize = 4   ; inch.
    lmargin = 8
    rmargin = 6
    ticklen = -0.01
    sgopen, 0, xsize=fig_xsize, ysize=fig_xsize, /inch
    tpos = sgcalcpos(1,lmargin=lmargin, rmargin=rmargin, xchsz=xchsz, ychsz=ychsz)
    ypad = 0.5  ; y-pad space between panels.
    sgclose, /wdelete

    color_start = 50
    color_end = 250
    color_table = 40
    colors = round(smkarthm(color_start,color_end, nprobe, 'n'))
    for ii=0, nprobe-1 do colors[ii] = sgcolor(colors[ii], ct=color_table)

    ; TIlt angle panels.
    tilt_panel = dictionary()
    tilt_panel.aspect_ratio = 0.2
    tilt_panel_xsize = tpos[2]-tpos[0]
    tilt_panel_ysize = tilt_panel_xsize*tilt_panel.aspect_ratio*nprobe+ychsz*ypad*(nprobe-1)

    if keyword_set(plot_bmag) then begin
        panel_vars = probes+'_db_mag'
        panel_title = '|B|-|B|_T89'
    endif else begin
        panel_vars = probes+'_b_tilt'
        panel_title = 'Tilt angle in SM'
    endelse


    yticks = 2
    yminor = 5
    tilt_var_labels = strarr(nprobe)
    foreach tvar, panel_vars, ii do begin
        data = get_var_data(tvar)
        index = lazy_where(data-mean(data), [-1,1]*5*stddev(data))
        yrange = minmax(data[index])
        yrange = yrange-(yrange mod yticks)+[0,2]
        options, tvar, 'yrange', yrange
        options, tvar, 'yticks', yticks
        options, tvar, 'yminor', yminor
        options, tvar, 'colors', sgcolor('black')
        options, tvar, 'yticklen', ticklen
        options, tvar, 'xticklen', ticklen/tilt_panel.aspect_ratio
        options, tvar, 'labflag', -1
        options, tvar, 'ystyle', 1
        tilt_var_labels[ii] = get_setting(tvar, 'labels')
        if tilt_var_labels[ii] eq '' then tilt_var_labels[ii] = strupcase(strmid(tvar,0,strpos(tvar,'_')))
        options, tvar, 'labels', ''
    endforeach



    ; Orbit panel.
    pos_vars = probes+'_r_gsm'
    orbit_x_range = [-1,1]
    foreach tvar, pos_vars do orbit_x_range = minmax([orbit_x_range,minmax((get_var_data(tvar))[*,0])])
    orbit_y_range = [-1,1]
    foreach tvar, pos_vars do orbit_y_range = minmax([orbit_y_range,minmax((get_var_data(tvar))[*,1])])

    orbit_panel = dictionary()
    orbit_panel.xminor = 5
    orbit_panel.xrange = [ceil(max(orbit_x_range)),floor(min(orbit_x_range))]
    orbit_panel.xtickv = make_bins(orbit_panel.xrange, orbit_panel.xminor, /inner)
    orbit_panel.xticks = n_elements(orbit_panel.xtickv)-1
    orbit_panel.xtitle = 'GSM X (Re)'
    orbit_panel.yminor = 5
    orbit_panel.yrange = [ceil(max(orbit_y_range)),floor(min(orbit_y_range))]
    orbit_panel.ytickv = make_bins(orbit_panel.yrange, orbit_panel.yminor, /inner)
    orbit_panel.yticks = n_elements(orbit_panel.ytickv)-1
    orbit_panel.ytitle = 'GSM Y (Re)'
    orbit_panel.aspect_ratio = abs(double(orbit_panel.yrange[1]-orbit_panel.yrange[0])/(orbit_panel.xrange[1]-orbit_panel.xrange[0]))

    orbit_panel_xsize = tpos[2]-tpos[0]
    orbit_panel_ysize = orbit_panel_xsize*orbit_panel.aspect_ratio
    orbit_panel.xticklen = ticklen/orbit_panel.aspect_ratio
    orbit_panel.yticklen = ticklen

    tmargin = 2
    bmargin = 5
    fig_ysize = (tilt_panel_ysize+bmargin*2*ychsz+orbit_panel_ysize+tmargin*ychsz)*fig_xsize


;---Plot.
    if keyword_set(test) then begin
        file = test
        magnify = 1.2
    endif else begin
        file = (keyword_set(plot_bmag))? join_path([shomedir(),id+'_bmag.pdf']): $
            join_path([shomedir(),id+'_tilt_angle.pdf'])
        magnify = 1
    endelse

    sgopen, file, xsize=fig_xsize, ysize=fig_ysize, magnify=magnify
    poss = sgcalcpos(2, lmargin=lmargin, rmargin=rmargin, tmargin=tmargin, bmargin=bmargin, ypans=[tilt_panel_ysize,orbit_panel_ysize], ypad=bmargin)

    tilt_panel.pos = poss[*,0]
    orbit_panel.pos = poss[*,1]

    middle_time = mean(time_range)
    psym = 6
    label_size = 0.8

    ; The tilt angle.
    poss = sgcalcpos(nprobe, position=tilt_panel.pos, ypad=ypad)
    ; Sort by MLT.
    mlts = fltarr(nprobe)
    for ii=0, nprobe-1 do mlts[ii] = get_var_data(probes[ii]+'_mlt',at=middle_time)
    index = sort(mlts)
    panel_vars = panel_vars[index]
    tilt_var_labels = tilt_var_labels[index]
    probes = probes[index]

    ;tplot, panel_vars, trange=time_range, position=poss, /novtitle
    for ii=0, nprobe-1 do begin
        tpos = poss[*,ii]
        nouttick = (ii ne nprobe-1)? 1: 0
        tplot, panel_vars[ii], trange=time_range, position=tpos, /novtitle, nouttick=nouttick, /noerase
        tx = tpos[2]+xchsz*1*label_size
        ty = 0.5*(tpos[3]+tpos[1])-ychsz*0.1*label_size
        xyouts, tx,ty,/normal, tilt_var_labels[ii], color=colors[ii], charsize=label_size
    endfor
    tpos = poss[*,0]
    tx = tpos[0]
    ty = tpos[3]+ychsz*0.3*label_size
    xyouts, tx,ty,/normal, panel_title

    ; The positions.
    plot, orbit_panel.xrange, orbit_panel.yrange, $
        xstyle=5, ystyle=5, $
        /nodata, /noerase, position=orbit_panel.pos, _extra=orbit_panel.tostruct()
    for ii=0, nprobe-1 do begin
        rgsm = get_var_data(probes[ii]+'_r_gsm')
        plots, rgsm[*,0], rgsm[*,1], color=colors[ii], linestyle=1
        rgsm = get_var_data(probes[ii]+'_r_gsm', at=middle_time)
        plots, rgsm[0], rgsm[1], color=colors[ii], psym=psym, symsize=label_size
        tmp = convert_coord(rgsm[0:1], /data, /to_normal)
        tx = tmp[0]+xchsz*1*label_size
        ty = tmp[1]
        xyouts, tx,ty,/normal, tilt_var_labels[ii], color=colors[ii], charsize=label_size
    endfor

    ; Add earth and lines.
    tmp = 50
    tmp = findgen(tmp)/(tmp-1)*2*!dpi
    xs = cos(tmp)
    ys = sin(tmp)
    polyfill, xs<0, ys, /line_fill, orientation=45
    plots, xs, ys
    foreach r, [5,10] do oplot, xs*r, ys*r, linestyle=1

    plots, orbit_panel.xrange, [0,0], linestyle=1
    plots, [0,0], orbit_panel.yrange, linestyle=1

    plot, orbit_panel.xrange, orbit_panel.yrange, $
        xstyle=1, ystyle=1, $
        /nodata, /noerase, position=orbit_panel.pos, _extra=orbit_panel.tostruct()

    if keyword_set(test) then stop
    sgclose

end

probes = ['rbsp'+['b'], $
    'g'+['13','15'], $
    'th'+['a','d','e']]
time_range = time_double(['2014-08-28/09:50','2014-08-28/11:00'])

probes = ['rbsp'+['a'], $
    'g'+['13','15'], $
    'th'+['a','d','e']]
time_range = time_double(['2014-08-27/09:20','2014-08-27/10:30'])

probes = ['rbsp'+['a'], $
    'g'+['13','15'], $
    'th'+['a','d','e']]
time_range = time_double(['2014-12-22/01:30','2014-12-22/03:00'])


probes = ['rbsp'+['a','b'], $
    'g'+['13','15'], $
    'th'+['a','d','e']]
time_range = time_double(['2012-12-22/01:30','2014-12-22/03:30'])




; This is the checked events, including probes and time range.
event_list = hash()

    event_list['2014_0828_09'] = dictionary($
        'time_range', time_double(['2014-08-28/09:50','2014-08-28/11:00']), $
        'probes', ['rbspb','g13','g15','tha','thd','the'])
    event_list['2016_1013_12'] = dictionary($
        'time_range', time_double(['2016-10-13/12:00','2016-10-13/13:30']), $
        'probes', ['rbspa','rbspb','g13','g14','g15','thd'])
    event_list['2017_0328_09'] = dictionary($
        'time_range', time_double(['2017-03-28/09:00','2017-03-28/11:00']), $
        'probes', ['rbspa','rbspb','tha','thd','the'])
    event_list['2017_0301_22'] = dictionary($
        'time_range', time_double(['2017-03-01/21:00','2017-03-02/00:30']), $
        'probes', ['rbspa','rbspb','g13','g15','tha','the'])
    event_list['2013_1015_04'] = dictionary($
        'time_range', time_double(['2013-10-15/05:30','2013-10-15/07:30']), $
        'probes', ['rbspa','rbspb','g15','tha','thd','the'])


    
;event_list['2016_0502'] = dictionary($
;    'time_range', time_double(['2016-05-02/23:30','2016-05-03/05:00']), $
;    'probes', ['g13','g14','g15'])
;event_list['2016_0509'] = dictionary($
;    'time_range', time_double(['2016-05-10/01:00','2016-05-10/07:00']), $
;    'probes', ['g13','g14','g15'])
;event_list['2016_0605'] = dictionary($
;    'time_range', time_double(['2016-06-05/23:30','2016-06-06/02:00']), $
;    'probes', ['g13','g14','g15','rbspa','rbspb'])
;event_list['2016_0902'] = dictionary($
;    'time_range', time_double(['2016-09-02/03:30','2016-09-02/06:30']), $
;    'probes', ['g13','g14','g15'])
;event_list['2016_0929'] = dictionary($
;    'time_range', time_double(['2016-09-28/09:00','2016-09-28/11:30']), $
;    'probes', ['rbspb','g13','g14','g15'])
;event_list['2016_1013'] = dictionary($
;    'time_range', time_double(['2016-10-13/12:00','2016-10-13/13:30']), $
;    'probes', ['rbspa','rbspb','g13','g14','g15','thd'])
;event_list['2017_0302'] = dictionary($
;    'time_range', time_double(['2017-03-01/22:00','2017-03-01/23:30']), $
;    'probes', ['rbspa','rbspb','g13','the'])


;    event_list['2012_1001_02'] = dictionary($
;        'time_range', time_double(['2012-10-01/01:30','2012-10-01/03:30']), $
;        'probes', ['g'+['13','15'], 'th'+['a','d','e']])
;    event_list['2013_0321_07'] = dictionary($
;        'time_range', time_double(['2013-03-21/06:00','2013-03-21/08:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a','b'], 'th'+['a']])
;    event_list['2013_0518_10'] = dictionary($
;        'time_range', time_double(['2013-05-18/09:00','2013-05-18/11:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a'], 'th'+['a','d','e']])
;    event_list['2013_0601_02'] = dictionary($
;        'time_range', time_double(['2013-06-01/01:00','2013-06-01/04:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a','b'], 'th'+['a','d','e']])
;    event_list['2013_0607_05'] = dictionary($
;        'time_range', time_double(['2013-06-07/04:00','2013-06-07/07:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a','b'], 'th'+['d','e']])
;    event_list['2013_0710_03'] = dictionary($
;        'time_range', time_double(['2013-07-10/02:00','2013-07-10/04:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a'], 'th'+['a','d','e']])
;    event_list['2013_0710_05'] = dictionary($
;        'time_range', time_double(['2013-07-10/05:30','2013-07-10/08:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['b'], 'th'+['a','e']])
;    event_list['2013_0816_03'] = dictionary($
;        'time_range', time_double(['2013-08-16/02:00','2013-08-16/04:30']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a','b'], 'th'+['a','d','e']])
;    event_list['2013_1002_07'] = dictionary($
;        'time_range', time_double(['2013-10-02/04:00','2013-10-02/07:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a','b'], 'th'+['a','d','e']])
;    event_list['2014_1028_06'] = dictionary($
;        'time_range', time_double(['2014-10-28/06:00','2014-10-28/07:00']), $
;        'probes', ['g'+['13','15'], 'th'+['a','d','e']])
;    event_list['2014_1028_07'] = dictionary($
;        'time_range', time_double(['2014-10-28/07:00','2014-10-28/08:00']), $
;        'probes', ['g'+['15'], 'th'+['a','d','e']])
    foreach key, event_list.keys() do event_list[key].id = key

; This is for the shock event.
;    event_list['2014_0912_16'] = dictionary($
;        'time_range', time_double(['2014-09-12/15:50','2014-09-12/16:00']), $
;        'probes', ['g'+['13','15'], 'rbsp'+['a','b'], 'th'+['a','d','e']])

foreach event, event_list do $
    foreach plot_bmag, [0] do test_bfield_tilt_angle, event=event, plot_bmag=plot_bmag
end
