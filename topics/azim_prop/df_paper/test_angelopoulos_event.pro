;+
; Plot in the time-MLT plane, and color-code the tilt angle data, for a given time.
; Radial distance range is 4 to 12 Re.
;-

function data2color, data, range=range, ct=ct, reverse_ct=reverse_ct, index_color=index_color

    if n_elements(range) ne 2 then begin
        range = [-1,1]*90
    endif
    ;the_range = range
    ;index_color = bytscl(data, min=the_range[0], max=the_range[1])


;    ; make data logarithmic.
    linear_range = 1
    signs = sign(data)
    the_data = alog10(abs(data)>linear_range)*signs
    the_range = [-1,1]*alog10(max(abs(range)))
    index_color = bytscl(the_data, min=the_range[0], max=the_range[1])


;    ; tanh.
;    the_range = [-1,1]*90
;    index_color = bytscl(tanh(data/30), min=tanh(-3), max=tanh(3))


    if keyword_set(reverse_ct) then index_color = 256-index_color

    if n_elements(ct) eq 0 then ct = 70
    loadct, ct, rgb_table=rgb0
    rgb0 = float(rgb0)
    true_colors = 256L*(256L*rgb0[index_color,2]+rgb0[index_color,1])+rgb0[index_color,0]

    return, true_colors
end

pro azim_prop_gen_storm_keox_ewo, time_range, $
    dis_range=dis_range, mlt_range=mlt_range

    test = 0

;---Preparation.
    ; Cut off data outside the range.
    if n_elements(dis_range) ne 2 then dis_range = [5,20]
    if n_elements(mlt_range) ne 2 then mlt_range = [-1,1]*9
    if n_elements(xsm_range) ne 2 then xsm_range = [-20,5]

    ; All spacecraft.
    mission_info = dictionary()
    mission_info['rbspa'] = dictionary($
        'routine_name', 'rbsp', $
        'prefix', 'rbspa_', $
        'probe', 'a', $
        'short_name', 'rba')
    mission_info['rbspb'] = dictionary($
        'routine_name', 'rbsp', $
        'prefix', 'rbspb_', $
        'probe', 'b', $
        'short_name', 'rbb')
    mission_info['g13'] = dictionary($
        'routine_name', 'goes', $
        'prefix', 'g13_', $
        'probe', '13', $
        'short_name', 'g13')
    mission_info['g14'] = dictionary($
        'routine_name', 'goes', $
        'prefix', 'g14_', $
        'probe', '14', $
        'short_name', 'g14')
    mission_info['g15'] = dictionary($
        'routine_name', 'goes', $
        'prefix', 'g15_', $
        'probe', '15', $
        'short_name', 'g15')
    mission_info['tha'] = dictionary($
        'routine_name', 'themis', $
        'prefix', 'tha_', $
        'probe', 'a', $
        'short_name', 'tha')
;    mission_info['thc'] = dictionary($
;        'routine_name', 'themis', $
;        'prefix', 'thc_', $
;        'probe', 'c', $
;        'short_name', 'thc')
;    mission_info['thb'] = dictionary($
;        'routine_name', 'themis', $
;        'prefix', 'thb_', $
;        'probe', 'b', $
;        'short_name', 'thb')
    mission_info['thd'] = dictionary($
        'routine_name', 'themis', $
        'prefix', 'thd_', $
        'probe', 'd', $
        'short_name', 'thd')
    mission_info['the'] = dictionary($
        'routine_name', 'themis', $
        'prefix', 'the_', $
        'probe', 'e', $
        'short_name', 'the')
    mission_probes = mission_info.keys()
    mission_probes = mission_probes[sort(mission_probes)]


;---Load data.
    ; Clear memory.
    del_data, '*'
    have_data_flags = bytarr(n_elements(mission_probes))
    foreach mission_probe, mission_probes, ii do begin
        routine_name = mission_info[mission_probe].routine_name
        probe = mission_info[mission_probe].probe
        call_procedure, routine_name+'_read_bfield', time_range, probe=probe
        call_procedure, routine_name+'_read_orbit', time_range, probe=probe
        flag = 1
        prefix = mission_info[mission_probe].prefix
        if tnames(prefix+'b_gsm') eq '' then flag = 0
        if tnames(prefix+'r_gsm') eq '' then flag = 0
        have_data_flags[ii] = flag
    endforeach

    index = where(have_data_flags eq 1, nmission_probe)
    if nmission_probe eq 0 then begin
        errmsg = handle_error('No spacecraft has data ...')
        return
    endif
    mission_probes = mission_probes[index]


;---Calculate derived data.
    deg = 180d/!dpi
    ndim = 3
    par = 2
    time_step = 60.
    foreach mission_probe, mission_probes do begin
        prefix = mission_info[mission_probe].prefix
        r_var = prefix+'r_gsm'
        get_data, r_var, times, rgsm, limits=lim

        ; Calculate MLT.
        rmag = cotran(rgsm, times, 'gsm2mag')
        mlon = atan(rmag[*,1],rmag[*,0])*deg
        mlt = mlon2mlt(mlon, times)
        mlt_var = prefix+'mlt'
        store_data, mlt_var, times, mlt
        add_setting, mlt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'hr', $
            short_name: 'MLT' }

        ; Calculate distance.
        dis = snorm(rgsm)
        dis_var = prefix+'dis'
        store_data, dis_var, times, dis
        add_setting, dis_var, /smart, {$
            display_type: 'scalar', $
            unit: 'Re', $
            short_name: 'R'}

        ; Calculate R SM.
        rsm = cotran(rgsm, times, 'gsm2sm')
        rsm_var = prefix+'r_sm'
        store_data, rsm_var, times, rsm
        add_setting, rsm_var, /smart, {$
            display_type: 'vector', $
            unit: 'Re', $
            short_name: 'R', $
            coord: 'SM', $
            coord_labels: ['x','y','z']}

        ; Calculate B model.
        ntime = n_elements(times)
        b0gsm = fltarr(ntime,ndim)
        for ii=0, ntime-1 do begin
            tilt = geopack_recalc(times[ii])
            ; in-situ position
            rx = rgsm[ii,0]
            ry = rgsm[ii,1]
            rz = rgsm[ii,2]
            ; in-situ B field.
            geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            geopack_t89, par, rx,ry,rz, dbx,dby,dbz
            b0gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor
        bmod_var = prefix+'bmod_gsm'
        store_data, bmod_var, times, b0gsm
        add_setting, bmod_var, /smart, {$
            display_type: 'vector', $
            unit: 'nT', $
            short_name: 'B!S!UT89!N!S', $
            coord: 'GSM', $
            coord_labels: ['x','y','z']}


        ; Interpolate to common times.
        common_times = make_bins(time_range, time_step)
        vars = tnames(prefix+'*')
        foreach var, vars do interp_time, var, common_times

        ; Calculate tilt angle
        vars = prefix+['b','bmod']
        foreach var, vars do begin
            get_data, var+'_gsm', times, bgsm
            bsm = cotran(bgsm, times, 'gsm2sm')
            tilt = atan(bsm[*,2],snorm(bsm[*,0:1]))*deg
            the_var = var+'_tilt'
            store_data, the_var, times, tilt
            add_setting, the_var, /smart, {$
                display_type: 'scalar', $
                unit: 'deg', $
                short_name: 'Tilt'}
        endforeach

        ; Calculate the tilt angle with model subtracted.
        tilt_var = prefix+'db_tilt'
        sys_subtract, prefix+'b_tilt', prefix+'bmod_tilt', to=tilt_var
        add_setting, tilt_var, /smart, {$
            display_type: 'scalar', $
            unit: 'deg', $
            short_name: 'Tilt'}
    endforeach


    omni_read_index, time_range


;---Generate the plot.
    secofday = 86400d
    duration = total(time_range*[-1,1])
    nday = duration/secofday
    ysize = 6.
    scale = 1.8
    xsize = ysize*scale*nday ; normalize by duration.
    tilt_range = [-1,1]*64
    psym = 8
    symsize = 0.5
    label_size = 0.7
    usersym, [1,1,-1,-1,1]*0.1, [-1,1,1,-1,-1]*1
    ct = 70
    smooth_width = 3600./5

    file = join_path([homedir(),'fig_azim_prop_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_storm_overview2.pdf'])
    if keyword_set(test) then file = 0
    sgopen, file, xsize=xsize, ysize=ysize


    npanel = 4
    poss = sgcalcpos(npanel, xchsz=xchsz, ychsz=ychsz, ypans=[1,1,2,2], tmargin=2, rmargin=8, lmargin=12)
    labels = letters(npanel)+'. '+['Dst','AE','UT-MLT','UT-X!USM!N']
    for ii=0, npanel-1 do xyouts, xchsz*2, poss[3,ii]-ychsz*0.8, /normal, labels[ii]


;---The ewogram.
    tpos = poss[*,2]
    xrange = [0,duration]
    yrange = mlt_range
    zrange = tilt_range
    xticklen = -0.02
    yticklen = xticklen*(tpos[3]-tpos[1])
    yminor = 3
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    xminor = 8  ; hour
    xtickv = make_bins(time_range[0]+xrange, 3600*xminor, /inner)-time_range[0] ; make time line up at hours.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)+' '
    ytitle = 'MLT (hr)'

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
        ystyle=5, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, $
        noerase=1, nodata=1

    foreach mission_probe, mission_probes do begin
        prefix = mission_info[mission_probe].prefix
        if tnames(prefix+'db_tilt') eq '' then continue
        get_data, prefix+'db_tilt', times, tilt
        get_data, prefix+'mlt', times, mlt
        ntime = n_elements(times)
        xxs = times-time_range[0]
        yys = mlt
        index = where_pro(yys, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        get_data, prefix+'r_sm', time, rsm
        index = where_pro(rsm[*,0], '][', xsm_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        index = where(finite(tilt,/nan), count)
        if count ne 0 then xxs[index] = !values.f_nan


        zzs = tilt
        zzs -= smooth(zzs,80,/edge_truncate, /nan)
        ; scale with mlt.
        mean = mean(zzs,/nan)
        zzs = (zzs-mean)*sqrt(abs(mlt)+1)+mean
        ;zzs -= smooth(zzs, smooth_width, /edge_truncate, /nan)
        ;zzs = deriv(tilt)
        zzs = data2color(zzs, range=zrange, ct=ct, /reverse_ct, index_color=index_color)

        if count ne 0 then yys[index] = !values.f_nan
        for ii=0, ntime-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    plot, xrange, yrange, position=tpos, $
        xstyle=1, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
        ystyle=1, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, $
        noerase=1, nodata=1
    plots, xrange, [0,0], linestyle=1, color=sgcolor('black')

;---Color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    ztitle = 'Detrended tilt angle (deg)!Cweighted by (1+|MLT|)!U0.5!N'
    ztickv = [-80,-20,-5,0,5,20,80]
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks, zcharsize=label_size



;---The keogram.
    tpos = poss[*,3]
    xrange = [0,duration]
    yrange = xsm_range
    zrange = tilt_range
    xticklen = -0.02
    yticklen = xticklen*(tpos[3]-tpos[1])
    yminor = 4
    ytickv = make_bins(yrange,yminor, /inner)
    yticks = n_elements(ytickv)-1
    xminor = 8  ; hour
    xtickv = make_bins(time_range[0]+xrange, 3600*xminor, /inner)-time_range[0] ; make time line up at hours.
    xticks = n_elements(xtickv)-1
    xtickn = strarr(xticks+1)
    for ii=0, xticks do begin
        the_time = time_range[0]+xtickv[ii]
        xtickn[ii] = time_string(the_time,tformat='hh:mm')
        date = time_string(the_time,tformat='YYYY-MM-DD')
        if ii eq 0 then begin
            xtickn[ii] += '!C'+date
            continue
        endif
        if the_time mod secofday ne 0 then continue
        xtickn[ii] += '!C'+date
    endfor
    ytitle = 'X!DSM!N (Re)'

    plot, xrange, yrange, position=tpos, $
        xstyle=5, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
        ystyle=5, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, $
        noerase=1, nodata=1

    foreach mission_probe, mission_probes do begin
        prefix = mission_info[mission_probe].prefix
        if tnames(prefix+'db_tilt') eq '' then continue
        get_data, prefix+'db_tilt', times, tilt
        get_data, prefix+'r_sm', times, rsm
        ntime = n_elements(times)
        xxs = times-time_range[0]
        yys = rsm[*,0]

        index = where_pro(yys, '][', yrange, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        get_data, prefix+'mlt', time, mlt
        index = where_pro(mlt, '][', mlt_range, count=count)
        if count ne 0 then xxs[index] = !values.f_nan
        index = where(finite(tilt,/nan), count)
        if count ne 0 then xxs[index] = !values.f_nan


        zzs = tilt
        zzs -= smooth(zzs,80,/edge_truncate, /nan)
        ; scale with mlt.
        mean = mean(zzs,/nan)
        zzs = (zzs-mean)*sqrt(abs(mlt)+1)+mean
        ;zzs -= smooth(zzs, smooth_width, /edge_truncate, /nan)
        ;zzs = deriv(tilt)
        zzs = data2color(zzs, range=zrange, ct=ct, /reverse_ct, index_color=index_color)

        if count ne 0 then yys[index] = !values.f_nan
        for ii=0, ntime-1 do plots, xxs[ii], yys[ii], color=zzs[ii], psym=psym, symsize=symsize
    endforeach

    plot, xrange, yrange, position=tpos, $
        xstyle=1, xticklen=xticklen, xticks=xticks, xtickv=xtickv, xminor=xminor, xtickname=xtickn, $
        ystyle=1, yticklen=yticklen, yticks=yticks, ytickv=ytickv, yminor=yminor, ytitle=ytitle, $
        noerase=1, nodata=1
    plots, xrange, [0,0], linestyle=1, color=sgcolor('black')

;---Color bar.
    cbpos = tpos[[2,1,2,3]]+[1,0,2,0]*xchsz
    ztitle = 'Detrended tilt angle (deg)!Cweighted by (1+|MLT|)!U0.5!N'
    ztickv = [-80,-20,-5,0,5,20,80]
    ztickv = [-64,-16,-4,0,4,16,64]
    ztickn = string(ztickv,format='(I0)')
    ztickv = alog10(abs(ztickv)>1)*sign(ztickv)
    zrange = alog10(abs(zrange))*[-1,1]
    zticks = n_elements(ztickv)-1
    sgcolorbar, 256-findgen(256), ct=ct, zrange=zrange, position=cbpos, ztitle=ztitle, ztickv=ztickv, ztickn=ztickn, zticks=zticks, zcharsize=label_size



;---Dst and AE.
    vars = ['dst','ae']
    options, vars, 'yticklen', yticklen
    options, vars, 'ystyle', 1
    options, vars, 'yminor', 5
    ystep = [50,500]
    constant = [0,500]
    foreach var, vars, ii do begin
        get_data, var, times, data
        ytickv = make_bins(data, ystep[ii])
        yticks = n_elements(ytickv)-1
        options, var, 'ytickv', ytickv
        options, var, 'yticks', yticks
        options, var, 'yrange', minmax(ytickv)
        options, var, 'constant', constant[ii]
    endforeach

    tplot, vars, trange=time_range, position=poss[*,0:1], /nouttick, /noerase
    if keyword_set(test) then stop
    sgclose

end


;---The ten storms.
; Read from /Volumes/GoogleDrive/My Drive/works/works/azim_prop/ten_event_paper/data/storm_list.txt
;    storm_list = list()
;    storm_list.add, time_double(['2013-08-27/08:41:00','2013-08-28/16:26:00'])
;    storm_list.add, time_double(['2013-10-01/17:40:00','2013-10-03/10:41:00'])
;    storm_list.add, time_double(['2013-10-08/12:13:00','2013-10-09/16:11:00'])
;    storm_list.add, time_double(['2014-08-26/20:51:00','2014-08-29/05:29:00'])
;    storm_list.add, time_double(['2014-09-12/10:30:00','2014-09-13/12:08:00'])
;    storm_list.add, time_double(['2016-10-12/21:37:04','2016-10-14/21:58:56'])
;    storm_list.add, time_double(['2016-10-25/01:02:56','2016-10-28/11:58:56'])
;    storm_list.add, time_double(['2016-10-28/13:14:00','2016-10-29/23:58:00'])
;    storm_list.add, time_double(['2017-03-01/08:49:04','2017-03-02/14:28:00'])
;    storm_list.add, time_double(['2017-03-27/00:34:56','2017-03-28/17:36:00'])


    ;storm_list.add, time_double(['2014-08-28/08:00','2014-08-29/00:00'])
;    storm_list.add, time_double(['2007-03-23/00:00','2007-03-24/00:00'])
;    storm_list.add, time_double(['2009-02-27/00:00','2009-02-28/00:00'])
    foreach time_range, storm_list do azim_prop_gen_storm_keox_ewo, time_range

end
