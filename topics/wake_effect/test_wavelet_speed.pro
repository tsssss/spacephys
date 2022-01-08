;+
; Calculate the frequency spectrum at full frequency and selected frequency.
;-


;---Load data.
    test = 0
    duration = 180.
    secofday = 86400d

    date = '1998-09-14'
    mission = 'polar'
    var = 'po_e_uv'

    date = '2013-11-14'
    mission = 'rbsp'
    probe = 'a'
    var = 'rbsp'+probe+'_e_uv'

    full_time_range = time_double(date)+[0,secofday]
    ntest_time = round((full_time_range[1]-full_time_range[0])/duration)
    test_times = smkarthm(full_time_range[0], duration, ntest_time, 'x0')

    case mission of
        'polar': spinrate = polar_info('spin_rate')
        'rbsp': spinrate = rbsp_info('spin_rate')
        'themis': spinrate = themis_info('spin_rate')
    endcase

    data_var = var+'_section'
    psd_var = data_var+'_psd_cwt'
    base = 2d^0.25
    selected_freqs = base^(findgen(11)-1)/spinrate
    selected_periods = 1d/selected_freqs
    wavelet_info = wavelet_info()
    t2s = wavelet_info[6]
    selected_scales = selected_periods*t2s
    if tnames(var) ne '' then begin
        get_data, var, times
        index = lazy_where(times,time_range, count=count)
        if count eq 0 then store_data, var, /delete
    endif
    if tnames(var) eq '' then begin
        if mission eq 'polar' then begin
            polar_read_euv, full_time_range
        endif else if prefix eq 'rbsp' then begin
            probe = strmid(var,4,1)
            rbsp_read_euv, full_time_range, probe=probe
        endif else if prefix eq 'themis' then begin
            probe = strmid(var,2,1)
            themis_read_euv, full_time_range, probe=probe
        endif
    endif


;---Plot.
    base_name = get_prefix(data_var)+'eu_'+time_string(full_time_range[0],tformat='YYYY_MMDD')+'.pdf'
    file = join_path([srootdir(),'plot',base_name])
    magnify = 1
    if keyword_set(test) then begin
        file = 0
        magnify = 2
    endif
    sgopen, file, xsize=4, ysize=3, magnify=magnify
    poss = sgcalcpos(1, tmargin=2, bmargin=4, rmargin=2, lmargin=8, xchsz=xchsz, ychsz=ychsz)
    tpos = poss[*,0]

    xticklen = -0.02
    yticklen = -0.02
    unit = 'mV/m'
    xtitle = 'Period (sec)'
    xticklen = xticklen
    xlog = 1
    ytitle = '|'+unit+'|/Hz!U1/2!N'
    yticklen = yticklen
    ylog = 1
    yrange = [-3,3]
    ytickv = make_bins(yrange,1)
    yrange = 10.^minmax(ytickv)
    yticks = n_elements(ytickv)-1
    ytickn = strarr(yticks+1)
    for ii=0, yticks do begin
        case ytickv[ii] of
            0: ytickn[ii] = '1'
            else: ytickn[ii] = '10!U'+string(ytickv[ii],format='(I0)')
        endcase
    endfor
    if yticks ge 6 then begin
        target = 3
        for ii=0, yticks do begin
            if ii mod target eq 0 then continue
            ytickn[ii] = ' '
        endfor
    endif
    ytickv = 10.^ytickv
    yminor = 10

    xrange = [0.5,20]
    xtickv = [1.,10]

    label_size = 0.7
    psym = 3    ; dot.
    symsize = 0.1

    plot, xrange, yrange, /nodata, $
        xlog=xlog, xstyle=1, xticklen=xticklen, xtitle=xtitle, xrange=xrange, xtickv=xtickv, $
        ylog=ylog, ystyle=1, yticklen=yticklen, ytitle=ytitle, yrange=yrange, $
        yticks=yticks, ytickv=ytickv, yminor=yminor, ytickname=ytickn, $
        /noerase, position=tpos

    tx = tpos[0]
    ty = tpos[3]+ychsz*0.35
    xyouts, tx,ty,/normal, 'a. Amplitude '+strjoin(strupcase((strsplit(data_var,'_',/extract))[0]),' ')

    foreach number, [1d,2,3] do begin
        period = spinrate/number
        plots, period+[0,0], yrange, linestyle=1
        msg = (number eq 1)? 'f  ': string(number,format='(I0)')+'f'
        tmp = convert_coord(period,yrange[1], /data, /to_normal)
        tx = tmp[0]
        ty = tmp[1]-ychsz*1
        xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size
        if number eq 1 then xyouts, tx,ty,/normal, alignment=0.1, 'spin', charsize=label_size*0.8
    endforeach



;---Method 1, for sections.
    time = test_times[0]
    time_range = time+[0,duration]
    get_data, var, times, data, limits=lim
    index = lazy_where(times, time_range, count=count)
    if count eq 0 then message, 'No data in the given time ...'
    stplot_index, var, 0, newname=data_var
    times = times[index]
    data = data[index,0]
    store_data, data_var, times, data

;---Calculate the psd on selected frequencies.
    tic
    calc_psd, data_var, scales=selected_scales
    result1 = toc()

;---Method 2, as a whole.
    tic
    u_var = var+'_u'
    get_data, var, times, data
    store_data, u_var, times, data[*,0]
    calc_psd, u_var, scales=selected_scales
    result2 = toc()

    print, 'In sections (sec): '+string(round(result1*ntest_time))
    print, 'As a whole  (sec): '+string(round(result2))
    if keyword_set(test) then stop
    sgclose
end
