;+
; Plot the psd using wavelet at [1,3,9] of the f_spin.
;-
pro plot_psd_at_freq, time, duration=duration, var=var, freqs=freqs, plot_freqs=plot_freqs

    if n_elements(duration) eq 0 then duration = 180.   ; sec.
    if n_elements(time) eq 0 then message, 'No input time ...'
    if tnames(var) eq '' then message, 'No input variable ...'

    xticklen = -0.02
    yticklen = -0.02
    unit = 'mV/m'

    test = 0

    time_range = time[0]+[0,duration]
    get_data, var, times, data, limits=lim
    index = lazy_where(times, time_range, count=count)
    if count eq 0 then message, 'No data in the given time ...'
    data_var = var+'_section'
    stplot_index, var, 0, newname=data_var
    times = times[index]
    data = data[index,0]
    yrange = sg_autolim(data)
    yrange = [-1,1]*max(abs(yrange))
    store_data, data_var, times, data
    options, data_var, 'labels', ''
    options, data_var, 'xticklen', xticklen*2
    options, data_var, 'yticklen', yticklen*0.5
    options, data_var, 'yrange', yrange
    options, data_var, 'ystyle', 1
    options, data_var, 'yticks', 4
    options, data_var, 'yminor', 5
    options, data_var, 'constant', 0
    options, data_var, 'ytitle', '('+unit+')'

    if n_elements(freqs) le 1 then begin
        scales = wavelet_make_scale(times)
    endif else begin
        scales = 1d/freqs
        scales = scales[sort(scales)]
    endelse

    calc_psd, data_var, scales=scales


;---Plot.
    base_name = get_prefix(data_var)+'eu_'+time_string(time,tformat='YYYY_MMDD_hhmm')+'.pdf'
    file = join_path([srootdir(),'plot',base_name])
    magnify = 1
    if keyword_set(test) then begin
        file = 0
        magnify = 2
    endif
    sgopen, file, xsize=8, ysize=2, magnify=magnify
    poss = sgcalcpos(1,2,xpans=[3,1.2], xpad=10, $
        tmargin=2, bmargin=4, rmargin=2, xchsz=xchsz, ychsz=ychsz)

    tpos = poss[*,0]
    tplot, data_var, position=tpos, trange=time_range
    tx = tpos[0]
    ty = tpos[3]+ychsz*0.35
    xyouts, tx,ty,/normal, 'a. Waveform '+strjoin(strupcase((strsplit(data_var,'_',/extract))[0]),' ')

    tpos = poss[*,1]
    psd_var = data_var+'_psd_cwt'
    get_data, psd_var, xs, ys
    ys = sqrt(ys)
    xs = 1d/xs
    xtitle = 'Period (sec)'
    xticklen = xticklen
    xlog = 1
    ytitle = '|'+unit+'|'
    yticklen = yticklen
    ylog = 1

    if ylog eq 1 then begin
        yrange = minmax(alog10(ys))
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
    endif
    plot, xs, ys, $
        xlog=xlog, xstyle=1, xticklen=xticklen, xtitle=xtitle, $
        ylog=ylog, ystyle=1, yticklen=yticklen, ytitle=ytitle, yrange=yrange, $
        yticks=yticks, ytickv=ytickv, ytickname=ytickn, yminor=yminor, $
        /noerase, position=tpos

    label_size = 0.7
    if n_elements(plot_freqs) ne 0 then begin
        periods = 1d/plot_freqs
        numbers = round(max(periods)/periods)
        foreach period, periods, ii do begin
            plots, period+[0,0], yrange, linestyle=1
            tmp = convert_coord(period, yrange[1], /data, /to_normal)
            tx = tmp[0]
            ty0 = tmp[1]
            ty = ty0-ychsz

            if period eq max(periods) then begin
                msg = 'T!Dspin'
                alignment = 0.25
                tchsize = label_size
            endif else begin
                msg = 'T/'+string(numbers[ii],format='(I0)')
                alignment = 0.5
                tchsize = label_size
            endelse
            if numbers[ii] eq 2 then ty = ty-ychsz*label_size
            xyouts, tx,ty,/normal,alignment=alignment, charsize=tchsize, msg
        endforeach
    endif

    tx = tpos[0]
    ty = tpos[3]+ychsz*0.35
    xyouts, tx,ty,/normal, 'b. Amplitude'

    if keyword_set(test) then stop
    sgclose

end



    duration = 180.
    secofday = 86400d


;---RBSP.
;    date = '2012-11-14'
;    times = ['22:23','09:23', '10:00','11:00','12:00']
;    probe = 'a'
;    time_range = time_double(date)+[0,secofday]
;    ;rbsp_read_euv, time_range, probe=probe
;    spinrate = rbsp_info('spin_rate')
;    plot_freqs = [1.,2,3,9]/spinrate
;    var = 'rbsp'+probe+'_e_uv'
;    foreach time, times do $
;        plot_psd_at_freq, time_double(date+'/'+time), $
;        duration=duration, var=var, plot_freqs=plot_freqs


;---Polar.
    ; 1998-09-14/03:00, washing machine mode?
    ; 1998-09-14/21:00, washing machine mode?
    ; weak wake, 01:30,18:30.
    ; no wake. 02:30,21:30,22:30,23:30
    ; wave. 06:30,07:30,08:30,09:30,10:30,19:30
    ; marginally wake. 13:30,14:30,15:30

    date = '1998-09-14'
    times = ['22:30','02:00','02:30','21:30','22:30','23:30','06:30','07:30','08:30','09:30','10:30','13:30','15:30']
    time_range = time_double(date)+[0,secofday]
    ;polar_read_euv, time_range
    spinrate = polar_info('spin_rate')
    plot_freqs = [1.,2,3,9]/spinrate
    var = 'po_e_uv'
    foreach time, times do $
        plot_psd_at_freq, time_double(date+'/'+time), $
        duration=duration, var=var, plot_freqs=plot_freqs
end
