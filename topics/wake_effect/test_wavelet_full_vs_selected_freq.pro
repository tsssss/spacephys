;+
; Calculate the frequency spectrum at full frequency and selected frequency.
;-


;---Load data.
    test = 0
    duration = 180.
    secofday = 86400d

    date = '1998-09-14'
    test_times = ['22:30','02:00','02:30','21:30','22:30','23:30','06:30','07:30','08:30','09:30','10:30','13:30','15:30']
    ; peaks at 2f.
    test_times = ['21:03','21:00']
    var = 'po_e_uv'
    mission = 'polar'

    date = '2012-11-14'
    test_times = ['22:23','09:23', '10:00','11:00','12:00']
    date = '2013-05-01'
    test_times = ['19:24']
    probe = 'b'
    mission = 'rbsp'
    var = 'rbsp'+probe+'_e_uv'

    case mission of
        'polar': spinrate = polar_info('spin_rate')
        'rbsp': spinrate = rbsp_info('spin_rate')
        'themis': spinrate = themis_info('spin_rate')
    endcase
    full_time_range = time_double(date)+[0,secofday]

    foreach time, test_times do begin
        xticklen = -0.02
        yticklen = -0.02
        unit = 'mV/m'

        time_range = time_double(date+'/'+time)+[0,duration]
        if tnames(var) ne '' then begin
            get_data, var, times
            index = where_pro(times,time_range, count=count)
            if count eq 0 then store_data, var, /delete
        endif
        if tnames(var) eq '' then begin
            prefix = strmid(get_prefix(var),0,2)
            if prefix eq 'po' then begin
                polar_read_euv, full_time_range
            endif else if prefix eq 'rb' then begin
                probe = strmid(var,4,1)
                rbsp_read_euv, full_time_range, probe=probe
            endif else if prefix eq 'th' then begin
                probe = strmid(var,2,1)
                themis_read_euv, full_time_range, probe=probe
            endif
        endif
        get_data, var, times, data, limits=lim
        index = where_pro(times, time_range, count=count)
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


    ;---Plot.
        base_name = get_prefix(data_var)+'eu_'+time_string(time_range[0],tformat='YYYY_MMDD_hhmm')+'.pdf'
        file = join_path([srootdir(),'plot',base_name])
        magnify = 1
        if keyword_set(test) then begin
            file = 0
            magnify = 2
        endif
        sgopen, file, xsize=8, ysize=2, magnify=magnify
        poss = sgcalcpos(1,2,xpans=[3,1.2], xpad=10, $
            tmargin=2, bmargin=4, rmargin=2, xchsz=xchsz, ychsz=ychsz)

    ;---Calculate the psd on full frequency.
        scales = wavelet_make_scale(times)
        calc_psd, data_var, scales=scales


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
        ytitle = '|'+unit+'|/Hz!U1/2!N'
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
        plot, xs, ys, /nodata, $
            xlog=xlog, xstyle=1, xticklen=xticklen, xtitle=xtitle, $
            ylog=ylog, ystyle=1, yticklen=yticklen, ytitle=ytitle, yrange=yrange, $
            yticks=yticks, ytickv=ytickv, ytickname=ytickn, yminor=yminor, $
            /noerase, position=tpos
        oplot, xs,ys, color=sgcolor('silver')

        tx = tpos[0]
        ty = tpos[3]+ychsz*0.35
        xyouts, tx,ty,/normal, 'b. Amplitude'


    ;---Calculate spectrum at given frequency.
        label_size = 0.7
        
        spin_freq = 1d/spinrate
        base = 2d^(1d/4)
        freq_range = [0.5,8]/spinrate
        selected_freqs = smkgmtrc(min(freq_range),max(freq_range),base,'dx')
        nselected_freq = n_elements(selected_freqs)
        selected_periods = 1d/selected_freqs
        wavelet_info = wavelet_info()
        t2s = wavelet_info[7]
        selected_scales = selected_periods*t2s
        calc_psd, data_var, scales=selected_scales
        get_data, psd_var, xs, ys
        ys = sqrt(ys)
        xs = 1d/xs
        plots, xs,ys,/data, psym=1, symsize=label_size*0.5, color=sgcolor('blue')
        
        foreach number, [1d,2,3,4] do begin
            period = spinrate/number
            plots, period+[0,0], yrange, linestyle=1
            msg = (number eq 1)? 'f  ': string(number,format='(I0)')+'f'
            tmp = convert_coord(period,yrange[1], /data, /to_normal)
            tx = tmp[0]
            ty = tmp[1]-ychsz*1
            xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=label_size
            if number eq 1 then xyouts, tx,ty,/normal, alignment=0.1, 'spin', charsize=label_size*0.8
        endforeach
        
        if keyword_set(test) then stop
        sgclose
    endforeach
end
