;+
; Plot Poynting flux stuff.
;-

;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 0

;---Settings.
    einfo = get_var_data(info_var)
    time_range = einfo.efield['time_range']
    time_range = time_double(['2014-08-28/10:00','2014-08-28/10:30'])
    scale_info = {s0:10d, s1:1000, dj:1d/8, ns:0d}

    fac = ['||','east','north']
    rgb = sgcolor(['red','green','blue'])
    probes = ['a','d','e']
    pf_unit = 'mW/m!U2!N'


;---Calculate pflux.
    ; Need position, E and B field in a same coordinate.
    ; Get S in FAC, E/B ratio, PSD of the E and B.
    foreach probe, probes do begin
        pre0 = 'th'+probe+'_'

        ; Uniform time.
        vars = pre0+['b0_gsm','db_gsm','r_gsm']
        foreach var, vars do interp_time, var, to=pre0+'e_gsm'

        ; Convert E/B fields to FAC.
        define_fac, pre0+'b0_gsm', pre0+'r_gsm'
        vars = pre0+['db_gsm','e_gsm']
        foreach var, vars do to_fac, var

        vars = pre0+['db','e']+'_fac'
        short_names = ['dB','dE']
        foreach var, vars, ii do begin
            add_setting, var, /smart, {$
                display_type: 'vector', $
                short_name: short_names[ii], $
                coord: 'FAC', $
                coord_labels: fac, $
                colors: rgb}
        endforeach

        ; Calculate pflux.
        stplot_calc_pflux_mor, pre0+'e_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scale_info, /power
        rename_var, pre0+'pf_fac', to=pre0+'pf_power_fac'
        stplot_calc_pflux_mor, pre0+'e_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scale_info;, /power

        scales = wavelet_make_scale(s0=scale_info.s0, sj=scale_info.s1, dj=scale_info.dj)
        vars = pre0+['db','e']
        foreach var, vars do begin
            get_data, var+'_gsm', times, data
            index = lazy_where(times, time_range)
            store_data, var+'_mag', times[index], snorm(data[index,*])
            calc_psd, var+'_mag', scales=scales
        endforeach

        tvar = pre0+'pf_fac_mor_spec_1'
        get_data, tvar, uts, dat
        index = lazy_where(uts, time_range, count=nrec)
        spsd = total(dat[index,*], 1)/nrec


        get_data, pre0+'db_mag_psd_cwt', freqs, bpsd
        get_data, pre0+'e_mag_psd_cwt', freqs, epsd
        eb_ratio = sqrt(epsd/bpsd)*1e3  ; in km/s.
        add_setting, pre0+'pf_fac', /smart, {$
            display_type: 'vector', $
            unit: 'mW/m!U2!N', $
            short_name: 'S', $
            coord: 'FAC', $
            coord_labels: fac, $
            colors: rgb, $
            freq: freqs, $
            freq_unit: 'Hz', $
            constant: 0., $
            e_psd: epsd, $
            e_psd_unit: '(mV/m)!U2!N/Hz', $
            b_psd: bpsd, $
            b_psd_unit: '(nT)!U2!N/Hz', $
            s_psd: spsd, $
            eb_ratio: eb_ratio, $
            eb_ratio_unit: 'km/s'}
    endforeach
    



    foreach probe, probes do begin
        pre0 = 'th'+probe+'_'

        get_data, pre0+'pf_fac', times, pffac
        cmap = get_var_data(pre0+'cmap', at=times)
        pfpara = pffac[*,0]*cmap
        get_data, pre0+'pf_power_fac', times, pfpfac
        pfppara = pfpfac[*,0]*cmap
        store_data, pre0+'pf_para', times, [[pfpara],[pfppara]], limits={$
            ytitle: 'S!D||!N @100km!C(mW/m!U2!N)', constant: 0, colors:[0,6], yticks:2, yminor:4}

        if keyword_set(test) then begin
            file = test
            magnify = 2
        endif else begin
            file = sparentdir(srootdir())+'/plot/fig_case2_pflux_detail_th'+probe+'.pdf'
            magnify = 1
        endelse
        sgopen, file, xsize=4, ysize=3, magnify=magnify

        poss = sgcalcpos(2,2, xpans=[3,1], ypans=[1,1], xpad=1, xchsz=xchsz,ychsz=ychsz, rmargin=2)
        ticklen = -0.02
        str_mu ='!9'+string(109b)+'!X'
        unit = '('+str_mu+'W/m!U2!N)'

        tvar = pre0+'pf_fac'
        freqs = get_setting(tvar, 'freq')*1e3
        freq_range = minmax(freqs)
        ylogrange = alog10(minmax(freq_range)) & ylogrange = round(ylogrange)
        freq_range = 10.^ylogrange
        eb_ratio = get_setting(tvar, 'eb_ratio')

        cb_pos = poss[[0,3,2,3],0,0]+[0,0.5,0,1]*ychsz
        tvar = pre0+'pf_para_mor_spec'
        get_data, pre0+'pf_fac_mor_spec_1', times, data, limits=lim & data *= 1e3
        store_data, tvar, times, data, freqs, limits=lim
        options, tvar, 'zposition', cb_pos
        options, tvar, 'zhorizontal', 1
        options, tvar, 'ystyle', 1
        options, tvar, 'ylog', 1
        options, tvar, 'yrange', freq_range
        options, tvar, 'ytitle', 'Freq (mHz)'
        options, tvar, 'ztitle', 'S '+unit
        options, tvar, 'zrange', [-1,1]*max(abs(data))*0.2

        spec_colortable = 66
        device, decomposed=0
        loadct2, spec_colortable
        tplot, pre0+['pf_para_mor_spec','pf_para'], position=reform(poss[*,0,*]), trange=time_range, /novtitle
        ; Add label.
        get_data, pre0+'pf_para', limits=lim
        colors=lim.colors
        tpos = poss[*,0,1]
        tx = tpos[0]+xchsz*0.5
        if probe eq 'd' then tx+= 0.6*(tpos[2]-tpos[0])
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,/normal, 'Instantaneous', charsize=0.7, color=colors[0]
        ty = tpos[3]-ychsz*2*0.8
        xyouts, tx,ty,/normal, 'Period-averaged', charsize=0.7, color=colors[1]
        device, decomposed=1
        ; Add DF time.
        case probe of
            'a': df_time = time_double(['2014-08-28/10:20:39', '2014-08-28/10:24:21'])
            'd': df_time = time_double(['2014-08-28/10:12:15', '2014-08-28/10:13:54'])
            'e': df_time = time_double(['2014-08-28/10:20:00', '2014-08-28/10:22:30'])
        endcase
        txs = dblarr(2)
        plot, time_range, [0,1], /nodata, /noerase, xstyle=5,ystyle=5, position=tpos
        for ii=0,1 do begin
            tmp = convert_coord([df_time[ii],0],/data,/to_normal)
            txs[ii] = tmp[0]
        endfor
        ty = tpos[3]
        thick = (size(file,/type) eq 7)? 8:4
        plots, txs,ty+[0,0],/normal, thick=thick, color=sgcolor('red')

        tpos = poss[*,1,0]
        txs = eb_ratio
        xlogrange = alog10(minmax(txs)) & xlogrange = [floor(xlogrange[0]),ceil(xlogrange[1])]
        xrange = 10.^xlogrange
        xminor = 10
        xticks = (xlogrange[1]-xlogrange[0])
        xlogtickv = smkarthm(xlogrange[0],xlogrange[1],xticks+1,'n')
        xtickv = 10.^xlogtickv
        xtickn = '10!U'+string(xlogtickv,format='(I0)')
        xticklen = ticklen
        yrange = freq_range
        plot, txs, freqs, position=tpos, $
            xstyle=9, xlog=1, xtitle='', xrange=xrange, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, xticklen=xticklen, $
            ystyle=1, ylog=1, yrange=yrange, ytickformat='(A1)', $
            /noerase, /nodata
        color1 = sgcolor('gold')
        oplot, txs, freqs, color=color1
        tx = mean(tpos[[0,2]])
        ty = tpos[1]-ychsz*3
        xyouts, tx,ty,/normal, 'R!DE/B!N (km/s)', alignment=0.5, color=color1

        ; 2nd axis for pflux.
        color2 = sgcolor('black')
        tvar = pre0+'pf_fac'
        txs = get_setting(tvar, 's_psd')*1e3
        xtitle = 'S '+unit
        xrange = sg_autolim([-1,1]*abs(max(txs,/absolute)))
        xminor = 5
        xticks = 2
        xtickv = smkarthm(xrange[0],xrange[1],xticks+1,'n')
        xticklen = ticklen
        axis, xaxis=1, /save, xstyle=1, xlog=0, xtitle=xtitle, xrange=xrange, xminor=xminor, xticks=xticks, xtickv=xtickv, xticklen=xticklen, color=color2
        plots, txs, freqs, color=color2
        plots, [0,0], yrange, color=color2, linestyle=1

        if keyword_set(test) then stop
        sgclose
    endforeach


end
