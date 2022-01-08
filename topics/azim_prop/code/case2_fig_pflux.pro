;+
; Plot E, B, S in FAC, and the power spectrum of E, B, S.
;
; de_var, the FAC dE in mV/m.
; db_var, the FAC dB in nT.
; pf_var, the FAC S in mW/m^2.
;-


;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 1


;---Settings.
    time_range = time_double(['2014-08-28/10:00','2014-08-28/10:30'])

    test_times = time_double(['2014-08-28/10:10','2014-08-28/10:20'])
    plot_time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])


    fac = ['||','east','north']
    rgb = sgcolor(['red','green','blue'])
    probes = ['d','e','a']
    pf_unit = 'mW/m!U2!N'
    ct = 66
    va_color = sgcolor('silver')
    nprobe = n_elements(probes)
    fig_labels = letters(nprobe)

    xticklen_chsz = -0.15
    yticklen_chsz = -0.40



;---Calculate Alfven speed.
    v_alfven = list()
    va0 = 22.0  ; km/s, for B in nT, n in cc, m in atomic mass
    foreach probe, probes do begin
        prefix = 'th'+probe+'_'
        bmag = snorm(get_var_data(prefix+'b_gsm',in=time_range))
        if check_if_update(prefix+'ele_n',time_range) then themis_read_density, time_range, probe=probe
        den = get_var_data(prefix+'ele_n',in=time_range)
        va = va0*mean(bmag)/sqrt(mean(den)*1)   ; 1 for H+.
        v_alfven.add, va
    endforeach
    v_alfven = v_alfven.toarray()




;---Plot.
    root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
    foreach probe, probes, probe_id do begin
        pre0 = 'th'+probe+'_'

        get_data, pre0+'pf_fac', times, pffac
        cmap = get_var_data(pre0+'cmap', at=times)
        pfpara = pffac[*,0]*cmap
        store_data, pre0+'pf_para', times, pfpara, limits={$
            ytitle: 'S!D||!N @100km!C(mW/m!U2!N)', constant: 0, yticks:2, yminor:4}

        if keyword_set(test) then begin
            file = test
            magnify = 2
        endif else begin
            file = join_path([root_dir,'case2_fig_pflux_th'+probe+'.pdf'])
            magnify = 1
        endelse
        sgopen, file, xsize=5, ysize=3, magnify=magnify
        margins = [8,4,2,4]
        poss = sgcalcpos(2,2, xpans=[3,1], ypans=[1,1], xpad=2, xchsz=xchsz,ychsz=ychsz, margins=margins)
        ticklen = -0.02
        str_mu ='!9'+string(109b)+'!X'
        unit = '('+str_mu+'W/m!U2!N)'

        tvar = pre0+'pf_fac'
        freqs = get_setting(tvar, 'freq')*1e3
        freq_range = minmax(freqs)
        ;ylogrange = alog10(minmax(freq_range)) & ylogrange = round(ylogrange)
        ;freq_range = 10.^ylogrange
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
        options, tvar, 'color_table', ct

        tvar = pre0+'pf_para'
        get_data, tvar, times, data
        yrange = minmax(data)*0.1
        yrange = [floor(yrange[0]),ceil(yrange[1])]*10
        options, tvar, 'yrange', yrange
        options, tvar, 'yminor', 5
        options, tvar, 'ytitle', '(mW/m!U2!N)'


        tpos = reform(poss[*,0,0])
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        vars = pre0+['pf_para_mor_spec','pf_para']
        options, vars, 'xticklen', xticklen
        options, vars, 'yticklen', yticklen
        options, vars, 'zhorizontal', 1
        tplot, vars, position=reform(poss[*,0,*]), trange=plot_time_range, /novtitle
        timebar, test_times, linestyle=1
        foreach var, vars, var_id do begin
            tpos = poss[*,0,var_id]
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*0.9
            msg = fig_labels[probe_id]+'-'+string(var_id+1,format='(I0)')+'.'
            xyouts, tx,ty,/normal, msg
            if var_id eq 1 then begin
                msg = 'S!D||!N @100 km'
                tx = tpos[2]-xchsz*0.5
                ty = tpos[3]-ychsz*0.9
                xyouts, tx,ty,/normal, msg, alignment=1
            endif
        endforeach


    ;---E/B ratio.
        tpos = poss[*,1,0]
        xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
        yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
        txs = eb_ratio
        xlogrange = alog10(minmax(txs)) & xlogrange = [floor(xlogrange[0]),ceil(xlogrange[1])]
        xrange = 10.^xlogrange
        xminor = 10
        xticks = (xlogrange[1]-xlogrange[0])
        xlogtickv = smkarthm(xlogrange[0],xlogrange[1],xticks+1,'n')
        xtickv = 10.^xlogtickv
        xtickn = '10!U'+string(xlogtickv,format='(I0)')
        xtitle = 'R!DE/B!N (km/s)'
        yrange = freq_range

        ; Set up coord.
        plot, txs, freqs, position=tpos, $
            xstyle=9, xlog=1, xtitle=xtitle, xrange=xrange, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, xticklen=xticklen, $
            ystyle=9, ylog=1, yrange=yrange, ytickformat='(A1)', yticklen=yticklen, $
            /noerase, /nodata

        ; Add v_A.
        polyfill, va/[1,3,3,1,1],yrange[[0,0,1,1,0]],/data, color=va_color
        tmp = (convert_coord(va,yrange[1], /data,/to_normal)+$
            convert_coord(va/3,yrange[1], /data,/to_normal))*0.5
        tx = tmp[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, 'v!DA', alignment=0.5;, color=va_color

        oplot, txs, freqs

        ; Draw box.
        plot, txs, freqs, position=tpos, $
            xstyle=1, xlog=1, xtitle=xtitle, xrange=xrange, xminor=xminor, xticks=xticks, xtickv=xtickv, xtickn=xtickn, xticklen=xticklen, $
            ystyle=1, ylog=1, yrange=yrange, ytickformat='(A1)', yticklen=yticklen, $
            /noerase, /nodata

        ; Add label.
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*0.9
        msg = fig_labels[probe_id]+'-3.'
        xyouts, tx,ty,/normal, msg


        ; Add s/c.
        tpos = poss[*,1,1]
        msg = strupcase('th'+probe)
        tx = (tpos[0]+tpos[2])*0.5
        ty = tpos[1]+ychsz*0.5
        xyouts, tx,ty,/normal, alignment=0.5, msg, charsize=1.5

        if keyword_set(test) then stop
        sgclose
    endforeach


end
