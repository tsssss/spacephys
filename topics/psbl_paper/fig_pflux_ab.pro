
; plot vsc and bmag.
test = 0

;---Plot settings.
    sgopen, 0, xsize=1, ysize=1, /inch, xchsz=abs_xchsz, ychsz=abs_ychsz
    sgclose, /wdelete

    nypanel = 5
    panel_ysize = 0.8   ; inch.
    panel_aspect_ratio = 3
    panel_xsize = panel_ysize*panel_aspect_ratio
    panel_ypad = 0.4
    upper_left_ysize = panel_ysize*nypanel+panel_ypad*(nypanel-1)*abs_ychsz
    upper_left_xsize = panel_xsize
    lower_left_ysize = panel_ysize

    xpad = 20
    ypad = 4
    margins = [12,4,8,2]
    fig_xsize = upper_left_xsize*2+(xpad+total(margins[[0,2]]))*abs_xchsz
    fig_ysize = upper_left_ysize+lower_left_ysize+(ypad+total(margins[[1,3]]))*abs_ychsz


    plot_file = join_path([srootdir(),'fig_pflux_ab.pdf'])
    if keyword_set(test) then plot_file = test
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch

    pos = margins*[abs_xchsz,abs_ychsz,abs_xchsz,abs_ychsz]/[fig_xsize,fig_ysize,fig_xsize,fig_ysize]
    pos = [pos[0:1],1-pos[2:3]]
    poss = sgcalcpos(2, ypans=[upper_left_ysize,lower_left_ysize], ypad=ypad, position=pos)
    upper_pos = poss[*,0]
    lower_pos = poss[*,1]
    poss = sgcalcpos(1,2, xpans=[1,1], xpad=xpad, position=upper_pos)
    upper_left_pos = poss[*,0]
    upper_right_pos = poss[*,1]
    upper_left_poss = sgcalcpos(nypanel, ypad=panel_ypad, position=upper_left_pos)
    upper_right_poss = sgcalcpos(nypanel, ypad=panel_ypad, position=upper_right_pos)

    nxpanel = 6
    panel_xpads = fltarr(nxpanel-1)+1 & panel_xpads[2] = 10
    lower_pos[[0,2]] += [-1,1]*abs_xchsz*2/fig_xsize
    lower_poss = sgcalcpos(1,nxpanel, xpad=panel_xpads, position=lower_pos)

;    for ii=0, nypanel-1 do begin
;        tpos = upper_left_poss[*,ii]
;        plot, [0,1],[0,1], /nodata, position=tpos, /noerase
;        tpos = upper_right_poss[*,ii]
;        plot, [0,1],[0,1], /nodata, position=tpos, /noerase
;    endfor
;    for ii=0, nxpanel-1 do begin
;        tpos = lower_poss[*,ii]
;        plot, [0,1],[0,1], /nodata, position=tpos, /noerase
;    endfor


    id = '2013_0607_event_info'
    if tnames(id) eq '' then _2013_0607_load_data3


;---Settings.
    plot_time_range = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    probes = ['a','b']
    spinrate = 11
    perp = '!9'+string(94b)+'!X'
    labfac = ['||',perp+',West',perp+',North']
    rgb = sgcolor(['red','green','blue'])
    xyz = ['x','y','z']
    dr0 = 1d/16
    default_model = 't89'

    ticklen = -0.02
    labelsize = 0.8

    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 4
    tplot_options, 'constant', 0
    tplot_options, 'labflag', -1
    tplot_options, 'ticklen', ticklen

;---Constants.
    deg = constant('deg')
    rad = constant('rad')
    re = constant('re')
    re1 = 1d/re
    r0 = 100d/re+1

    strmu = tex2str('mu')
    strsigma = tex2str('Sigma')
    fill_color = sgcolor('silver')
    label_size = constant('label_size')


;---Calculate pflux.
    filter = [0.25,100]    ; sec.
    scaleinfo = {s0:0.25d, s1:1200d, dj:1d/8, ns:0d}    ; pflux calc.
    foreach pre0, 'rb'+probes+'_' do begin
    ;---Calculate the Morlet wavelet.
        stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', $
            scaleinfo=scaleinfo

    ;---Get the pflux within the filter.
        get_data, pre0+'pf_fac', times
        nrec = n_elements(times)
        ndim = 3
        pffac = fltarr(nrec,ndim)
        for ii=0, ndim-1 do begin
            var = pre0+'pf_fac_mor_spec_'+string(ii+1,format='(I0)')
            get_data, var, times, dat, val
            index = where(val ge filter[0] and val le filter[1])
            pffac[*,ii] = total(dat[*,index], 2)
        endfor
        store_data, pre0+'pf_fac_tmp', times, pffac, limits={colors:rgb, labels:['','','']}

    ;---pflux frequency spectra, normalized version.
        tvar = pre0+'pf_fac'
        get_data, tvar, uts
        nrec = n_elements(uts)
        dr0 = sdatarate(uts)

        ; change unit to (uW/m^2), convert y from period to freq
        tvar = pre0+'pf_fac_mor_spec_1'
        get_data, tvar, uts, dat, val, limits = lim

        dat *= 1e3
        zr = [-1,1]*5
        pfunit = strmu+'W/m!U2!N'
        tvar = pre0+'pf_para_mor_spec_tmp'
        store_data, tvar, uts, dat, 1d/val, limits = lim
        options, tvar, 'zrange', zr
        options, tvar, 'ztitle', 'S!D'+labfac[0]+'!N ('+pfunit+')'
        options, tvar, 'yrange', [1e-3,4]
        options, tvar, 'ytickv', [1e-2,1e-1,1e0]
        options, tvar, 'yticks', 2
        options, tvar, 'yminor', 10
        options, tvar, 'ytitle', 'Freq (Hz)'

    ;----E/B spectrogram.
        vars = pre0+['de_fac','db_fac']
        ndim = 3

        ; settings for wavelet transform.
        s0 = 4d*dr0
        dj = 1d/8
        s1 = 2000
        j1 = floor(alog(s1/s0)/alog(2)/dj)
        s1 = s0*2d^(dj*j1)
        ns = j1+1
        w0 = 6d
        cdelta = 0.776d
        psi0 = !dpi^(-0.25)

        foreach tvar, vars do begin
            get_data, tvar, uts, dat
            dat = snorm(dat)
            mor = wavelet(dat, dr0, /pad, $
                s0=s0, dj=dj, j=j1, mother='Morlet', param=w0, period=ps, scale=ss)
            psd = abs(mor)^2
            idx = where(uts ge plot_time_range[0] and uts le plot_time_range[1], tnrec)
            psd = psd[idx,*]
            gws = total(psd,1)/tnrec^2
            ngws = (gws/ss)*(dr0*dj/cdelta)*tnrec
            store_data, tvar+'_tmp', ps, [[gws],[ngws]]
        endforeach
    endforeach


    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    thick = (size(plot_file,/type) eq 7)? 100: 4

;--plot pflux 3d.
    poss = sgcalcpos(5,position=posu)
    foreach pre0, 'rb'+probes+'_' do begin
        if pre0 eq 'rba_' then begin
            tvar = pre0+'pf_para_mor_spec_tmp'
            options, tvar, 'ztitle', '( '+pfunit+')'

            tvar = pre0+'de_fac'
            options, tvar, 'ystyle', 1
            options, tvar, 'ytitle', '(mV/m)'
            options, tvar, 'yrange', [-1,1]*70
            options, tvar, 'ytickv', [-1,0,1]*50
            options, tvar, 'yticks', 2
            options, tvar, 'yminor', 5
            options, tvar, 'constant', 0
            options, tvar, 'labels', 'dE!D'+labfac


            tvar = pre0+'db_fac'
            options, tvar, 'ystyle', 1
            options, tvar, 'ytitle', '(nT)'
            options, tvar, 'yrange', [-1,1]*18
            options, tvar, 'ytickv', [-1,0,1]*15
            options, tvar, 'yticks', 2
            options, tvar, 'yminor', 5
            options, tvar, 'constant', 0
            options, tvar, 'labels', 'dB!D'+labfac


            tvar = pre0+'n_combine'
            options, tvar, 'ytitle', '(cm!U-3!N)'
            options, tvar, 'labels', ['n!DVsc','n!Dele']

            tvar = pre0+'pf_fac_tmp'
            get_data, pre0+'pf_fac_tmp', times, data
            cmap = get_var_data(pre0+'cmap', in=plot_time_range, limits=lim)
            model_index = where(lim.labels eq default_model)
            cmap0 = mean(cmap[*,model_index])
            store_data, tvar, times, data*cmap0, limits={$
                ytitle:'Map@!C100km!C(mW/m!U2!N)', colors:rgb, cmap0:cmap0}

            options, tvar, 'yticks', 2
            options, tvar, 'yminor', 5
            options, tvar, 'ystyle', 1
            options, tvar, 'ytickv', [0,1,2]*25
            options, tvar, 'yrange', [-1,4]*12.5
        endif else begin
            tvar = pre0+'pf_para_mor_spec_tmp'
            options, tvar, 'ztitle', '( '+pfunit+')'

            tvar = pre0+'de_fac'
            options, tvar, 'ystyle', 1
            options, tvar, 'ytitle', '(mV/m)'
            options, tvar, 'yrange', [-1,1]*120
            options, tvar, 'ytickv', [-1,0,1]*120
            options, tvar, 'yticks', 2
            options, tvar, 'yminor', 5
            options, tvar, 'constant', 0
            options, tvar, 'labels', 'dE!D'+labfac


            tvar = pre0+'db_fac'
            options, tvar, 'ystyle', 1
            options, tvar, 'ytitle', '(nT)'
            options, tvar, 'yrange', [-1,1]*30
            options, tvar, 'ytickv', [-1,0,1]*30
            options, tvar, 'yticks', 2
            options, tvar, 'yminor', 5
            options, tvar, 'constant', 0
            options, tvar, 'labels', 'dB!D'+labfac


            tvar = pre0+'n_combine'
            options, tvar, 'ytitle', '(cm!U-3!N)'
            options, tvar, 'labels', ['n!DVsc','n!Dele']

            tvar = pre0+'pf_fac_tmp'
            get_data, pre0+'pf_fac', times, data
            get_data, pre0+'map_coef1', tuts, cmap
            cmap = get_var_data(pre0+'cmap', in=plot_time_range, limits=lim)
            model_index = where(lim.labels eq default_model)
            cmap0 = mean(cmap[*,model_index])
            store_data, tvar, times, data*cmap0, limits={$
                ytitle:'Map@!C100km!C(mW/m!U2!N)', colors:rgb, cmap0:cmap0}

            options, tvar, 'yticks', 2
            options, tvar, 'yminor', 5
            options, tvar, 'ystyle', 1
            options, tvar, 'ytickv', [0,1,2]*80
            options, tvar, 'yrange', [-1,4]*40
        endelse
    endforeach

    panel_label_xshift = 11

    foreach probe, probes do begin
        pre0 = 'rb'+probe+'_'
        label_pre = probe
        the_poss = (probe eq 'a')? upper_left_poss: upper_right_poss

        vars = pre0+['n_combine','de_fac','db_fac','pf_fac_tmp']
        figlabs = label_pre+['-1. Density', '-2. dE FAC', '-3. dB FAC', '-4. S FAC']
        poss = the_poss[*,0:nypanel-2]
        tplot, vars, position=poss, /noerase, /nouttick, trange=plot_time_range
        for i=0,nypanel-2 do xyouts, poss[0,i]-xchsz*panel_label_xshift, poss[3,i]-ychsz*0.5, /normal, figlabs[i]
        tpos = poss[*,3]
        xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1.2, /normal, 'Filtered in '+sgnum2str(1d3/filter[1],ndec=0)+'mHz-'+sgnum2str(1d/filter[0],ndec=0)+'Hz'
        labels = 'S!D'+labfac
        colors = rgb
        lengths = [0,2,6]-2
        for i=0,2 do xyouts, tpos[0]+xchsz*(20+lengths[i]), tpos[3]-ychsz*1.2, /normal, labels[i], color=colors[i]
        tvar = pre0+'pf_fac_tmp'
        get_data, tvar, limits=lims
        plot, [0,1], [0,1], position=poss[*,3], /nodata, xstyle=5, ystyle=5, /noerase
        axis, yaxis=1, yrange=lims.yrange/lims.cmap0, ystyle=1, ytitle='In-situ!C(mW/m!U2!N)', yticks=lims.yticks, ytickv=lims.ytickv/lims.cmap0, yminor=lims.yminor, yticklen=-0.01, ytickformat='(F3.1)'

        device, decomposed=0
        loadct2, 66
        vars = pre0+'pf_para_mor_spec_tmp'
        tpos = the_poss[*,nypanel-1]
        tplot, vars, position=tpos, /noerase, trange=plot_time_range, /novtitle
        xyouts, tpos[0]-xchsz*panel_label_xshift, tpos[3]-ychsz*0.5, /normal, label_pre+'-5. S!D||!N'
        device, decomposed=1

    ;---Lower panels.
        red = sgcolor('red')
        poss = (probe eq 'a')? lower_poss[*,0:2]: lower_poss[*,3:5]
        ;for ii=0,2 do plot, [0,1],[0,1], position=the_poss[*,ii], /nodata, /noerase

        vars = pre0+['de_fac','db_fac']
        get_data, vars[0]+'_tmp', ps, edat
        get_data, vars[1]+'_tmp', ps, bdat

        yrange = [0.1,10000]
        ytickv = 10^smkarthm(-1,4,1,'dx')
        ytickn = '10!U'+string(alog10(ytickv),format='(I0)') & ytickn[1:*:2]=' '
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytitle = 'Period (sec)'

        ; convert to frequency.
        ps = 1d/ps
        yrange = minmax(1d/[0.1,10000])
        ytickv = 10^smkarthm(-4,1,1,'dx')
        ytickn = '10!U'+string(alog10(ytickv),format='(I0)') & ytickn[1:*:2] = ' '
        yticks = n_elements(ytickv)-1
        yminor = 10
        ytitle = 'Freq (Hz)'

        ; E spectrogram.
        tpos = poss[*,0]
        xrange = 5*[1e-5,10]
        xtickv = 10^smkarthm(-4,1,1,'dx')
        xtickn = '10!U'+string(alog10(xtickv),format='(I0)') & xtickn[1:*:2]=' '
        xticks = n_elements(xtickv)-1
        xminor = 10
        xtitle = '|E| PS (mV/m)!U2!N'
        plot, edat[*,0], ps, /noerase, position=tpos, $
            /xlog, xstyle=1, xrange=xrange, xtitle=xtitle, $
            xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xticklen=ticklen, $
            /ylog, ystyle=1, yrange=yrange, ytitle=ytitle, $
            ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, yticklen=ticklen
        xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.25, /normal, label_pre+'-6.'
        ; add line fit.
        idx = where(ps ge min(ps), tns)
        tps = ps[idx]
        res = linfit(alog10(tps), alog10(edat[idx,0]))
        tys = 10^(alog10(tps)*res[1]+res[0])
        plots, tys, ps[idx], color=red, linestyle=0
        tmp = convert_coord(tps[tns/2], tys[tns/2], /data, /to_normal)
        xyouts, tmp[0]+xchsz*1, tmp[1]+ychsz*0.5, color=red, /normal, 'f!U '+sgnum2str(res[1],ndec=2)


        ; B spectrogram.
        tpos = poss[*,1]
        xrange = 5*[1e-8,10]
        xtickv = 10^smkarthm(-6,1,1,'dx')
        xtickn = '10!U'+string(alog10(xtickv),format='(I0)') & xtickn[1:*:2]=' '
        xticks = n_elements(xtickv)-1
        xminor = 10
        xtitle = '|B| PS (nT)!U2!N'
        plot, bdat[*,0], ps, /noerase, position=tpos, $
            /xlog, xstyle=1, xrange=xrange, xtitle=xtitle, $
            xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xticklen=ticklen, $
            /ylog, ystyle=1, yrange=yrange, ytitle='', $
            ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, yticklen=ticklen, ytickformat='(A1)'
        xyouts, tpos[0]+xchsz*0.5, tpos[1]+ychsz*0.25, /normal, label_pre+'-7.'
        ; add line fit.
        idx = where(ps ge min(ps), tns)
        tps = ps[idx]
        res = linfit(alog10(tps), alog10(bdat[idx,0]))
        tys = 10^(alog10(tps)*res[1]+res[0])
        plots, tys, ps[idx], color=red, linestyle=0
        tmp = convert_coord(tps[tns/2], tys[tns/2], /data, /to_normal)
        xyouts, tmp[0]+xchsz*1, tmp[1]+ychsz*0.5, color=red, /normal, 'f!U '+sgnum2str(res[1],ndec=2)


        ; E/B in km/s.
        tpos = poss[*,2]
        xrange = 5*[1e2,1e4]
        xtickv = 10^smkarthm(3,4,1,'dx')
        xtickn = '10!U'+string(alog10(xtickv),format='(I0)')
        xticks = n_elements(xtickv)-1
        xminor = 10
        xtitle = 'E/B (km/s)'
        xyouts, tpos[2]-xchsz*0.5, tpos[1]+ychsz*0.25, /normal, label_pre+'-8.', alignment=1

        plot, xrange, yrange, /nodata, /noerase, position=tpos, $
            /xlog, xstyle=5, xrange=xrange, xtickformat='(A1)', $
            /ylog, ystyle=5, yrange=yrange, ytickformat='(A1)'

        mu0 = 4*!dpi*1e-7
        sigma_pedersens = [1,10]
        vpedersens = 1e-3/(mu0*sigma_pedersens)
        vplot = minmax(vpedersens>xrange[0])
        polyfill, vplot[[0,0,1,1,0]],yrange[[0,1,1,0,0]], color=fill_color
        tmp1 = convert_coord(vplot[0],yrange[1], /data, /to_normal)
        tmp2 = convert_coord(vplot[1],yrange[1], /data, /to_normal)
        ty = tpos[3]+ychsz*0.3
        tx = mean([tmp1[0],tmp2[0]])
        xyouts, tx, ty, /normal, '1/'+strmu+'!D0!N'+strsigma+'!DP!N', alignment=0.5


        event_info = get_var_data(id)
        the_info = (probe eq 'a')? event_info.rbspa: event_info.rbspb
        vplot = the_info.va*[1d,1/sqrt(3)]
        polyfill, vplot[[0,0,1,1,0]],yrange[[0,1,1,0,0]], color=fill_color
        tmp1 = convert_coord(vplot[0],yrange[1], /data, /to_normal)
        tmp2 = convert_coord(vplot[1],yrange[1], /data, /to_normal)
        ty = tpos[3]+ychsz*0.3
        tx = mean([tmp1[0],tmp2[0]])
        xyouts, tx, ty, /normal, 'v!DA!N', alignment=0.5

        plot, sqrt(edat[*,0]/bdat[*,0])*1e3, ps, /noerase, position=tpos, $
            /xlog, xstyle=1, xrange=xrange, xtitle=xtitle, $
            xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, xticklen=ticklen, $
            /ylog, ystyle=1, yrange=yrange, ytitle='', $
            ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn, yticklen=ticklen, ytickformat='(A1)'
    endforeach


    sgclose

end
