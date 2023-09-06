;+
; Pflux detais on E, B, S, E/B ratio, ewogram.
;
; Test referee's suggestion on using Eo and Bw.
;-


    event_info = _2015_0218_02_load_data()


;---Settings.
    time_range = event_info['time_range']
    probe = event_info['probe']
    prefix = event_info['prefix']
    dr0 = event_info['field_time_step']
    e_var = prefix+'edot0_fac'
    b_var = prefix+'b1_fac'
    pf_var = prefix+'pfdot0_fac'
test = 1

    perp = '!9'+string(94b)+'!X'
    labfac = ['||',perp+',west',perp+',out']
    rgb = constant('rgb')
    xyz = ['x','y','z']
    label_size = 0.8

    ticklen = -0.02
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 4
    tplot_options, 'labflag', -1
    tplot_options, 'ticklen', ticklen
    tplot_options, 'zticklen', -0.1
    tplot_options, 'thick', 1



;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1

    strmu = '!9'+string(109b)+'!X'
    strsigma = '!9'+string(83b)+'!X'
    red = sgcolor('red')
    blue = sgcolor('blue')



    filter = event_info['pflux_filter']
    scaleinfo = event_info['pflux_scale_info']
    copy_data, pf_var, prefix+'pf_fac_tmp'


;---EWOgram and KEOgram.
    mlt_range = [-2.5,-0.5]+24
    mlat_range = [61.5,64.5]
    zrange = [60,180]
    zstep = 30
    ztickv = make_bins(zrange,zstep)
    zticks = n_elements(ztickv)-1
    zminor = 6
    asi_ct = 33

    ewo_var = 'thg_asf_ewo'
    get_data, ewo_var, uts, data, val
    index = where(val le 0, count)
    if count ne 0 then store_data, ewo_var, uts, data, val+24
    yrange = mlt_range
    ystep = 1
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 5
    options, ewo_var, 'yrange', yrange
    options, ewo_var, 'yticks', yticks
    options, ewo_var, 'ytickv', ytickv
    options, ewo_var, 'yminor', yminor
    options, ewo_var, 'zrange', zrange
    options, ewo_var, 'zticks', zticks
    options, ewo_var, 'ytitle', 'MLT!C(h)'
    options, ewo_var, 'color_table', asi_ct

    keo_var = 'thg_asf_keo'
    yrange = mlat_range
    ystep = 1
    ytickv = make_bins(yrange, ystep, /inner)
    yticks = n_elements(ytickv)-1
    yminor = 4
    options, keo_var, 'yrange', yrange
    options, keo_var, 'yticks', yticks
    options, keo_var, 'ytickv', ytickv
    options, keo_var, 'yminor', yminor
    options, keo_var, 'zrange', zrange
    options, keo_var, 'zticks', zticks
    options, keo_var, 'ytitle', 'MLat!C(deg)'



;---pflux spectrogram, normalized version.
    ; change unit to (uW/m^2), convert y from period to freq
    tvar = pf_var+'_mor_spec_1'
    get_data, tvar, uts, dat, val, limits = lim
    dat *= 1e3
    zr = [-1,1]*0.4
    pfunit = strmu+'W/m!U2!N'
    tvar = prefix+'pf_para_mor_spec_tmp'
    store_data, tvar, uts, dat, val, limits = lim
    options, tvar, 'zrange', zr
    options, tvar, 'ztitle', 'S!D'+labfac[0]+'!N ('+pfunit+')'
    options, tvar, 'zticks', 2
    options, tvar, 'yrange', [1e-3,4]
    options, tvar, 'ytickv', [1e-2,1e-1,1e0]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 10
    options, tvar, 'ytitle', 'Freq!C(Hz)'
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    options, tvar, 'ylog', 1
    options, tvar, 'no_color_scale', 0


;----E/B spectrogram.
    options, e_var, 'labels', 'dE!D'+labfac
    options, e_var, 'colors', rgb
    options, b_var, 'labels', 'dB!D'+labfac
    options, b_var, 'colors', rgb
    fac_text_labels = ['para','west','out']
    stplot_split, e_var, newnames=e_var+'_'+fac_text_labels
    stplot_split, b_var, newnames=b_var+'_'+fac_text_labels
    spec_vars = [e_var+'_out',b_var+'_west']
;    spec_vars = [e_var+'_west',b_var+'_out']
    spec_vars = [e_var,b_var]+'_mag'
    foreach var, [e_var,b_var] do begin
        get_data, var, uts, vec
        store_data, var+'_mag', uts, snorm(vec[*,1:2])
    endforeach
    ndim = 3

    ; settings for wavelet transform.
    s0 = 4d*dr0
    dj = 1d/8
    s1 = 2000
s1 = 600d
    j1 = floor(alog(s1/s0)/alog(2)/dj)
    s1 = s0*2d^(dj*j1)
    ns = j1+1
    w0 = 6d
    cdelta = 0.776d
    psi0 = !dpi^(-0.25)

    suff = '_'+[fac_text_labels,'mag']
    all_spec_vars = [e_var+suff,b_var+suff]
    foreach tvar, all_spec_vars do begin
        get_data, tvar, uts, dat
        index = where_pro(uts, '[]', time_range)
        uts = uts[index]
        dat = dat[index,*]
;        dat = snorm(dat)
        index = where(finite(dat,/nan), count)
        if count ne 0 then dat[index] = 0
        mor = wavelet(dat, dr0, /pad, $
            s0=s0, dj=dj, j=j1, mother='Morlet', param=w0, period=ps, scale=ss)
        psd = abs(mor)^2
        idx = where(uts ge time_range[0] and uts le time_range[1], tnrec)
        psd = psd[idx,*]
        gws = total(psd,1)/tnrec^2
        ngws = (gws/ss)*(dr0*dj/cdelta)*tnrec
        store_data, tvar+'_tmp', ps, [[gws],[ngws]]
    endforeach


    plot_file = join_path([srootdir(),'2015_0218_fig_pflux_detail2.pdf'])
    if test eq 1 then plot_file = 0

    sgopen, plot_file, xsize=5, ysize=6, /inch

    posu = [0.25,0.28,0.85,0.98]
    posl = [0.15,0.08,0.92,0.20]

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    thick = (size(plot_file,/type) eq 7)? 100: 4

;--plot pflux 3d.
    tvar = prefix+'efw_density'
    options, tvar, 'yrange', [1,20]
    options, tvar, 'ylog', 1
    options, tvar, 'yminor', 10
    options, tvar, 'labels', 'From!C  UH line'

    tvar = prefix+'pf_para_mor_spec_tmp'
    options, tvar, 'ztitle', '( '+pfunit+')'
    options, tvar, 'color_table', 66


    tvar = e_var+['0','1','2']
    stplot_split, e_var, newnames=tvar, colors=rgb, labels='dE!D'+labfac
    options, tvar, 'ystyle', 1
    options, tvar, 'ytitle', '(mV/m)'
    options, tvar, 'yrange', [-1,1]*20
    options, tvar, 'ytickv', [-1,0,1]*15
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 3
    options, tvar, 'constant', 0

    tvar = b_var
    options, tvar, 'ystyle', 1
    options, tvar, 'ytitle', '(nT)'
    options, tvar, 'yrange', [-1,1]*20
    options, tvar, 'ytickv', [-1,0,1]*15
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 3
    options, tvar, 'constant', 0
    options, tvar, 'labels', 'dB!D'+labfac
    options, tvar, 'colors', rgb


    tvar = prefix+'pf_fac_tmp'
    get_data, prefix+'pf_fac_tmp', times, data
    get_data, prefix+'cmap', tuts, cmap
    cmap0 = mean(cmap[where(tuts ge time_range[0] and tuts le time_range[1])])
    store_data, tvar, times, data*cmap0, limits={$
        ytitle:'Map@!C100km!C(mW/m!U2!N)', colors:rgb}

    options, tvar, 'ystyle', 1
    options, tvar, 'yrange', [-15,5]
    options, tvar, 'ytickv', [-3,-1,1]*5
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'constant', 0




    ; First batch of vars.
    vars = [[ewo_var],e_var+['1','2'],prefix+['pf_fac_tmp']]
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar+1,position=posu, ypans=[2,fltarr(nvar)+1])
    letters = letters(nvar+2)
    figlabs = letters[0:nvar-1]+') '+['Ewogram', 'dE west', 'dE out', 'S']


    tpos = poss[*,0:nvar-1]

    ; Use one colarbar for ewo and keo.
;    cbpos = [tpos[2,0],tpos[1,1],tpos[2,0],tpos[3,0]]
    cbpos = tpos[*,0]
    cbpos[0] = cbpos[2]+xchsz*1
    cbpos[2] = cbpos[0]+xchsz*1
    options, ewo_var, 'zposition', cbpos
    options, keo_var, 'no_color_scale', 1
    options, tnames('*'), 'charsize', 1


    ; Plot the first batch of vars.
    tplot, vars, position=tpos, /noerase, /nouttick, trange=time_range
    for i=0,nvar-1 do xyouts, tpos[0,i]-xchsz*12, tpos[3,i]-ychsz*0.5, /normal, figlabs[i]

    ; Add fpt to keo and ewo.
    model_setting = event_info['model_setting']
    model = model_setting['model']
;    foreach var, [keo_var,ewo_var] do begin
;        get_data, var, limits=lim
;        xrange = time_range
;        yrange = lim.yrange
;        tpos = reform(poss[*,where(vars eq var)])
;        plot, xrange, yrange, $
;            xstyle=5, xrange=xrange, $
;            ystyle=5, yrange=yrange, $
;            position=tpos, nodata=1, noerase=1
;
;        pos_var = (var eq ewo_var)? prefix+'fmlt_'+model: prefix+'fmlat_'+model
;        get_data, pos_var, xxs, yys
;        if var eq ewo_var then yys += 24
;        oplot, xxs, yys, color=sgcolor('black'), linestyle=2
;    endforeach

    ; Add lines to ewogram.
    tpos = poss[*,where(vars eq ewo_var)]
    get_data, ewo_var, limits=lim
    xrange = time_range
    yrange = lim.yrange
    plot, xrange, yrange, $
        xstyle=5, xrange=xrange, $
        ystyle=5, yrange=yrange, $
        position=tpos, nodata=1, noerase=1

    ; Add footpoint.
    fpt_color = sgcolor('red')
    pos_var = prefix+'fmlt_'+model
    get_data, pos_var, xxs, yys
    index = where(yys le 0, count)
    if count ne 0 then yys[index] += 24
    oplot, xxs, yys, linestyle=3, color=fpt_color
    tx = tpos[2]-xchsz*10
    ty = tpos[1]+ychsz*2
    xyouts, tx,ty,normal=1, 'RBSP-A';, color=fpt_color

    east_times = time_double('2015-02-18/'+['02:07','02:11'])
    east_mlts = [22.7,21.9]
    east_color = sgcolor('white')
    oplot, east_times, east_mlts, linestyle=3, color=east_color
    east_speed = total(east_mlts*[-1,1])/total(east_times*[-1,1])*15*60
    msg = tex2str('Omega')+' west!N!C'+strtrim(string(east_speed,format='(F5.1)'),2)+' deg/min'
    tx = tpos[0]+xchsz*4
    ty = tpos[3]-ychsz*4
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=east_color

    west_times = time_double('2015-02-18/'+['02:09:40','02:11:40'])
    west_mlts = [22.25,22.5]-0.07
;    west_mlts = [22.3,22.7]
    west_color = sgcolor('black')
 ;   west_color = sgcolor('purple')
    oplot, west_times, west_mlts, linestyle=3, color=west_color
    west_speed = total(west_mlts*[-1,1])/total(west_times*[-1,1])*15*60
    msg = tex2str('Omega')+' east2!N!C'+strtrim(string(west_speed,format='(F5.1)'),2)+' deg/min'
    tmp = convert_coord(west_times[-1],west_mlts[-1], data=1, to_normal=1)
    tx = tmp[0]+xchsz*5
    ty = tmp[1]+ychsz*0.5
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=west_color

    west_times = time_double('2015-02-18/'+['02:07:20','02:08:50'])-10
    west_mlts = [22.7,23.2]
    oplot, west_times, west_mlts, linestyle=3, color=west_color
    west_speed = total(west_mlts*[-1,1])/total(west_times*[-1,1])*15*60
    msg = tex2str('Omega')+' east1!N!C'+strtrim(string(west_speed,format='(F5.1)'),2)+' deg/min'
    tmp = convert_coord(west_times[-1],west_mlts[-1], data=1, to_normal=1)
    tx = tmp[0]+xchsz*5
    ty = tmp[1]+ychsz*0.5
    xyouts, tx,ty,normal=1, msg, alignment=0.5, color=sgcolor('white')


    ; Add labels for pflux.
    tpos = poss[*,nvar-1]
    xyouts, tpos[0]+xchsz*1, tpos[1]+ychsz*0.5, /normal, 'Filtered in '+sgnum2str(1d3/filter[1],ndec=0)+'mHz-'+sgnum2str(1d/filter[0],ndec=0)+'Hz'
    colors = rgb
    labels = 'S!D'+labfac
    lengths = [0,3,8]
    for i=0,2 do xyouts, tpos[0]+xchsz*(22+lengths[i]), tpos[1]+ychsz*0.5, /normal, labels[i], color=colors[i]
    tvar = prefix+'pf_fac_tmp'
    get_data, tvar, limits=lims
    plot, time_range, lims.yrange, position=poss[*,nvar-1], /nodata, xstyle=5, ystyle=5, /noerase
    tmp = convert_coord(time_range[1],0, /data, /to_normal)
    tx = tmp[0]-xchsz*0.5
    ty = tmp[1]
    xyouts, tx,ty+ychsz*0.2,/normal, 'To N-hem', alignment=1, charsize=label_size
    xyouts, tx,ty-ychsz*1.0,/normal, 'To S-hem', alignment=1, charsize=label_size
    axis, yaxis=1, yrange=lims.yrange/cmap0*1e3, ystyle=1, $
        ytitle='In-situ!C('+pfunit+')', yticks=lims.yticks, $
        ytickv=lims.ytickv/cmap0*1e3, yminor=lims.yminor, yticklen=-0.01, ytickformat='(F5.1)'

    ; pflux spec.
    spec_var = prefix+'pf_para_mor_spec_tmp'
    tpos = poss[*,nvar]
    cbpos = tpos
    cbpos[0] = cbpos[2]+xchsz*1
    cbpos[2] = cbpos[0]+xchsz*1
    options, spec_var, 'zposition', cbpos
    tplot, spec_var, position=tpos, /noerase, trange=time_range, /novtitle
    xyouts, tpos[0]-xchsz*12, tpos[3]-ychsz*0.5, /normal, letters[nvar]+') S!D||!N PSD'
    get_data, spec_var, limits=lim
    plot, time_range, lim.yrange, $
        xstyle=5, $
        ystyle=5, ylog=1, $
        nodata=1, noerase=1, position=tpos
    fb = 0.015  ; hz.
    oplot, time_range, fb+[0,0], linestyle=1
    tx = tpos[0]+xchsz*0.5
    ty = (convert_coord(0,fb, data=1, to_normal=1))[1]+ychsz*0.3
    msg = 'f!Dbead!N = '+string(fb,format='(F5.3)')+' Hz'
;    xyouts, tx,ty, normal=1, msg;, charsize=label_size
stop



;--plot spectrums.
    poss = sgcalcpos(1,3, position=posl, xpad=8)
    xticklen = -0.02
    yticklen = -0.02


;    xrange = [0.1,10000]
;    xtickv = 10^smkarthm(-1,4,1,'dx')
;    xtickn = '10!U'+string(alog10(xtickv),format='(I0)') & xtickn[1:*:2]=' '
;    xticks = n_elements(xtickv)-1
;    xminor = 10
;    xtitle = 'Period (sec)'


    ; convert to frequency.
    get_data, e_var+'_mag_tmp', ps
    fs = 1d/ps
    xrange = minmax(1d/[0.1,10000])
    xtickv = 10^smkarthm(-4,1,1,'dx')
    xtickn = '10!U'+string(alog10(xtickv),format='(I0)') & xtickn[1:*:2] = ' '
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtitle = 'Freq (Hz)'

    ; E spectrogram.
    tpos = poss[*,0]
    yrange = 5*[1e-6,0.1]
    ytickv = 10^smkarthm(-4,1,1,'dx')
    yrange = minmax([yrange,ytickv])
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)') & ytickn[1:*:2]=' '
    index = where(ytickn eq '10!U0', count)
    if count ne 0 then ytickn[index] = '1'
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytitle = '(mV/m)!U2!N'

    ; Draw power spectrum.
    plot, xrange, yrange, $
        xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, xtickname=xtickn, $
        xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, ylog=1, yrange=yrange, ytitle=ytitle, ytickname=ytickn, $
        yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
        noerase=1, position=tpos, nodata=1

    f_spin = 2d/10.8
    plots, f_spin+[0,0], yrange, linestyle=1
    tmp = convert_coord(f_spin, yrange[0], data=1, to_normal=1)
    tx = tmp[0]-xchsz*0.2
    ty = tmp[1]+ychsz*0.35
    msg = '2f!Dspin!N'
    xyouts, tx,ty,normal=1, msg, alignment=0.25;, charsize=label_size
    plots, fb+[0,0], yrange, linestyle=1

    e_vars = e_var+['_west','_out','_mag']+'_tmp'
    labels = ['E!D'+labfac[[1,2]],'|E|']
    e_vars = e_var+['_out','_mag']+'_tmp'
    labels = ['E!D'+labfac[2],'|E|']
    e_vars = e_var+['_out','_west']+'_tmp'
    labels = ['E!D'+labfac[2],'E!D'+labfac[1]]
    colors = sgcolor(['black','green','red'])
    foreach var, e_vars, var_id do begin
        get_data, var, ps, dat
        plots, fs, dat, color=colors[var_id]
        tx = tpos[2]-xchsz*4
        ty = tpos[3]-ychsz*(var_id+0.9)
        msg = labels[var_id]
        xyouts, tx,ty,normal=1, msg, color=colors[var_id];, charsize=label_size
    endforeach

    ; Add label.
    xyouts, tpos[0]-xchsz*7, tpos[3]-ychsz*0.3, /normal, $
        letters[-1]+'-1) E PS'


;    ; add line fit.
;    idx = where(ps ge min(ps), tns)
;    tps = ps[idx]
;    res = linfit(alog10(tps), alog10(edat[idx,0]))
;    tys = 10^(alog10(tps)*res[1]+res[0])
;    plots, ps[idx], tys, color=red, linestyle=0
;    tmp = convert_coord(tps[tns/2], tys[tns/2], /data, /to_normal)
;    xyouts, tmp[0]+xchsz*1, tmp[1]+ychsz*0.5, color=red, /normal, 'f!U '+sgnum2str(res[1],ndec=2)


    ; B spectrogram.
    tpos = poss[*,1]
    yrange = 5*[1e-8,10]
    ytickv = 10^smkarthm(-6,1,1,'dx')
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)') & ytickn[1:*:2]=' '
    index = where(ytickn eq '10!U0', count)
    if count ne 0 then ytickn[index] = '1'
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytitle = '(nT)!U2!N'

    ; Draw power spectrum.
    plot, xrange, yrange, $
        xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, xtickname=xtickn, $
        xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, ylog=1, yrange=yrange, ytitle=ytitle, ytickname=ytickn, $
        yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
        noerase=1, position=tpos, nodata=1

    b_vars = b_var+['_out','_west','_mag']+'_tmp'
    labels = ['B!D'+labfac[[2,1]],'|B|']
    b_vars = b_var+['_west','_mag']+'_tmp'
    labels = ['B!D'+labfac[[1]],'|B|']
    b_vars = b_var+['_west']+'_tmp'
    labels = ['B!D'+labfac[[1]]]
    foreach var, b_vars, var_id do begin
        get_data, var, ps, dat
        plots, fs, dat, color=colors[var_id]
        tx = tpos[2]-xchsz*4
        ty = tpos[3]-ychsz*(var_id+0.9)
        msg = labels[var_id]
        xyouts, tx,ty,normal=1, msg, color=colors[var_id];, charsize=label_size
    endforeach

    plots, fb+[0,0], yrange, linestyle=1

    ; Add label.
    xyouts, tpos[0]-xchsz*6, tpos[3]-ychsz*0.3, /normal, $
        letters[-1]+'-2) B PS'

;    ; add line fit.
;    idx = where(ps ge min(ps), tns)
;    tps = ps[idx]
;    res = linfit(alog10(tps), alog10(bdat[idx,0]))
;    tys = 10^(alog10(tps)*res[1]+res[0])
;    plots, ps[idx], tys, color=red, linestyle=0
;    tmp = convert_coord(tps[tns/2], tys[tns/2], /data, /to_normal)
;    xyouts, tmp[0]+xchsz*1, tmp[1]+ychsz*0.5, color=red, /normal, 'f!U '+sgnum2str(res[1],ndec=2)


    ; E/B in km/s.
    tpos = poss[*,2]
    yrange = 1*[50,5e4]
    ytickv = 10^smkarthm(2,4,1,'dx')
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)')
 ;   ytickn[-1] = ' '
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytitle = '(km/s)'

    ; Draw power spectrum.
    plot, xrange, yrange, $
        xstyle=5, xlog=1, xrange=xrange, xtitle=xtitle, xtickname=xtickn, $
        xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=5, ylog=1, yrange=yrange, ytitle=ytitle, ytickname=ytickn, $
        yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
        noerase=1, position=tpos, nodata=1

    f_spin = 2d/10.8
    plots, f_spin+[0,0], yrange, linestyle=1

;    ; Add Pedersen conductivity.
;    mu0 = 4*!dpi*1e-7
;    sigma_pedersens = [1,10]
;    vpedersens = 1e-3/(mu0*sigma_pedersens)
;    vplot = minmax(vpedersens>yrange[0])
;    polyfill, xrange[[0,1,1,0,0]],vplot[[0,0,1,1,0]], color=sgcolor('silver')
;    tmp1 = convert_coord(xrange[1],vplot[0], /data, /to_normal)
;    tmp2 = convert_coord(xrange[1],vplot[1], /data, /to_normal)
;    ty = ((tmp1[1]+tmp2[1])*0.5-ychsz*0.5)>tmp1[1]
;    tx = tmp1[0]+xchsz*0.5
;    xyouts, tx, ty, /normal, '1/'+strmu+'!D0!N'+strsigma+'!DP!N!C1-10 S'


    labels = ['E!D'+labfac[2]+'!N/B!D'+labfac[1],'|E|/|B|']
    for var_id=0,0 do begin
        e_var = e_vars[var_id]
        b_var = b_vars[var_id]
        get_data, e_var, ps, edat
        get_data, b_var, ps, bdat
        plots, fs, sqrt(edat[*,0]/bdat[*,0])*1e3, color=colors[var_id]
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*(var_id+0.9)
        msg = labels[var_id]
        xyouts, tx,ty,normal=1, msg, color=colors[var_id];, charsize=label_size
    endfor

    ; Draw power spectrum.
    plot, xrange, yrange, $
        xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, xtickname=xtickn, $
        xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
        ystyle=1, ylog=1, yrange=yrange, ytitle=ytitle, ytickname=ytickn, $
        yticklen=yticklen, ytickv=ytickv, yticks=yticks, yminor=yminor, $
        noerase=1, position=tpos, nodata=1

    ; Add label.
    xyouts, tpos[0]-xchsz*6, tpos[3]-ychsz*0.3, /normal, letters[-1]+'-3) E/B'



    ;foreach suf, ['','_p','_o'] do begin
    foreach suf, [''] do begin
        va = event_info['va'+suf]
        omega_i = event_info['fg'+suf]
        vi = event_info['vi'+suf]
        vf = event_info['vf'+suf]
        ebr = va*sqrt(1+(fs/omega_i*(vi/vf))^2)
        plots, fs, ebr, color=red

        index = n_elements(fs)*0.2
        tmp = convert_coord(fs[index], ebr[index], /data, /to_normal)
        xyouts, tmp[0]-xchsz*1.5, tmp[1]-ychsz*1.8, color=red, /normal, 'KAW'
    endforeach

;    plots, fb+[0,0], yrange, linestyle=1
    ;plots, minmax(ps), vao+[0,0], color=blue, linestyle=1

    ;plots, 1d/filter, vao*1.05, color=sgcolor('red')



    sgclose

end
