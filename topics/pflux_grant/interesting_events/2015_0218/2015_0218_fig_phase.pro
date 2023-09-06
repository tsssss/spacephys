;+
; Check the phase between E and B.
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
        store_data, tvar+'_mor', uts, mor, ps
    endforeach


    plot_file = join_path([srootdir(),'2015_0218_fig_phase.pdf'])
    if test eq 1 then plot_file = 0

    sgopen, plot_file, xsize=5, ysize=3, /inch, magnify=2

    posu = [0.25,0.28,0.85,0.98]
    posl = [0.15,0.08,0.92,0.20]

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size

    thick = (size(plot_file,/type) eq 7)? 100: 4


    ; Phase.
    get_data, e_var+'_out_mor', times, emor, ps
    get_data, b_var+'_west_mor', times, bmor, ps
    wxy_nj = emor*conj(bmor)
    cp_nj = atan(imaginary(wxy_nj), real_part(wxy_nj))*deg
    wave_time_range = time_double(['2015-02-18/02:08','2015-02-18/02:11'])
;    wave_time_range = time_double(['2015-02-18/02:12','2015-02-18/02:15'])
    index = where_pro(times, '[]', wave_time_range, count=count)
    avg_phase = total(cp_nj[index,*],1)/count

    fs = 1d/ps
    xrange = minmax(1d/[0.1,10000])
    xrange = [1e-3,4]
    xtickv = 10^smkarthm(-2,0,1,'dx')
    xtickn = '10!U'+string(alog10(xtickv),format='(I0)') & xtickn[1:*:2] = ' '
    xtickn = ['0.01','0.1','1']
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtitle = 'Freq (Hz)'

    yrange = [-1,1]*180
    ystep = 90d
    ytickv = make_bins(yrange,ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = 3
    ytitle = 'Phase (deg)'

    sgopen, 0, xsize=3, ysize=3
    tpos = sgcalcpos(1, margins=[8,4,2,1])
;    plot, xrange, yrange, $
;        xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, xtickname=xtickn, $
;        xticklen=xticklen, xtickv=xtickv, xticks=xticks, xminor=xminor, $
;        ystyle=1, ylog=0, yrange=yrange, ytitle=ytitle, yticklen=yticklen, $
;        ytickv=ytickv, yticks=yticks, yminor=yminor, $
;        noerase=1, position=tpos, nodata=1, charsize=1.2
;
;    plots, fs, avg_phase
;    foreach val, [-90,0,90] do oplot, xrange, [0,0]+val, linestyle=1
    
    
    plot, yrange, xrange, $
        ystyle=1, ylog=1, yrange=xrange, ytitle=xtitle, ytickname=xtickn, $
        yticklen=xticklen, ytickv=xtickv, yticks=xticks, yminor=xminor, $
        xstyle=1, xlog=0, xrange=yrange, xtitle=ytitle, xticklen=yticklen, $
        xtickv=ytickv, xticks=yticks, xminor=yminor, $
        noerase=1, position=tpos, nodata=1, charsize=1.2

    plots, avg_phase, fs
    foreach val, [-90,0,90] do oplot, [0,0]+val, xrange, linestyle=1
    
    if keyword_set(test) then stop
    sgclose

end
