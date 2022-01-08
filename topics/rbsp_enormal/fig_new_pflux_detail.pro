
; plot vsc and bmag.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
probes = ['a','b']
spinrate = 11
perp = '!9'+string(94b)+'!X'
labfac = ['||',perp+',West',perp+',North']
rgb = [6,4,2]
dr0 = 1d/16



device, decomposed = 0
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 4
tplot_options, 'constant', 0
tplot_options, 'labflag', -1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1

strmu = '!9'+string(109b)+'!X'
pre0 = 'rbspa_'




;----pflux spectrogram, normalized version.
tvar = pre0+'pf_fac'
get_data, tvar, uts
nrec = n_elements(uts)
dr0 = sdatarate(uts)

; change unit to (uW/m^2).
tvar = pre0+'pf_para_mor_spec'
get_data, tvar, uts, dat, val, limits = lim

dat *= 1e3
zr = [-1,1]*5
pfunit = strmu+'W/m!U2!N'
tvar = pre0+'pf_para_mor_spec_tmp'
store_data, tvar, uts, dat, val, limits = lim
options, tvar, 'zrange', zr
options, tvar, 'ztitle', 'S!D'+labfac[0]+'!N ('+pfunit+')'
options, tvar, 'yrange', [0.25,1000]


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
        idx = where(uts ge utr[0] and uts le utr[1], tnrec)
        psd = psd[idx,*]
        gws = total(psd,1)/tnrec^2
        ngws = (gws/ss)*(dr0*dj/cdelta)*tnrec
        store_data, tvar+'_tmp', ps, [[gws],[ngws]]
    endforeach


    ofn = 0
    ofn = shomedir()+'/fig_psbl_pflux_detail.pdf'

    sgopen, ofn, xsize=5, ysize=6, /inch

    posu = [0.25,0.35,0.85,0.95]
    posl = [0.15,0.10,0.95,0.25]
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
;--plot pflux 3d.
    poss = sgcalcpos(5,position=posu)

    device, decomposed=0
    loadct2, 43
    
    tvar = pre0+'pf_fac'
    options, tvar, 'ytickv', [0,0.1,0.2]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'yrange', [-0.05,0.2]
    options, tvar, 'ystyle', 1
    options, tvar, 'labels', 'S!D'+labfac
    
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
    ;options, tvar, 'labels', ['V!DSC!N','']
    
    vars = pre0+['n_combine','de_fac','db_fac','pf_fac']
    figlabs = ['a. Density', 'b. dE FAC', 'c. dB FAC', 'd. S FAC!C    in-situ']
    tpos = poss[*,0:3]
    tplot, vars, position=tpos, /noerase, /nouttick, trange=utr
    for i=0,3 do xyouts, tpos[0,i]-xchsz*12, tpos[3,i]-ychsz*0.5, /normal, figlabs[i]

    loadct2, 66
    vars = pre0+['pf_para_mor_spec_tmp']
    tpos = poss[*,4]
    tplot, vars, position=tpos, /noerase, trange=utr
    xyouts, tpos[0]-xchsz*12, tpos[3]-ychsz*0.5, /normal, 'e. S!D||!N PSD'
    
    ; add an arrow.
    get_data, pre0+'pf_fac', uts, dat
    idx = where(uts ge utr[0] and uts le utr[1])
    uts = uts[idx]
    dat = dat[idx]
    tmp = max(dat[*,0], /absolute, idx)
    tut = uts[idx]
    tx = poss[0,0]+(poss[2,0]-poss[0,0])*(tut-utr[0])/(utr[1]-utr[0])
    ty = poss[3,0]
    arrow, tx, ty+ychsz*0.8, tx, ty, /solid, /normal
    
    
;--plot spectrums.
    poss = sgcalcpos(1,3, position=posl, xpad=5)

    vars = pre0+['de_fac','db_fac']
    get_data, vars[0]+'_tmp', ps, edat
    get_data, vars[1]+'_tmp', ps, bdat
    
    
    xrange = [0.1,10000]
    xtickv = 10^smkarthm(-1,4,1,'dx')
    xtickn = '10!U'+string(alog10(xtickv),format='(I0)') & xtickn[1:*:2]=' '
    xticks = n_elements(xtickv)-1
    xminor = 10
    xtitle = 'Period (sec)'
    
    ; E spectrogram.
    tpos = poss[*,0]
    yrange = 5*[1e-5,10]
    ytickv = 10^smkarthm(-4,1,1,'dx')
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)') & ytickn[1:*:2]=' '
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytitle = '(mV/m)!U2!N'
    plot, ps, edat[*,0], /noerase, position=tpos, $
        /xlog, xstyle=1, xrange=xrange, xtitle=xtitle, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, $
        /ylog, ystyle=1, yrange=yrange, ytitle=ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn
    xyouts, tpos[0]-xchsz*5, tpos[3]-ychsz*0.5, /normal, $
        'f. |E|!U2!N'
    ; add line fit.
    idx = where(ps ge ps[0], tns)
    tps = ps[idx]
    res = linfit(alog10(tps), alog10(edat[idx,0]))
    tys = 10^(alog10(tps)*res[1]+res[0])
    plots, ps[idx], tys, color=6, linestyle=0
    tmp = convert_coord(tps[tns/2], tys[tns/2], /data, /to_normal)
    xyouts, tmp[0]+xchsz*0.5, tmp[1]-ychsz*1, color=6, /normal, $
        'f!U -'+sgnum2str(res[1],ndec=2)
        
    
    ; B spectrogram.
    tpos = poss[*,1]
    yrange = 5*[1e-8,10]
    ytickv = 10^smkarthm(-6,1,1,'dx')
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)') & ytickn[1:*:2]=' '
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytitle = '(nT)!U2!N'
    plot, ps, bdat[*,0], /noerase, position=tpos, $
        /xlog, xstyle=1, xrange=xrange, xtitle=xtitle, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, $
        /ylog, ystyle=1, yrange=yrange, ytitle=ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn
    xyouts, tpos[0]-xchsz*5, tpos[3]-ychsz*0.5, /normal, $
        'g. |B|!U2!N'
    ; add line fit.
    idx = where(ps ge ps[0], tns)
    tps = ps[idx]
    res = linfit(alog10(tps), alog10(bdat[idx,0]))
    tys = 10^(alog10(tps)*res[1]+res[0])
    plots, ps[idx], tys, color=6, linestyle=0
    tmp = convert_coord(tps[tns/2], tys[tns/2], /data, /to_normal)
    xyouts, tmp[0]+xchsz*0.5, tmp[1]-ychsz*1, color=6, /normal, $
        'f!U -'+sgnum2str(res[1],ndec=2)
    
    ; E/B in km/s.
    tpos = poss[*,2]
    yrange = 5*[1e2,1e4]
    ytickv = 10^smkarthm(3,4,1,'dx')
    ytickn = '10!U'+string(alog10(ytickv),format='(I0)')
    yticks = n_elements(ytickv)-1
    yminor = 10
    ytitle = '(km/s)'
    plot, ps, sqrt(edat[*,0]/bdat[*,0])*1e3, /noerase, position=tpos, $
        /xlog, xstyle=1, xrange=xrange, xtitle=xtitle, $
        xtickv=xtickv, xticks=xticks, xminor=xminor, xtickname=xtickn, $
        /ylog, ystyle=1, yrange=yrange, ytitle=ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, ytickname=ytickn
    xyouts, tpos[0]-xchsz*5, tpos[3]-ychsz*0.5, /normal, $
        'h. E/B'
        
    b0 = 240    ; nT.
    n0 = 0.5      ; cc.
    oratio = 0.67    ; 66.7% O.
    va_const = 21.8 ; nT/cc^0.5 -> km/s.
    vap = b0/sqrt(n0)*va_const
    vao = vap/sqrt(16*oratio+1*(1-oratio))
    
    plots, minmax(ps), vap+[0,0], color=6, linestyle=1
    tmp = convert_coord(max(ps), vap, /data, /to_normal)
    xyouts, tmp[0]-xchsz*5, tmp[1]+ychsz*0.5, color=6, /normal, $
        'v!DA!N all H!U+'
    
    plots, minmax(ps), vao+[0,0], color=2, linestyle=1
    tmp = convert_coord(min(ps), vao, /data, /to_normal)
    xyouts, tmp[0]+xchsz*0, tmp[1]-ychsz*1.5, color=2, /normal, $
        'v!DA!N '+sgnum2str(round(oratio*100))+'% O!U+'
        
    
    sgclose

end
