;+
; Use thm l1 data to get the flags, similar to what i did for rbsp
; L1 V are of low resolution 4 S/s.
; L1 E are at 8 S/s, but only have data for half day.
; So have to go back to use L2 EFS.
;-


; good day.
;utr0 = time_double(['2008-02-03','2008-02-04'])
;utr1 = time_double(['2008-02-03/04:50','2008-02-03/05:05']) ; this is example of good data.
;tprobe = 'd'

;; have perigee field.
;utr0 = time_double(['2008-01-16','2008-01-17'])
;tprobe = 'a'


; bad day.
utr0 = time_double(['2015-03-17','2015-03-18'])
tprobe = 'd'


;utr0 = time_double(['2008-01-21','2008-01-22'])
;tprobe = 'a'


;---Settings.
    pre0 = 'th'+tprobe+'_'

    spinrate = 3d


    thm_init
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'xticklen', -0.03
    tplot_options, 'yticklen', -0.005
    tplot_options, 'constant', 0
    
    device, decomposed=0
    loadct2, 43


;---Load data.
; thx_[euvw,emag,vsc].

    ; thx_eff.
    thm_load_efi, probe=tprobe, trange=utr0, datatype='eff', level='l1', /onthefly_edc_offset
    tvar = pre0+'eff'
    get_data, tvar, uts, edsl
    edsl[*,2] = 0
    store_data, tvar, uts, edsl
    options, tvar, 'colors', [6,4,2]

    

    ; thx_[pos_gsm,mlt,dis].
    posvar = ['Epoch','XYZ_GSM','DNEUTS','RADIUS','SM_LCT_T']
    dat = sread_thm_orbit(utr0, probes=tprobe, vars=posvar)
    if size(dat,/type) ne 8 then begin
        message, 'no pos data ...', /continue
        return
    endif
    
    tuts = sfmepoch(dat.epoch,'unix')
    store_data, pre0+'pos_gsm', tuts, [[dat.xyz_gsm],[dat.dneuts]], limits=$
        {ytitle:'(Re)', colors:[6,4,2,0], labels:['GSM '+['X','Y','Z'],'D NeutS!N'], constant:0}
    store_data, pre0+'mlt', tuts, dat.sm_lct_t, limits={labels:'MLT',ytitle:'(hr)',ystyle:1,yrange:[0,24],yticks:4,yminor:6,constant:[6,12,18]}
    store_data, pre0+'dis', tuts, dat.radius


    dat = sread_thm_eff_l1(utr0, probes=tprobe, newname=['euvw','utsec'])
    if size(dat,/type) ne 8 then begin
        message, 'no V data ...', /continue
        return
    endif
    
    
    uts = dat.utsec
    dr0 = 0.125
    nrec = n_elements(uts)
    
    euvw = dat.euvw
    
    thm_get_efi_cal_pars, uts, 'eff', tprobe, cal_pars=cp
    blen = cp.boom_length*cp.boom_shorting_factor
    gain = -1e3*cp.edc_gain
    eu = euvw[*,0] & eu = (eu-smooth(eu, 2*spinrate/dr0, /edge_truncate, /nan))/blen[0]*gain[0]
    ev = euvw[*,1] & ev = (ev-smooth(ev, 2*spinrate/dr0, /edge_truncate, /nan))/blen[1]*gain[1]
    ew = dblarr(nrec)
    store_data, pre0+'euvw', uts, [[eu],[ev]], limits=$
        {ytitle:'(mV/m)', colors:[6,4], labels:'E'+['u','v']}
    emag = sqrt(eu^2+ev^2)
    store_data, pre0+'emag', uts, emag, limits=$
        {ytitle:'(mV/m)', labels:['|E|']}


    ; thx_[vsc,vsvy].
    tvar = pre0+'vsc'
    dat = sread_thm_vaf_l1(utr0, probes=tprobe, newname=['vsvy','utsec'])
    if size(dat,/type) ne 8 then begin
        message, 'no V data ...', /continue
    endif else begin
        tuts = dat.utsec
        ; using V12, may consider more complex ways.
        vsc = mean(dat.vsvy[*,0:1], dimension=2)*1e-3   ; from mV to V.
        store_data, tvar, tuts, vsc, limits = {ytitle:'Vsc!C(V)', labels:'Vsc'}
    endelse


;    tplot, pre0+['euvw','emag','dis'], trange=utr0
;    stop
    

;---Coefficient to remove co-rotation E field near earth (<3 Re).
; thx_[e0,coef_e0].

    ; make a decaying window.
    mindis0 = 3 ; Re.
    get_data, pre0+'dis', tuts, dis
    dis = interpol(dis, tuts, uts)
    win0 = double(dis le mindis0)
    winwd = 9600/dr0
    winedge = exp(-findgen(winwd)/winwd*5)
    idx = where(win0-shift(win0,1) eq 1, cnt)   ; left edge.
    for i = 0, cnt-1 do begin
        j0 = idx[i]
        dj = (winwd-1)<(j0-0)
        win0[j0-dj:j0] = reverse(winedge[0:dj])
    endfor
    
    idx = where(win0-shift(win0,-1) eq 1, cnt)   ; right edge.
    for i = 0, cnt-1 do begin
        j0 = idx[i]
        dj = (winwd-1)<(nrec-1-j0)
        win0[j0:j0+dj] = (winedge[0:dj])
    endfor
    
    emagbg = dblarr(nrec)
    drec = (dis-mindis0)/(max(dis)-mindis0)>0
    drec = (drec^0.4*100+1)*spinrate/dr0*0.5
    for i = 0, nrec-1 do begin
        ti = i
        i1 = (ti-drec[i])>0
        i2 = (ti+drec[i])<(nrec-1)
        emagbg[ti] = min(emag[i1:i2])
    endfor
    ec = (1-emagbg/emag)
    store_data, pre0+'coef_e0', uts, ec
    store_data, pre0+'e0', uts, emagbg, limits=$
        {ytitle:'(mV/m)',labels:'|E| BG'}

    
    
    ; apply bg removal coefficient and update Euvw.
    get_data, pre0+'euvw', limits=lim
    eu = eu*ec
    ev = ev*ec
    store_data, pre0+'esvy', uts, [[eu],[ev]], limits=lim
    get_data, pre0+'eff', uts, edsl
    edsl[*,0] *= ec
    edsl[*,1] *= ec
    store_data, pre0+'eff', uts, edsl
    

    
    
    
    
;---Make flags to tell bad E field.
    
    ; get the upper and lower envelope for 30 sec cadence.
    dt0 = 60d   ; sec.
    drec = dt0/dr0
    perc = 0.95
    padt = 300d ; sec.

    tdat = eu
    fmin = dblarr(nrec)
    fmax = dblarr(nrec)
    for i=0, nrec-1, drec do begin
        ti = i
        i1 = (ti-drec)>0
        i2 = (ti+drec)<(nrec-1)
        f0s = tdat[i1:i2]
        tnrec = n_elements(f0s)
        f1s = abs(f0s)
        vmax = f1s[sort(f1s)] & vmax = vmax[tnrec*perc]
        f2s = f0s & f2s[where(f1s ge vmax)] = !values.d_nan
        fmin[i1:i2] = min(f2s,/nan)
        fmax[i1:i2] = max(f2s,/nan)
    endfor

    umin = fmin
    umax = fmax
    store_data, pre0+'eu_env', uts, [fmin+fmax], limits=$
        {yrange:[-1,1]*1}
    
    
    tdat = ev
    fmin = dblarr(nrec)
    fmax = dblarr(nrec)
    for i=0, nrec-1, drec do begin
        ti = i
        i1 = (ti-drec)>0
        i2 = (ti+drec)<(nrec-1)
        f0s = tdat[i1:i2]
        tnrec = n_elements(f0s)
        f1s = abs(f0s)
        vmax = f1s[sort(f1s)] & vmax = vmax[tnrec*perc]
        f2s = f0s & f2s[where(f1s ge vmax)] = !values.d_nan
        fmin[i1:i2] = min(f2s,/nan)
        fmax[i1:i2] = max(f2s,/nan)
    endfor
    
    vmin = fmin
    vmax = fmax
    store_data, pre0+'ev_env', uts, [fmin+fmax], limits=$
        {yrange:[-1,1]*1}
    
    dmax = smooth(umax-vmax,padt/dr0)
    dmin = smooth(umin-vmin,padt/dr0)
    store_data, pre0+'dmax', uts, dmax, limits={yrange:[-1,1]*5}
    store_data, pre0+'dmin', uts, dmin, limits={yrange:[-1,1]*5}
    ;store_data, pre0+'dmax', uts, umax-vmax, limits={yrange:[-1,1]*5}
    ;store_data, pre0+'dmin', uts, umin-vmin, limits={yrange:[-1,1]*5}

    flags = (abs(dmax) ge 5) or (abs(dmin) ge 5)
    idx = where(flags eq 1, cnt)
    drec = padt/dr0
    if cnt ne 0 then begin
        for i=0, cnt-1 do flags[(idx[i]-drec)>0:(idx[i]+drec)<(nrec-1)] = 1
        idx = where(flags eq 1)
        get_data, pre0+'eff', uts, edsl
        edsl[idx,*] = !values.d_nan
        store_data, pre0+'eff', uts, edsl
    endif
    store_data, pre0+'bade_flag', uts, flags, limits=$
        {labels:'1: bad E', yrange:[-0.5,1.5], ystyle:1, yticks:1, ytickv:[0,1], panel_size:0.2, ytitle:''}



;---E dot B = 0.
    thm_load_fgm, probe=tprobe, trange=utr0, coord='dsl', datatype='fgl'
    get_data, pre0+'eff', uts, edsl, limits=lim
    get_data, pre0+'fgl', tuts, bdsl & bdsl=sinterpol(bdsl,tuts, uts)
    options, pre0+'fgl', 'colors', [6,4,2]
    edsl[*,2] = -(edsl[*,0]*bdsl[*,0]+edsl[*,1]*bdsl[*,1])/bdsl[*,2]
    bmag = snorm(bdsl)
    idx = where(bmag le 5 or bdsl[*,2]/bmag le 0.25, cnt)
    if cnt ne 0 then edsl[idx,2] = !values.d_nan
    store_data, pre0+'e0gsm', uts, edsl, limits=$
        {ytitle:'(mV/m)', labels:'GSM E'+['x','y','z'], colors:[6,4,2]}
    thm_cotrans, pre0+'e0gsm', in_coord='dsl', out_coord='gsm'
    options, pre0+'e0gsm', 'ytitle', '(mV/m)'
    
    tvar = pre0+'bgsm'
    thm_cotrans, pre0+'fgl', in_coord='dsl', out_coord='gsm'
    get_data, pre0+'fgl', tuts, bgsm & store_data, pre0+'fgl', /delete
    store_data, tvar, tuts, bgsm, limits=$
        {ytitle:'(nT)', colors:[6,4,2], labels:'GSM B'+['x','y','z']}
    vars = tnames(pre0+'state*')
    store_data, vars, /delete
    
    
    options, pre0+'pos_gsm', 'panel_size', 0.6
    options, pre0+'mlt', 'panel_size', 0.4
    tdegap, pre0+'bgsm', /overwrite, /nowarning
    tdegap, pre0+'e0gsm', /overwrite, /nowarning
    
    
    ofn = 0
    sgopen, ofn, xsize=11, ysize=8.5, /inch
    
    device, decomposed=0
    loadct2, 43
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tplot_options, 'ymargin', [5,5]
    tplot_options, 'xmargin', [15,15]
    
    figlabs = ['a','b','c','d','e','f']+'.'
    vars = pre0+['euvw','bade_flag','e0gsm','bgsm','pos_gsm','mlt']
    nvar = n_elements(vars)
    
    ;poss = sgcalcpos(nvar)
    titl = 'TH-'+strupcase(tprobe)+' '+time_string(utr0[0],tformat='YYYY-MM-DD')+' E field preprocess: remove co-rotation E, bad E, and E.B=0'
    tplot, vars, trange=utr0, get_plot_position=poss
    xyouts, (poss[0,0]+poss[2,0])*0.5, poss[3,0]+ychsz*0.8, /normal, alignment=0.5, titl, charsize=1.2
    for i=0, nvar-1 do xyouts, poss[0,i]-xchsz*8, poss[3,i]-ychsz*0.5, figlabs[i]

    sgclose

end
