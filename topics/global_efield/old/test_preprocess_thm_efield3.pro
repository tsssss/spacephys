;+
; Use L2 EFS spinfit data at 3 sec cadance.
;-


; good day.
;utr0 = time_double(['2008-02-03','2008-02-04'])
;utr1 = time_double(['2008-02-03/04:50','2008-02-03/05:05']) ; this is example of good data.
;tprobe = 'd'

;; have perigee field.
;utr0 = time_double(['2008-01-16','2008-01-17'])
;tprobe = 'a'


; bad day: charging?
utr0 = time_double(['2015-03-17','2015-03-18'])
tprobe = 'd'

; bad day: 1.5sec tone.
utr0 = time_double(['2009-01-09','2009-01-10'])
tprobe = 'a'

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
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    uvw = ['u','v','w']


;---Load data.
    ; thx_[pos_gsm,mlt,dis].
    posvar = ['Epoch','XYZ_GSM','DNEUTS','RADIUS','SM_LCT_T']
    dat = sread_thm_orbit(utr0, probes=tprobe, vars=posvar)
    if size(dat,/type) ne 8 then begin
        message, 'no pos data ...', /continue
        return
    endif
    
    tuts = sfmepoch(dat.epoch,'unix')
    store_data, pre0+'pos_gsm', tuts, [[dat.xyz_gsm],[dat.dneuts]], limits=$
        {ytitle:'(Re)', colors:[6,4,2,0], labels:['GSM '+['X','Y','Z'],'D NeutS!N'], constant:[0,8,-8]}
    store_data, pre0+'mlt', tuts, dat.sm_lct_t, limits={labels:'MLT',ytitle:'(hr)',ystyle:1,yrange:[0,24],yticks:4,yminor:6,constant:[6,12,18]}
    store_data, pre0+'dis', tuts, dat.radius


    ; thx_[edsl1,emag].
    efsvar0s = pre0+['efs_dot0_dsl','efs_dot0_time']
    efsvar1s = ['efs_dsl','efs_utsec']
    efs = sread_thm_efi_l2(utr0, vars=efsvar0s, newname=efsvar1s, probes=tprobe)
    tvar = pre0+'e_dsl1'
    uts = efs.efs_utsec
    nrec = n_elements(uts)
    edsl = efs.efs_dsl & edsl[*,2] = !values.d_nan
    ;for i=0,1 do edsl[*,i] = edsl[*,i]-smooth(edsl[*,i],twin/spinrate, /edge_truncate)
    store_data, tvar, uts, edsl[*,0:1], limits={ytitle:'(mV/m)', colors:rgb[0:1], labels:'DSL E'+xyz[0:1]}
    emag = sqrt(edsl[*,0]^2+edsl[*,1]^2)
    store_data, pre0+'emag', uts, emag, limits={ytitle:'(mV/m)', labels:['|E|']}


    
    ; thx_[b_dsl].
    fgmvar0s = pre0+['fgs_dsl','fgs_time']
    fgmvar1s = ['fgs_dsl','fgs_utsec']
    fgm = sread_thm_fgm_l2(utr0, vars=fgmvar0s, newname=fgmvar1s, probes=tprobe)
    
    tvar = pre0+'b_dsl'
    tuts = fgm.fgs_utsec
    bdsl = fgm.fgs_dsl
    store_data, tvar, tuts, bdsl, limits={ytitle:'(nT)', colors:rgb, labels:'DSL B'+xyz}
    

    tplot, pre0+['b_dsl','e_dsl1'], trange=utr0
    stop

;---Remove co-rotation field.
    ; thx_[e0,coef_e0].
    mindis0 = 3 ; Re.
    get_data, pre0+'dis', tuts, dis
    dis = interpol(dis, tuts, uts)
    win0 = double(dis le mindis0)
    winwd = 9600/spinrate
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
    get_data, pre0+'e_dsl1', limits=lim
    edsl[*,0] *= ec
    edsl[*,1] *= ec
    store_data, pre0+'e_dsl2', uts, edsl[*,0:1], limits=lim

    vars = pre0+['b_dsl','e_dsl1','e_dsl2']
    tplot, vars, trange=utr0
    stop
    


;---Calc E dot0.
    ; thx_[b_dsl1,b0,db_dsl,e_dsl].
    bmin = 5d   ; nT.
    rmin = 0.25 ; percent.
    twin = 600d  ; sec. time window to remove offset.

    
    get_data, pre0+'b_dsl', tuts, bdsl
    bdsl = sinterpol(bdsl, tuts, uts)
    for i=0, 2 do bdsl[*,i] = smooth(bdsl[*,i], twin/spinrate, /edge_mirror)
    bmag = snorm(bdsl)
    idx = where(bmag gt bmin and abs(bdsl[*,2])/bmag gt rmin, cnt)
    if cnt ne 0 then edsl[idx,2] = -(edsl[idx,0]*bdsl[idx,0]+edsl[idx,1]*bdsl[idx,1])/bdsl[idx,2]
    store_data, pre0+'e_dsl', uts, edsl, limits={ytitle:'(mV/m)', colors:rgb, labels:'DSL E'+xyz}
    store_data, pre0+'b0_dsl', uts, bdsl, limits={ytitle:'(nT)', colors:rgb, labels:'DSL B'+xyz}
    
    
    vars = pre0+['b0_dsl','b_dsl','e_dsl','e_dsl1','e_dsl2','pos_gsm']
    tplot, vars, trange=utr0
    
    stop

    

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
    
end
