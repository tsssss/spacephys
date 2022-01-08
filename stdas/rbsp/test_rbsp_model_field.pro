;+
; Generate plots to show the effect of model field removal (E and B).
; The calculation follows rbsp_prep_field.pro
;-

;---Inputs.
    utr0 = time_double(['2013-06-25','2013-06-25:10:00'])
    utr0 = time_double(['2012-09-25','2012-09-25:10:00'])
    tprobe = 'a'


;---Constant.
    secofday = 86400d
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re    ; earth's radii.
    we = 2*!dpi/secofday  ; rad/sec, earth's rotation.


;---Settings.
    dr0 = 1d/16     ; data rate, 16 Hz.
    spinrate = 11   ; sec, 10.9378 sec.
    blen1s = [100d,100,12]  ; use the physical length.

    mindis0 = 4d    ; re.
    padrec = 4      ; 4 min for flag.
    tbuff = 0
    if n_elements(tprobe) eq 0 then tprobe = 'a'
    pre0 = 'rbsp'+tprobe+'_'

    ut0 = utr0[0]-(utr0[0] mod secofday)   
    utr = ut0+[0,secofday]
    utr1 = utr+[-1,1]*tbuff
    date = time_string(utr[0],tformat='YYYY_MMDD')
    
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    uvw = ['u','v','w']


    device, decomposed = 0
    loadct2, 43

    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'xticklen', -0.03
    tplot_options, 'yticklen', -0.005
    tplot_options, 'constant', 0
    tplot_options, 'xmargin', [18,10]


;---Load spice product.
; rbspx_[pos_gsm,mlt,lshell,mlat,dis,quvw2gsm]
    spice = sread_rbsp_spice_product(utr1, probe=tprobe)
    tuts = spice.ut_pos
    idx = uniq(tuts)
    uts = tuts[idx]
    dat = (spice.pos_gsm)[idx,*]
    store_data, pre0+'pos_gsm', uts, dat, limits={ytitle:'(Re)', colors:rgb, labels:'GSM '+xyz}
    dat = (spice.mlt)[idx,*]
    store_data, pre0+'mlt', uts, dat, limits={ytitle:'(hr)',labels:'MLT'}
    dat = (spice.lshell)[idx,*]
    store_data, pre0+'lshell', uts, dat, limits={ytitle:'(Re)',labels:'L-shell'}
    dat = (spice.mlat)[idx,*]
    store_data, pre0+'mlat', uts, dat, limits = {ytitle:'MLat (deg)'}
    dat = (spice.dis)[idx,*]
    store_data, pre0+'dis', uts, dat, limits = {ytitle:'Dist (Re)'}
    tuts = spice.ut_cotran
    idx = uniq(tuts)
    uts = tuts[idx]
    dat = spice.q_uvw2gsm[idx,*]
    store_data, pre0+'quvw2gsm', uts, dat
    

    
    
;---Load vsvy and calc euvw (instead of load euvw).
; rbspx_[e_uvw], interpolated to uniform times.
    dat = sread_rbsp_efw_l2(utr1, probes=tprobe, type='vsvy')
    if size(dat,/type) ne 8 then begin
        message, 'no V data ...', /continue
        return
    endif
    
    tuts = sfmepoch(dat.epoch, 'unix')
    ebuts = smkarthm(min(tuts),max(tuts),dr0, 'dx')
    nrec = n_elements(ebuts)
    vsvy = sinterpol(dat.vsvy, tuts, ebuts)
    ; calibrate vsvy, tried rbsp_efw_get_cal_params, but really no effect
    ; since offset is 0, gain is not used.  
    ; the following is based on linear regression on v1,v2 and eu, v3,v4 and ev.
    eu = (vsvy[*,0]-vsvy[*,1])/blen1s[0]*1e3   ; V -> V/m -> mV/m.
    ev = (vsvy[*,2]-vsvy[*,3])/blen1s[1]*1e3
    ;ew = (vsvy[*,4]-vsvy[*,5])/blen1s[2]*1e3
    ew = dblarr(nrec)
    ; remove dc-offset.
    nspin = 1
    tnrec = nspin*spinrate/dr0
    eu = eu-smooth(eu, tnrec, /edge_truncate, /nan)
    ev = ev-smooth(ev, tnrec, /edge_truncate, /nan)
    sdespike, eu
    sdespike, ev
    store_data, pre0+'e_uvw', ebuts, [[eu],[ev],[ew]], limits = $
        {ytitle:'(mV/m)', labels:'E'+uvw, colors:rgb}

;---Calc model E/B fields. E vxB, E corotation, B model from EMFISIS.
; rbspx_[bmod_gsm,db_gsm,ecor_gsm,evxb_gsm,emod_gsm,de_uvw].
    get_data, pre0+'pos_gsm', posuts, rgsm
    nrec = n_elements(posuts)
    ; calc velocity.
    vgsm = dblarr(nrec,3)
    for i=0, 2 do vgsm[*,i] = deriv(posuts,rgsm[*,i])*re
    ; get model B field from EMFISIS.
    bmodgsm = dblarr(nrec,3)
    dat = sread_rbsp_emfisis(utr1, coord='gsm', probe=tprobe)
    if size(dat,/type) ne 8 then begin
        message, 'no B data ...', /continue
        return
    endif
    tuts = sfmepoch(dat.epoch,'unix',/tt2000)
    bmodgsm = sinterpol(dat.mag, tuts, posuts)
    ; get model B field from T89.
    bmodt89gsm = dblarr(nrec,3)
    posets = stoepoch(posuts,'unix')
    get_data, pre0+'pos_gsm', posuts, rgsm
    par = 2
    for i=0, nrec-1 do begin
        tet = posets[i]
        xp = rgsm[i,0] & yp = rgsm[i,1] & zp = rgsm[i,2]
        geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
        geopack_igrf_gsm, xp,yp,zp, bxp,byp,bzp
        geopack_t89, par, xp,yp,zp, tbx,tby,tbz
        bmodt89gsm[i,*] = [bxp,byp,bzp]+[tbx,tby,tbz]
    endfor
    delbmag = abs(snorm(bmodt89gsm)-snorm(bmodgsm))
    delbmag0 = 200d
    flags = delbmag gt delbmag0
    idx = where(flags eq 1, cnt)
    padrec = 5  ; equiv. to 5 min.
    tmp = flags[1:nrec-1]-flags[0:nrec-2]
    tmp = tmp[1:nrec-2]*tmp[1:nrec-3]
    idx = where(tmp ne 0, cnt)+1
    if cnt ne 0 then begin
        for i=0, cnt-1 do begin
            i0 = (idx[i]-padrec)>0
            i1 = (idx[i]+padrec)<(nrec-1)
            flags[i0:i1] = 1
        endfor
    endif
    idx = where(flags eq 0)
    bmodgsm = sinterpol(bmodgsm[idx,*], posuts[idx], posuts)
    dbgsm = sinterpol(dat.mag, tuts, ebuts)-sinterpol(bmodgsm, posuts, ebuts)
    store_data, pre0+'badb_flag', posuts, flags
    nrec = n_elements(posuts)
    for i=0, 2 do begin
        tmp = dbgsm[*,i]
        sdespike, tmp
        dbgsm[*,i] = tmp
    endfor
    store_data, pre0+'bmod_gsm', posuts, bmodgsm, limits={ytitle:'(nT)', colors:rgb, labels:'GSM Mod B'+xyz}
    store_data, pre0+'db_gsm', ebuts, dbgsm, limits={ytitle:'(nT)', colors:rgb, labels:'GSM dB'+xyz}
    ; calc E vxB.
    evxb = scross(vgsm,bmodgsm)*1e-3    ; no minus sign! v is the s/c motion.
    store_data, pre0+'evxb_gsm', posuts, evxb, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM VxB E'+xyz}
    ; calc corotation E.
    rgei = sgse2gei(sgsm2gse(rgsm,posets),posets)
    vcorgei = [[-rgei[*,1]],[rgei[*,0]],[dblarr(nrec)]]*(re*we)
    vcorgsm = sgse2gsm(sgei2gse(vcorgei,posets),posets)
    ecor = -scross(vcorgsm,bmodgsm)*1e-3    ; has minus sign! v is the plasma motion
    store_data, pre0+'ecor_gsm', posuts, ecor, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM Cor E'+xyz}



    ; the total model E field.
    emodgsm = evxb+ecor
    ; project to spin plane.
    get_data, pre0+'quvw2gsm', tuts, quvw2gsm
    quvw2gsm = qslerp(quvw2gsm, tuts, posuts)
    muvw2gsm = transpose(qtom(quvw2gsm))
    wgsm = transpose(reform(muvw2gsm[2,*,*]))
    emodmag = snorm(emodgsm)
    emodmag *= sin(sang(emodgsm,wgsm))
    ; the measured emag for spin-plane efield.
    get_data, pre0+'e_uvw', tuts, euvw
    emag = snorm(euvw) & emag = smooth(emag, spinrate/dr0, /nan, /edge_truncate)
    emag = sinterpol(emag, tuts, posuts)
    store_data, pre0+'emag', posuts, [[emag],[emodmag]], limits={ytitle:'(mV/m)',labels:['|E|','|E|mod'],colors:[6,2]}
    store_data, pre0+'demag', posuts, emag-emodmag, limits={ytitle:'(mV/m)', labels:'|E|-|E|mod'}
    
    ; correction using linear regression.
    tys = emodmag & txs = emag
    tmp = linfit(txs, tys, yfit=tfs)
    ybar = mean(tys)
    ssr = total((tys-tfs)^2)
    sst = total((tys-ybar)^2)
    rval = sqrt(1-ssr/sst)
    ec = (rval ge 0.9)? tmp[1]: 1   ; ec should be the shorting factor???
    print, ec
    
    get_data, pre0+'e_uvw', ebuts, euvw
    for i=0,2 do euvw[*,i] *= ec
    store_data, pre0+'emod_gsm', posuts, emodgsm, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM Mod E'+xyz}
    ; match the datarate of model E to EFW.
    get_data, pre0+'emod_gsm', tuts, emodgsm
    emodgsm = sinterpol(emodgsm, tuts, ebuts)
    get_data, pre0+'quvw2gsm', tuts, quvw2gsm
    quvw2gsm = qslerp(quvw2gsm, tuts, ebuts)
    muvw2gsm = transpose(qtom(quvw2gsm))
    dat = emodgsm
    emoduvw = [$    ; gsm to uvw.
        [dat[*,0]*muvw2gsm[0,0,*] + dat[*,1]*muvw2gsm[0,1,*] + dat[*,2]*muvw2gsm[0,2,*]],$
        [dat[*,0]*muvw2gsm[1,0,*] + dat[*,1]*muvw2gsm[1,1,*] + dat[*,2]*muvw2gsm[1,2,*]],$
        [dat[*,0]*muvw2gsm[2,0,*] + dat[*,1]*muvw2gsm[2,1,*] + dat[*,2]*muvw2gsm[2,2,*]]]
    deuvw = euvw-emoduvw & deuvw[*,2] = 0
    store_data, pre0+'de_uvw', ebuts, deuvw, limits={ytitle:'(mV/m)', colors:rgb, labels:'dE'+uvw}

    ; remove remaining bg field around perigee.
    ; make a decaying window.
    mindis0 = 3.0  ; Re. The center of the decaying window.
    deldis0 = 0.1  ; Re. The half-width of the decaying window.
    get_data, pre0+'dis', tuts, dis
    dis = interpol(dis, tuts, uts)
    ; norm to 50:1, 1-> 1 spin.
    get_data, pre0+'de_uvw', ebuts, deuvw
    get_data, pre0+'dis', posuts, dis
    dr1 = sdatarate(posuts)
    nrec = n_elements(posuts)
    trate = 1
    tuts = posuts
    repeat begin
        ttuts = smkarthm(min(posuts),max(posuts), dr1*trate, 'dx')
        dis = sinterpol(dis, tuts, ttuts)
        tuts = ttuts
        drec = atan((dis-mindis0)/deldis0)
        drec = round((drec-min(drec))/(max(drec)-min(drec))*10+1)
        trate /= double(n_elements(posuts))/total(drec)
        drate = abs(n_elements(posuts)-total(drec))
    endrep until (drate le spinrate/dr0)
    tnrec = n_elements(drec)
    idx0s = lonarr(tnrec+1)
    for i=0, tnrec-1 do idx0s[i+1] = idx0s[i]+drec[i]
    idx1s = (idx0s[0:tnrec-1]+idx0s[1:tnrec])*0.5
    get_data, pre0+'demag', posuts, emag
    emagbg = interpol(emag[idx1s],posuts[idx1s],posuts)
;    emagbg = dblarr(tnrec)
;    for i=0, tnrec-1 do begin
;        i0 = idx0s[i]>0
;        i1 = idx0s[i+1]<(nrec-1)
;        if i0 ge i1 then continue
;        emagbg[i] = min(emag[i0:i1])
;    endfor
;    idx1s = (idx0s[0:tnrec-1]+idx0s[1:tnrec])*0.5
    store_data, pre0+'demag2', posuts, emag-emagbg
    stop




    ; remove remaining bg field around perigee.
    ; make a decaying window.
    mindis0 = 3.0  ; Re. The center of the decaying window.
    deldis0 = 0.2  ; Re. The half-width of the decaying window.
    get_data, pre0+'dis', tuts, dis
    dis = interpol(dis, tuts, uts)
    ; norm to 50:1, 1-> 1 spin.
    get_data, pre0+'de_uvw', ebuts, deuvw
    nrec = n_elements(ebuts)
    get_data, pre0+'dis', posuts, dis
    dr1 = sdatarate(posuts)
    tuts = smkarthm(min(posuts),max(posuts), dr1*0.5, 'dx')
    dis = sinterpol(dis, posuts, tuts)
    drec = atan((dis-mindis0)/deldis0)
    drec = ((drec-min(drec))/(max(drec)-min(drec))*50+1)
    store_data, pre0+'tmp', posuts, drec
    print, n_elements(ebuts)/total(drec)
    stop

    ; convert drec to time.
    duts = (max(ebuts)-min(ebuts))/total(drec)*drec
    duts *= spinrate*5/min(duts)
    stop
    
    tnrec = n_elements(tuts)-1
    emagbg = dblarr(tnrec)
    for i=0, tnrec-1 do begin
        idx = where(ebuts ge tuts[i] and ebuts le tuts[i+1])
        emagbg[i] = min(emag[idx])
    endfor
    store_data, pre0+'demag2', tuts, emagbg

    ylim, pre0+'demag', [-1,1]*5    
    tplot, pre0+['tmp','demag','demag2']
    
stop
    ; find apogee and perigee.
    df = dis[1:nrec-1]-dis[0:nrec-2]
    idx = where(df le 0, cnt)
    nodes = [1,df[0:nrec-3]*df[1:nrec-2],1] ; negative for node.
    idx0 = where(nodes le 0, nnode)
    flags = bytarr(nnode)       ; 1 for apogee, 0 for perigee.
    for i=0, nnode-1 do begin
        tdat = dis[(idx0[i]-5)>0:(idx0[i]+5)<(nrec-1)]
        if dis[idx0[i]] eq max(tdat) then flags[i] = 1 else $
            if dis[idx0[i]] eq min(tdat) then flags[i] = 0 else flags[i] = -1
    endfor
    
    idx = where(flags eq 1, cnt)
    idx1 = idx0[idx]
stop
    
    
    
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
    
    emag = sqrt(total(deuvw[*,0:1]^2,2))
    emagbg = dblarr(nrec)
    drec = (dis-mindis0)/(max(dis)-mindis0)>0
    drec = (drec^0.4*100+1)*spinrate/dr0*0.5
    for i = 0, nrec-1 do begin
        ti = i
        i1 = (ti-drec[i])>0
        i2 = (ti+drec[i])<nrec-1
        emagbg[ti] = min(emag[i1:i2])
    endfor
    ec = (1-emagbg/emag)
    ; apply amplitude correction.
    for i=0, 1 do deuvw[*,i] *= ec
    store_data, pre0+'de_uvw', ebuts, deuvw, limits = {ytitle:'(mV/m)', colors:rgb, labels:'dE'+uvw}
    demag2 = smooth(snorm(deuvw[*,0:1]),spinrate/dr0)
    demag2 = sinterpol(demag2, ebuts, posuts)
    get_data, pre0+'demag', posuts, dat
    store_data, pre0+'demag', posuts, [[dat],[demag2]], limits={ytitle:'(mV/m)', colors:[6,2], labels:'d|E|'+['raw','smth']}
    

;---Generate plot.
    emax0 = 300d

    get_data, pre0+'emag', tmp, dat
    erng = sg_autolim(dat)
    emax = max(abs(erng))<emax0
    tvar = pre0+['e_uvw']
    ylim, pre0+'e_uvw', [-1,1]*emax
    ylim, pre0+'emag', [0,emax]
    
    get_data, pre0+'demag', tmp, dat
    drec = 60
    erng = sg_autolim(dat[drec:-drec,*])
    emax = max(abs(erng))<emax0
    ylim, pre0+'demag', [-1,1]*emax
    ylim, pre0+'de_uvw', [-1,1]*emax
    vars = pre0+['e_uvw','emag','demag','de_uvw']
    tplot, vars
    stop



end
