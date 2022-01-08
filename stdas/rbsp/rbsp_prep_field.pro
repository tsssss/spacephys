;+
; on daily bases, preprocess E/B fields for certain probe.
; 
; The data quantities are: (rbx_)
;   'ut'. Uniform times at 16 sample/sec.
;   'db_gsm'. The perturbation B field, which equals total B minus model B.
;   'de_gsm'. The perturbation E field after removing model E, using Ew=0.
;   'dedot0_gsm'. Simlar to 'de_gsm', but using the Ew from E.B=0.
;   'bade_flags'. 1 for bad E, including flags for Eu, Ev, Euv comparison, SDT, eclipse.
;   'bade_flag'. The overall flag to be applied to model and perturbation E.
;   'badb_flag'. The flag to be applied to model and perturbation B.
;   'bw_ratio'. The ratio between Bw/|B| in UVW.
;   'ut_mod'. The times for model E and B fields.
;   'bmod_gsm'. The model B field calculated from EMFISIS total B.
;   'evxb_gsm'. The E field due to Vsc x total B.
;   'ecor_gsm'. The E field due to corotation.
;
; One set of time is at 16 sample/sec. The related variables are:
;   rbsp_[de_gsm,dedot0_gsm,db_gsm]
; One set of time is at 1 min. The related variables are:
;   rbsp_[bmod_gsm,evxb_gsm,ecor_gsm,bade_flags,bade_flag,badb_flag,bw_ratio]
; 
;-


pro rbsp_prep_field, ut0, probes = tprobe, pad_nspin=pad_nspin, gen_plot=gen_plot

;---Check inputs.
    if n_elements(tprobe) eq 0 then message, 'no probe ...'
    if n_elements(pad_nspin) eq 0 then pad_nspin = 4    ; pad #spin for flag.


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
    tbuff = 3600d   ; 1 hour of data buffer.
    if n_elements(tprobe) eq 0 then tprobe = 'a'
    pre0 = 'rbsp'+tprobe+'_'

    rootdir = sdiskdir('Research')+'/sdata/rbsp'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    cdfdir = rootdir+'/rbsp'+tprobe+'/ebfield/'+time_string(ut0,tformat='YYYY')
    if file_test(cdfdir,/directory) eq 0 then file_mkdir, cdfdir
    figdir = sdiskdir('Research')+'/splot/rbsp_ebfield/rbsp'+tprobe+'/'+time_string(ut0,tformat='YYYY')
    if file_test(figdir,/directory) eq 0 then file_mkdir, figdir
    
    utr = ut0-(ut0 mod secofday)+[0,secofday]
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
    store_data, pre0+'badb_flag', posuts, flags
    dbgsm = sinterpol(dat.mag, tuts, ebuts)-sinterpol(bmodgsm, posuts, ebuts)
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
    ; correction using linear regression.
    txs = emodmag & tys = emag
    tmp = linfit(txs, tys, yfit=tfs)
    ybar = mean(tys)
    ssr = total((tys-tfs)^2)
    sst = total((tys-ybar)^2)
    rval = sqrt(1-ssr/sst)
    ec = (rval ge 0.9)? tmp[1]: 1
    for i=0,2 do emodgsm[*,i] *= ec
    store_data, pre0+'emod_gsm', posuts, emodgsm, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM Mod E'+xyz}
    ; match the datarate of model E to EFW.
    get_data, pre0+'e_uvw', ebuts, euvw
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
    mindis0 = 3 ; Re.
    get_data, pre0+'de_uvw', ebuts, deuvw
    nrec = n_elements(ebuts)
    get_data, pre0+'dis', tuts, dis
    dis = interpol(dis, tuts, ebuts)
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


;---Make flags to tell bad E field.
; rbspx_bade_flag.
    ; check Eu and Ev.
    ; rbspx_[deu0,dev0]_[mat,amp,flag].
    get_data, pre0+'de_uvw', ebuts, dat
    nrec = n_elements(ebuts)
    store_data, pre0+'deu', ebuts, dat[*,0]
    store_data, pre0+'dev', ebuts, dat[*,1]
    vars = pre0+['deu','dev']
    foreach tvar, vars do stplot_analysis_spinplane_efield, tvar, spinrate=spinrate
    store_data, vars, /delete
    
    ; compare Eu and Ev.
    ; rbspx_flag_deuv.
    maxde0 = 2d ; mV/m.
    get_data, pre0+'de_uvw', ebuts, deuvw
    eu = deuvw[*,0]
    ev = deuvw[*,1]
    ampu = sqrt(eu^2+shift(eu,spinrate*0.25/dr0)^2)
    ampv = sqrt(ev^2+shift(ev,spinrate*0.25/dr0)^2)
    ampu = smooth(ampu, spinrate/dr0*0.2, /nan, /edge_truncate) ; remove small spikes.
    ampv = smooth(ampv, spinrate/dr0*0.2, /nan, /edge_truncate)
    ampu = smooth(ampu, spinrate/dr0, /nan, /edge_truncate)     ; remove spin tone.
    ampv = smooth(ampv, spinrate/dr0, /nan, /edge_truncate)
    store_data, pre0+'amp', ebuts, [[ampu],[ampv]], limits={ytitle:'(mV/m)',labels:['u','v'],colors:[6,2]}

    ; check the amplitude of high-freq waves.    
    tf0s = abs(eu-smooth(eu, spinrate/dr0*0.2, /nan, /edge_truncate))
    drec = spinrate/dr0
    tf1s = [] & tuts = []
    i1 = 0
    while i1 lt nrec-1 do begin
        i2 = i1+drec<(nrec-1)
        if i2-i1 lt 1 then continue
        tdat = tf0s[i1:i2]
        tdat = tdat[sort(tdat)]
        tf1s = [tf1s,max(tdat[0:(i2-i1)*0.8])]
        tuts = [tuts,ebuts[(i1+i2)*0.5]]
        i1 = i2
    endwhile
    henv = interpol(tf1s, tuts, ebuts, /nan)
    henv = smooth(henv, spinrate/dr0*10, /nan, /edge_truncate)+1    ; make it always loose the criteria.  

    damp = ampu-ampv
    flags = (abs(damp) ge maxde0*henv)
    store_data, pre0+'flag_deuv', ebuts, flags, limits={yrange:[-0.5,1.5]}

    ; check eclipse time.
    ; rbspx_flag_eclipse.
    tmp = sread_rbsp_eclipse_time(utr, probes = tprobe)
    flags = interpol(double(tmp.flags),tmp.uts,ebuts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'flag_eclipse', ebuts, flags, limits = $
        {ytitle:'', yrange:[-0.2,1.2], yticks:1, yminor:0, ystyle:1, $
        labels:['1:eclipse']}
    ; check sdt time.
    ; rbspx_flag_sdt.
    tmp = sread_rbsp_efw_sdt_time(utr[0], probes = tprobe)
    flags = interpol(double(tmp.flags),tmp.uts,ebuts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'flag_sdt', ebuts, flags, limits = $
        {ytitle:'', yrange:[-0.2,1.2], yticks:1, yminor:0, ystyle:1, $
        labels:['1:sdt']}

    ; need to expand eclipse, sdt, and deuv flags.
    ; also interpolate them to posuts.
    vars = pre0+['flag_eclipse','flag_sdt','flag_deuv']
    nrec = n_elements(posuts)
    padrec = 5  ; equiv. to 5 min.
    foreach tvar, vars do begin
        get_data, tvar, tuts, dat
        flags = interpol(dat,tuts, posuts)
        idx = where(flags ne 0, cnt)
        if cnt ne 0 then flags[idx] = 1
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
        store_data, tvar, posuts, flags
    endforeach
    ; interpolate the other flags to posuts.
    ; remove the detailed flags.
    vars = pre0+['deu0_flag','dev0_flag']
    foreach tvar, vars do begin
        get_data, tvar, tuts, dat
        flags = interpol(dat[*,0],tuts, posuts)
        idx = where(flags ne 0, cnt)
        if cnt ne 0 then flags[idx] = 1
        store_data, tvar, posuts, flags
    endforeach
    
    ; combine all flags.
    ; pre0_bade_[flag,flags].
    nrec = n_elements(posuts)
    vars = pre0+['deu0_flag','dev0_flag','flag_deuv','flag_eclipse','flag_sdt']
    flaglabs = ['eu','ev','deuv','eclipse','sdt']
    nflag = n_elements(vars)
    flagcols = findgen(nflag)+1
    flags = fltarr(nrec,nflag)
    foreach tvar, vars, i do begin
        get_data, tvar, tuts, tdat
        flags[*,i] = tdat
        idx = where(flags[*,i] ne 0, cnt)
        if cnt ne 0 then flags[idx,i] = 1-0.05*(i+1)
    endforeach
    store_data, pre0+'bade_flags', posuts, flags, $
        limits={yrange:[-0.2,1.2],labels:flaglabs,colors:flagcols, $
        ystyle:1, yticks:1, ytickv:[0,1], panel_size:0.4, ytitle:''}
    ; combine flags.
    flags = total(flags,2) ne 0
;    tmp = flags[1:nrec-1]-flags[0:nrec-2]
;    tmp = tmp[1:nrec-2]*tmp[1:nrec-3]
;    idx = where(tmp ne 0, cnt)+1
;    if cnt ne 0 then begin
;        for i=0, cnt-1 do begin
;            i0 = (idx[i]-padrec)>0
;            i1 = (idx[i]+padrec)<(nrec-1)
;            flags[i0:i1] = 1
;        endfor
;    endif
    store_data, pre0+'bade_flag', posuts, flags, limits=$
        {labels:'1: bad E', yrange:[-0.2,1.2], ystyle:1, yticks:1, ytickv:[0,1], panel_size:0.4, ytitle:''}


;---Trim data to utr0.
    vars = pre0+['de_uvw','e_uvw','db_gsm','bade_flag','badb_flag', $
        'bade_flags','evxb_gsm','ecor_gsm','bmod_gsm']
    foreach tvar, vars do begin
        get_data, tvar, ebuts, dat
        idx = where(ebuts ge utr[0] and ebuts le utr[1])
        store_data, tvar, ebuts[idx], dat[idx,*]
    endforeach
    
        
;---rotate uvw into gsm, calc E.B=0.
; rbspx_[de_gsm,dedot0_gsm,bw_ratio].
    ; rotate from uvw to gsm.
    get_data, pre0+'de_uvw', ebuts, dat
    get_data, pre0+'quvw2gsm', tuts, quvw2gsm
    quvw2gsm = qslerp(quvw2gsm, tuts, ebuts)
    muvw2gsm = transpose(qtom(quvw2gsm))
    ex = dat[*,0]*muvw2gsm[0,0,*] + dat[*,1]*muvw2gsm[1,0,*] + dat[*,2]*muvw2gsm[2,0,*]
    ey = dat[*,0]*muvw2gsm[0,1,*] + dat[*,1]*muvw2gsm[1,1,*] + dat[*,2]*muvw2gsm[2,1,*]
    ez = dat[*,0]*muvw2gsm[0,2,*] + dat[*,1]*muvw2gsm[1,2,*] + dat[*,2]*muvw2gsm[2,2,*]
    store_data, pre0+'de_gsm', ebuts, [[ex],[ey],[ez]], limits = $
        {ytitle:'(mV/m)', colors:rgb, labels:'GSM E'+xyz}
    
    ; rotate B from gsm to uvw.
    get_data, pre0+'de_uvw', tuts, deuvw
    get_data, pre0+'bmod_gsm', tuts, dat
    dat = sinterpol(dat, tuts, ebuts)
    bmag = snorm(dat)
    bu = dat[*,0]*muvw2gsm[0,0,*] + dat[*,1]*muvw2gsm[0,1,*] + dat[*,2]*muvw2gsm[0,2,*]
    bv = dat[*,0]*muvw2gsm[1,0,*] + dat[*,1]*muvw2gsm[1,1,*] + dat[*,2]*muvw2gsm[1,2,*]
    bw = dat[*,0]*muvw2gsm[2,0,*] + dat[*,1]*muvw2gsm[2,1,*] + dat[*,2]*muvw2gsm[2,2,*]
    deuvw[*,2] = -(deuvw[*,0]*bu + deuvw[*,1]*bv)/bw
    store_data, pre0+'de_uvw', ebuts, deuvw
    dat = deuvw
    ex = dat[*,0]*muvw2gsm[0,0,*] + dat[*,1]*muvw2gsm[1,0,*] + dat[*,2]*muvw2gsm[2,0,*]
    ey = dat[*,0]*muvw2gsm[0,1,*] + dat[*,1]*muvw2gsm[1,1,*] + dat[*,2]*muvw2gsm[2,1,*]
    ez = dat[*,0]*muvw2gsm[0,2,*] + dat[*,1]*muvw2gsm[1,2,*] + dat[*,2]*muvw2gsm[2,2,*]
    store_data, pre0+'dedot0_gsm', ebuts, [[ex],[ey],[ez]], limits = $
        {ytitle:'(mV/m)', colors:rgb, labels:'GSM Edot0'+xyz}
    get_data, pre0+'bade_flag', posuts
    store_data, pre0+'bw_ratio', posuts, interpol(bw/bmag,ebuts,posuts)
    
    
;---Save data to cdf.
    if keyword_set(gen_cdf) then begin
        cdfbfn = 'rbsp'+tprobe+'_ebfield_'+time_string(ut0,tformat='YYYY_MMDD')+'_v01.cdf'
        cdffn = cdfdir+'/'+cdfbfn
        if file_test(cdffn) ne 0 then file_delete, cdffn

        ginfo = {$
            title: 'RBSP E/B fields from EFW and EMFISIS',$
            text: 'Generated by Sheng Tian at the University of Minnesota'}
        scdfwrite, cdffn, gattribute=ginfo

        ;---Variables on ebuts.
        utname = 'ut'
        get_data, pre0+'db_gsm', ebuts
        ainfo = {$
            fieldnam: 'UT time', $
            units: 'sec', $
            var_type: 'support_data'}
        scdfwrite, cdffn, utname, value=ebuts, attribute=ainfo, cdftype='CDF_DOUBLE'

        vname = 'db_gsm'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'dB GSM', $
            units: 'nT', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=transpose(dat), attribute=ainfo, dimensions=[3], dimvary=[1]

        vname = 'de_gsm'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'dE GSM', $
            units: 'mV/m', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=transpose(dat), attribute=ainfo, dimensions=[3], dimvary=[1]

        vname = 'dedot0_gsm'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'dE dot0 GSM', $
            units: 'mV/m', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=transpose(dat), attribute=ainfo, dimensions=[3], dimvary=[1]


        utname = 'ut_mod'
        get_data, pre0+'bmod_gsm', posuts
        ainfo = {$
            fieldnam: 'UT time', $
            units: 'sec', $
            var_type: 'support_data'}
        scdfwrite, cdffn, utname, value=posuts, attribute=ainfo, cdftype='CDF_DOUBLE'

        vname = 'bmod_gsm'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'B model GSM', $
            units: 'nT', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=transpose(dat), attribute=ainfo, dimensions=[3], dimvary=[1]

        vname = 'evxb_gsm'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'Evxb GSM', $
            units: 'mV/m', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=transpose(dat), attribute=ainfo, dimensions=[3], dimvary=[1]

        vname = 'ecor_gsm'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'Ecor GSM', $
            units: 'mV/m', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=transpose(dat), attribute=ainfo, dimensions=[3], dimvary=[1]

        vname = 'bade_flag'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'Total bad E flag', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=dat, attribute=ainfo

        vname = 'badb_flag'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'Total bad E flag', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=dat, attribute=ainfo

        vname = 'bw_ratio'
        get_data, pre0+vname, tmp, dat
        ainfo = {$
            fieldnam: 'Bw/|B|', $
            var_type: 'data', $
            depend_0: utname}
        scdfwrite, cdffn, vname, value=dat, attribute=ainfo
    endif


;---Generate survey plot.
    if keyword_set(gen_plot) then begin
        figfn = figdir+'/'+$
            'rbsp'+tprobe+'_ebfield_'+time_string(ut0,tformat='YYYY_MMDD')+'.pdf'

        maxe = 100d ; mV/m.
        mine = 5d  ; mV/m.
    
        pre0 = 'rbsp'+tprobe+'_'
    
        ;---Apply flags to dE UVW and GSM.
        get_data, pre0+'bade_flag', posuts, flags
        get_data, pre0+'e_uvw', ebuts
        flags = interpol(flags, posuts, ebuts)
        idx = where(flags ne 0, cnt)
        if cnt ne 0 then flags[idx] = 1
        idx = where(flags eq 1, cnt)
        if cnt ne 0 then begin
            vars = pre0+['de_'+['uvw','gsm'],'dedot0_gsm']
            foreach tvar, vars do begin
                get_data, tvar, tmp, dat
                dat[idx,*] = !values.d_nan
                store_data, tvar, tmp, dat
            endforeach
        endif
        
        get_data, pre0+'badb_flag', posuts, flags
        flags = interpol(flags, posuts, ebuts)
        idx = where(flags ne 0, cnt)
        if cnt ne 0 then flags[idx] = 1
        idx = where(flags eq 1, cnt)
        if cnt ne 0 then begin
            vars = pre0+['db_gsm']
            foreach tvar, vars do begin
                get_data, tvar, tmp, dat
                dat[idx,*] = !values.d_nan
                store_data, tvar, tmp, dat
            endforeach
        endif    
        ;---Apply flags to Ew.
;        bwratio0 = 0.25d
;        bmag0 = 5d
;        get_data, pre0+'bw_ratio', posuts, bwratio
;        get_data, pre0+'bmod_gsm', posuts, bmodgsm
;        bmag = sinterpol(snorm(bmodgsm), posuts, ebuts)
;        bwratio = sinterpol(bwratio, posuts, ebuts)
;        idx = where(bmag lt bmag0 or abs(bwratio) lt bwratio0, cnt)
;        if cnt ne 0 then begin
;            vars = pre0+['de_uvw']
;            foreach tvar, vars do begin
;                get_data, tvar, tmp, dat
;                dat[idx,2] = !values.d_nan
;                store_data, tvar, tmp, dat
;            endforeach
;        endif
    
        ;---Options.
        stplot_split, pre0+'de_gsm', newnames=pre0+'de'+xyz, colors=rgb, labels='GSM dE'+xyz
        get_data, pre0+'de_gsm', tmp, degsm
        tmp = mine>max(sg_autolim(degsm))<maxe
        tvar = pre0+['de'+xyz,'de_uvw','e_uvw']
        options, tvar, 'yrange', [-1,1]*tmp
    
        tvar = pre0+'bade_flags'
        options, tvar, 'panel_size', 0.6
        
        get_data, pre0+'db_gsm', ebuts, dbgsm
        get_data, pre0+'dis', posuts, dis
        dis = interpol(dis, posuts, ebuts)
        tmp = max(sg_autolim(dbgsm[where(dis ge mindis0),*]))
        tvar = pre0+'db_gsm'
        options, tvar, 'yrange', [-1,1]*tmp
    
    
        ;---Plot.
        figfn = 0
        sgopen, figfn, xsize=8.5, ysize=11, /inch
    
        device, decomposed=0
        loadct2, 43
    
        ; down-sample to make the plot smaller.
        vars = pre0+['e_uvw','bade_flags','de_uvw','db_gsm','de'+xyz]
        foreach tvar, vars do begin
            get_data, tvar, uts, dat, limits=lim
            if tvar eq pre0+'de_uvw' then dat[*,2] = 0
            store_data, tvar+'tmp', uts[0:*:16], dat[0:*:16,*], limits=lim
        endforeach
        tplot, vars+'tmp', trange=ut0+[0,secofday]

        sgclose
stop
        vars = pre0+['de'+xyz]
        store_data, vars, /delete
    endif
    

;---Cleanup.
    tvar = ['','_'+['mat','amp','flag']]
    vars = pre0+['deu0'+tvar,'dev0'+tvar]
    store_data, vars, /delete
    
    vars = pre0+['vsvy','vel_gsm',$
        'e_mag','coef_e0','e0','deu','dev',$
        'deu0*','dev0*','amp']
    store_data, vars, /delete
    
    vars = pre0+['flag_'+['deuv','eclipse','sdt'],$
        'bade_flag','quvw2gsm','vsc','coef_e0']
    store_data, vars, /delete
    
    vars = pre0+['evxb_gsm','ecor_gsm','emod_gsm','emod_uvw']
    store_data, vars, /delete
end

probes = ['a','b']
utr0 = time_double(['2012-09-25','2015-12-31'])
utr0 = time_double(['2015-11-23','2015-12-31'])
utr0 = time_double(['2012-11-23','2012-11-23'])


;utr0 = time_double(['2012-10-14','2012-10-14'])  ; large E around 12, yes it's wake.
;utr0 = time_double(['2013-03-22','2013-03-22'])  ; large E on the edge of eclipse.
;utr0 = time_double(['2013-03-15','2013-03-15'])  ; E on the edge of eclipse.
;utr0 = time_double(['2013-04-01','2013-04-01'])  ; E on the edge of eclipse.
;utr0 = time_double(['2013-03-17','2013-03-17'])  ; large E during storms.
;utr0 = time_double(['2013-01-20','2013-01-20'])  ; large E around 12, yes it's wake.
;utr0 = time_double(['2013-06-10','2013-06-10'])  ; large E around 05:15, flags are fine.
;utr0 = time_double(['2013-06-29','2013-06-29'])  ; large E around 07 and 04:30, yes looks good.
;utr0 = time_double(['2014-04-28','2014-04-28'])  ; large E after 12, trouble day. may need to check slope of the 
;utr0 = time_double(['2013-03-29','2013-03-29'])  ; large E around 07:30, yes looks fine.
;utr0 = time_double(['2013-03-30','2013-03-30'])  ; large E around 04:50, looks fine.
;utr0 = time_double(['2015-05-14','2015-05-14'])  ; large E around 02:20, real?
;
;utr0 = time_double(['2013-01-12','2013-01-12'])  ; steps in ec.
;utr0 = time_double(['2013-11-27','2013-11-28'])  ; should reduce 2fsp.
;utr0 = time_double(['2012-11-14','2012-11-15']) 
;utr0 = time_double(['2012-09-25','2012-09-26'])
;utr0 = time_double(['2013-11-27','2013-11-27']) ; trouble day.
;utr0 = time_double(['2015-03-16','2015-03-16']) ; trouble day, 180 deg phase diff.
;utr0 = time_double(['2013-03-17','2013-03-18']) ; trouble day.
;utr0 = time_double(['2012-10-09','2012-10-09']) ; trouble day.
;utr0 = time_double(['2012-09-25','2012-09-26'])
;utr0 = time_double(['2015-01-10','2015-01-10']) ; contain sdt.

;utr0 = time_double(['2013-06-07','2013-06-07']) ; good wave & bad data.
;probes = ['b']

    

;---Loop through each day.
    uts = smkarthm(utr0[0], utr0[1], 86400, 'dx')
    foreach tut, uts do foreach tprobe, probes do $
        rbsp_prep_field, tut, probes=tprobe, /gen_plot


end