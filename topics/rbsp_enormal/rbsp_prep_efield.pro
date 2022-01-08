;+
; on daily bases, preprocess efield for certain probe.
; subtract background E field near earth (<3 Re)
; remove bad E field when Eu and Ev do not match and
; when they are asymetric around 0.
; 
; calc E GSM from V, save the quaternion of rotation matrix uvw2gsm, 
; from which wsc_gsm can be extracted.
;-


pro rbsp_prep_efield, ut0, probes = tprobe, pad_nspin=pad_nspin


    if n_elements(tprobe) eq 0 then message, 'no probe ...'
    if n_elements(pad_nspin) eq 0 then pad_nspin = 4    ; pad #spin for flag.


;---Constant.
    secofday = 86400d
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    we = 2*!dpi/secofday  ; rad/sec.


;---Settings.
    dr0 = 1d/16     ; data rate, 16 Hz.
    spinrate = 11   ; sec, 10.9378 sec.
    blen1s = [100d,100,12]  ; use the physical length.

    tbuff = 3600d   ; 1 hour of data buffer.
    if n_elements(tprobe) eq 0 then tprobe = 'a'
    pre0 = 'rbsp'+tprobe+'_'

    rootdir = shomedir()+'/rbsp_prep_efield'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    
    utr = ut0-(ut0 mod secofday)+[0,secofday]
    utr1 = utr+[-1,1]*tbuff
    ut2016 = time_double('2015-10-01')
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
    
    
;---load spice product.
; rbspx_[pos_gsm,mlt,lshell,dis]

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
    
    
;---load vsvy and calc esvy (instead of load esvy).
; rbspx_[vsc,vsvy,e_uvw], interpolated to uniform times.
    
    tvar = pre0+'vsc'    
    dat = sread_rbsp_efw_l2(utr1, probes=tprobe, type='vsvy')
    if size(dat,/type) ne 8 then begin
        message, 'no V data ...', /continue
        return
    endif
    
    tuts = sfmepoch(dat.epoch, 'unix')
    uts = smkarthm(min(tuts),max(tuts),dr0, 'dx')
    nrec = n_elements(uts)
    vsvy = sinterpol(dat.vsvy, tuts, uts)
    ; calibrate vsvy, tried rbsp_efw_get_cal_params, but really no effect
    ; since offset is 0, gain is not used.  
    ; end up using V12, may consider more complex ways.
    vsc = mean(vsvy[*,0:1], dimension = 2)
    if utr[0] gt ut2016 then vsc = mean(vsvy[*,2:3], dimension = 2)
    store_data, tvar, uts, vsc, limits = {ytitle:'(V)', labels:'Vsc'}
    store_data, pre0+'vsvy', uts, vsvy[*,0:3], limits = $
        {ytitle:'(V)', colors:[1,2,3,4], labels:'V'+['1','2','3','4']}
  
    
    ; calc E uvw.
    tvar = pre0+'e_uvw'    
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
    store_data, tvar, uts, [[eu],[ev],[ew]], limits = $
        {ytitle:'(mV/m)', labels:'E'+uvw, colors:rgb}

;---Calc model E field: vxB velocity, and corotation efield, use model.
; rbspx_[bmod_gsm,emod_gsm,emod_uvw].
    get_data, pre0+'pos_gsm', uts, rgsm
    nrec = n_elements(uts)
    ; calc velocity.
    vgsm = dblarr(nrec,3)
    for i=0, 2 do vgsm[*,i] = deriv(uts,rgsm[*,i])
    vgsm *= re
    store_data, pre0+'vel_gsm', uts, vgsm, limits={ytitle:'(km/s)', colors:rgb, labels:'GSM V'+xyz}
    ; get model B field.
    ets = stoepoch(uts,'unix')
    bmodgsm = dblarr(nrec,3)
    emfisis = sread_rbsp_emfisis(utr1, coord='gsm', probe=tprobe)
    if size(emfisis,/type) ne 8 then begin  ; get model field.
        par = 2d
        for i=0, nrec-1 do begin
            tet = ets[i]
            geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
            xp = rgsm[i,0] & yp = rgsm[i,1] & zp = rgsm[i,2]
            geopack_igrf_gsm, xp,yp,zp, bx0,by0,bz0
            geopack_t89, par, xp,yp,zp, bx1,by1,bz1
            bmodgsm[i,*] = [bx0,by0,bz0]+[bx1,by1,bz1]
        endfor
    endif else begin
        tuts = sfmepoch(emfisis.epoch,'unix',/tt2000)
        bmodgsm = sinterpol(emfisis.mag, tuts, uts)
    endelse
    store_data, pre0+'bmod_gsm', uts, bmodgsm, limits={ytitle:'(nT)', colors:rgb, labels:'GSM Mod B'+xyz}

    ; calc vxB.
    evxb = scross(vgsm,bmodgsm)*1e-3
    store_data, pre0+'evxb_gsm', uts, evxb, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM VxB E'+xyz}

    ; calc corotation E.
    rgei = sgse2gei(sgsm2gse(rgsm,ets),ets)
    vcorgei = [[-rgei[*,1]],[rgei[*,0]],[dblarr(nrec)]]*(re*we)
    vcorgsm = sgse2gsm(sgei2gse(vcorgei,ets),ets)
    ecor = -scross(vcorgsm,bmodgsm)*1e-3
    store_data, pre0+'ecor_gsm', uts, ecor, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM Cor E'+xyz}

    ; the total model E field.
    emodgsm = evxb+ecor
    ; project to spin plane.
    get_data, pre0+'quvw2gsm', tuts, quvw2gsm
    quvw2gsm = qslerp(quvw2gsm, tuts, uts)
    muvw2gsm = transpose(qtom(quvw2gsm))
    wgsm = transpose(reform(muvw2gsm[2,*,*]))
    emodmag = snorm(emodgsm)
    emodmag *= sin(sang(emodgsm,wgsm))
    ; the measured emag for spin-plane efield.
    get_data, pre0+'e_uvw', tuts, euvw
    emag = snorm(euvw) & emag = smooth(emag, spinrate/dr0, /nan, /edge_truncate)
    emag = sinterpol(emag, tuts, uts)
    ; correction using linear regression.
    txs = emodmag & tys = emag
    tmp = linfit(txs, tys, yfit=tfs)
    ybar = mean(tys)
    ssr = total((tys-tfs)^2)
    sst = total((tys-ybar)^2)
    rval = sqrt(1-ssr/sst)
    ec = (rval ge 0.9)? tmp[1]: dblarr(nrec)+1
    for i=0,2 do emodgsm[*,i] *= ec
    store_data, pre0+'emod_gsm', uts, emodgsm, limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM Mod E'+xyz}
    store_data, pre0+'e_mag', uts, [[emag],[emodmag]], limits={colors:[6,2],labels:['real','model']}
    ; match the datarate of model E to EFW.
    get_data, pre0+'e_uvw', uts, euvw
    get_data, pre0+'emod_gsm', tuts, emodgsm
    emodgsm = sinterpol(emodgsm, tuts, uts)
    get_data, pre0+'quvw2gsm', tuts, quvw2gsm
    quvw2gsm = qslerp(quvw2gsm, tuts, uts)
    muvw2gsm = transpose(qtom(quvw2gsm))
    dat = emodgsm
    emoduvw = [$    ; gsm to uvw.
        [dat[*,0]*muvw2gsm[0,0,*] + dat[*,1]*muvw2gsm[0,1,*] + dat[*,2]*muvw2gsm[0,2,*]],$
        [dat[*,0]*muvw2gsm[1,0,*] + dat[*,1]*muvw2gsm[1,1,*] + dat[*,2]*muvw2gsm[1,2,*]],$
        [dat[*,0]*muvw2gsm[2,0,*] + dat[*,1]*muvw2gsm[2,1,*] + dat[*,2]*muvw2gsm[2,2,*]]]
    deuvw = euvw-emoduvw & deuvw[*,2] = 0
    store_data, pre0+'emod_uvw', uts, emoduvw, limits={ytitle:'(mV/m)', colors:rgb, labels:'Mod E'+uvw}
    store_data, pre0+'de_uvw', uts, deuvw, limits={ytitle:'(mV/m)', colors:rgb, labels:'dE'+uvw}


    ; remove remaining bg field around perigee.
    ; make a decaying window.
    mindis0 = 3 ; Re.
    get_data, pre0+'de_uvw', uts, deuvw
    nrec = n_elements(uts)
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
    store_data, pre0+'coef_e0', uts, ec
    store_data, pre0+'e0', uts, emagbg, limits = {ytitle:'|E| BG!C(mV/m)'}
    ; apply amplitude correction.
    for i=0, 1 do deuvw[*,i] *= ec
    store_data, pre0+'de_uvw', uts, deuvw, limits = {ytitle:'(mV/m)', colors:rgb, labels:'dE'+uvw}



;---Make flags to tell bad E field.
; rbspx_bade_flag.

    ; check Eu and Ev.
    ; rbspx_deu, rbspx_deu_[mat,amp,flag_spin10,flag_spin3,flag_wake].
    get_data, pre0+'de_uvw', uts, dat
    nrec = n_elements(uts)
    store_data, pre0+'deu', uts, dat[*,0]
    store_data, pre0+'dev', uts, dat[*,1]
    vars = pre0+['deu','dev']
    sclstrs = string([1,3,10],format='(I0)')    ; 1,1/3,1/10 spin period.
    ; rbspx_dex0_[mat,amp,flag].
    foreach tvar, vars do stplot_analysis_spinplane_efield, tvar, spinrate=spinrate
    
    ;---compare Eu and Ev.
;    get_data, pre0+'deu0_flag', flaguts
;    flagnrec = n_elements(flaguts)
;    maxde0 = 2d ;mV/m.
;    deamps = fltarr(flagnrec,2)
;    henv = fltarr(flagnrec,2)
;    pres = pre0+['deu0','dev0']
;    foreach pre1, pres, j do begin
;        get_data, pre1+'_amp', uts, dat
;        tdat = dat[*,2] ; 1 spin band.
;        for i=0, flagnrec-2 do begin
;            tidx = where(uts ge flaguts[i] and uts le flaguts[i+1])
;            deamps[i,j] = median(tdat[tidx])
;        endfor
;        tdat = dat[*,0] ; 1/10 spin band.
;        for i=0, flagnrec-2 do begin
;            tidx = where(uts ge flaguts[i] and uts le flaguts[i+1])
;            henv[i,j] = median(tdat[tidx])
;        endfor
;    endforeach
;    store_data, pre0+'amp', flaguts, deamps, limits={colors:[6,2], labels:['u','v']}
;    ampu = deamps[*,0]
;    ampv = deamps[*,1]
;    damp = ampu-ampv
;    henv = total(henv,2)*0.5+1
;    flags = (abs(damp) ge maxde0*henv)
;    store_data, pre0+'flag_deuv', flaguts, flags, limits={yrange:[-0.5,1.5]}

    maxde0 = 2d ; mV/m.
    padt = 300d ; sec.

    get_data, pre0+'de_uvw', uts, deuvw
    eu = deuvw[*,0]
    ev = deuvw[*,1]
    ampu = sqrt(eu^2+shift(eu,spinrate*0.25/dr0)^2)
    ampv = sqrt(ev^2+shift(ev,spinrate*0.25/dr0)^2)
    ampu = smooth(ampu, spinrate/dr0*0.2, /nan, /edge_truncate) ; remove small spikes.
    ampv = smooth(ampv, spinrate/dr0*0.2, /nan, /edge_truncate)
    ampu = smooth(ampu, spinrate/dr0, /nan, /edge_truncate)     ; remove spin tone.
    ampv = smooth(ampv, spinrate/dr0, /nan, /edge_truncate)
    store_data, pre0+'amp', uts, [[ampu],[ampv]], limits={ytitle:'(mV/m)',labels:['u','v'],colors:[6,2]}

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
        tuts = [tuts,uts[(i1+i2)*0.5]]
        i1 = i2
    endwhile
    henv = interpol(tf1s, tuts, uts, /nan)
    henv = smooth(henv, spinrate/dr0*10, /nan, /edge_truncate)+1    ; make it always loose the criteria.  

    damp = ampu-ampv
    ;damp = smooth(damp, spinrate/dr0, /nan, /edge_truncate)
    flags = (abs(damp) ge maxde0*henv)
    store_data, pre0+'damp', uts, damp
    store_data, pre0+'flag_deuv', uts, flags, limits={yrange:[-0.5,1.5]}

    
    
    
    ; combine all flags.
    ; pre0_bade_flag.
    get_data, pre0+'deu0_flag', flaguts
    flagnrec = n_elements(flaguts)
    vars = pre0+[['deu0','dev0']+'_flag']
    flags = bytarr(flagnrec)
    foreach tvar, vars do begin
        get_data, tvar, tuts, tdat
        flags = flags or tdat[*,0]
    endforeach
    ; convert to high res.
    flags = interpol(flags, flaguts, uts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'bade_flag', uts, flags, limits={yrange:[-0.5,1.5]}
        

    ; check eclipse time.
    tmp = sread_rbsp_eclipse_time(utr, probes = tprobe)
    flags = interpol(double(tmp.flags),tmp.uts,uts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'flag_eclipse', uts, flags, limits = $
        {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
        labels:['1:eclipse']}
        
    ; check sdt time.
    tmp = sread_rbsp_efw_sdt_time(utr[0], probes = tprobe)
    flags = interpol(double(tmp.flags),tmp.uts,uts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'flag_sdt', uts, flags, limits = $
        {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
        labels:['1:sdt']}
    

    ; combine flags.
    vars = pre0+['bade_flag','flag_eclipse','flag_sdt','flag_deuv']
    flags = bytarr(nrec)
    foreach tvar, vars do begin
        get_data, tvar, uts, dat
        flags = flags or dat
    endforeach
    padrec = pad_nspin*spinrate/dr0
    tmp = flags[1:nrec-1]-flags[0:nrec-2]
    tmp = tmp[1:nrec-2]*tmp[1:nrec-3]
    idx = where(tmp ne 0, cnt)+1
    if cnt ne 0 then begin
        for i=0, cnt-1 do begin
            i0 = (idx[i]-padrec)>0
            i1 = (idx[i]+padrec)<nrec-1
            flags[i0:i1] = 1
        endfor
    endif
    store_data, pre0+'bade_flag', uts, flags, limits=$
        {labels:'1: bad E', yrange:[-0.2,1.2], ystyle:1, yticks:1, ytickv:[0,1], panel_size:0.4, ytitle:''}
    
    get_data, pre0+'de_uvw', uts, deuvw
    idx = where(flags eq 1, cnt)
    if cnt ne 0 then deuvw[idx,*] = !values.d_nan
    store_data, pre0+'de_uvw', uts, deuvw



    ; trim data to utr0.    
    vars = pre0+['de_uvw','e_uvw',']
    foreach tvar, vars do begin
        get_data, tvar, uts, dat
        idx = where(uts ge utr[0] and uts le utr[1])
        store_data, tvar, uts[idx], dat[idx,*]
    endforeach
    
    

        
    ;---rotate uvw into gsm.
    ; rbspx_de_gsm.
        get_data, pre0+'de_uvw', uts, dat
        eu = dat[*,0]
        ev = dat[*,1]
        ew = dat[*,2]
        
        get_data, pre0+'quvw2gsm', tuts, quvw2gsm
        quvw2gsm = qslerp(quvw2gsm, tuts, uts)
        muvw2gsm = transpose(qtom(quvw2gsm))
        
        ex = eu*muvw2gsm[0,0,*] + ev*muvw2gsm[1,0,*] + ew*muvw2gsm[2,0,*]
        ey = eu*muvw2gsm[0,1,*] + ev*muvw2gsm[1,1,*] + ew*muvw2gsm[2,1,*]
        ez = eu*muvw2gsm[0,2,*] + ev*muvw2gsm[1,2,*] + ew*muvw2gsm[2,2,*]
        
        tvar = pre0+'de_gsm'
        store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
            {ytitle:'(mV/m)', colors:[6,4,2], labels:'GSM E'+['x','y','z']}
        
        vars = pre0+['vsvy','vel_gsm',$
            'e_mag','coef_e0','e0','deu','dev',$
            'deu0*','dev0*','amp']
        store_data, vars, /delete
        vars = pre0+['flag_'+['deuv','eclipse','sdt'],$
            'bade_flag','quvw2gsm','vsc','coef_e0']
        store_data, vars, /delete
        
        vars = pre0+['evxb_gsm','ecor_gsm','emod_gsm','emod_uvw']

end

probes = ['a']
utr0 = time_double(['2012-09-25','2016-12-31'])
utr0 = time_double(['2013-01-12','2013-01-12']) ; steps in ec.
;utr0 = time_double(['2013-11-27','2013-11-28']) ; should reduce 2fsp.
;utr0 = time_double(['2013-05-01','2013-05-02']) & probes = 'b'
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

;utr0 = time_double(['2013-06-07','2013-06-07'])

    

;---Loop through each day.
    uts = smkarthm(utr0[0], utr0[1], 86400, 'dx')
    foreach tut, uts do foreach tprobe, probes do $
        rbsp_prep_efield, tut, probes=tprobe


end
