;+
; rbspx_vb1. Raw VB1.
; rbspx_eb1. Raw EB1 calculated from VB1.
; rbspx_mb1. Raw MB1 loaded from cdf file.
; rbspx_[mb1,eb1]_uvw. Detrended EB1 and MB1, assume linear background.
;-

pro _2013_0607_0456_load_burst_data

    ; settings.
    utr0 = time_double('2013-06-07/04:52'+[':57.00',':57.90'])  ; the nice wave form, 50 Hz.
    ;utr0 = time_double('2013-06-07/04:55'+[':46.00',':48.00'])  ; p = 0.015 sec, 66 Hz.
    ;utr0 = time_double(['2013-06-07/04:55:59.80','2013-06-07/04:56:00.90'])
    utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; nice wave packet.
    ;utr0 = time_double(['2013-06-07/04:56:19.00','2013-06-07/04:56:20.40'])
    ;utr0 = time_double(['2013-06-07/04:56:25.00','2013-06-07/04:56:26.00'])

    ; chosen time, large S, during aurora, nice waveform.
    utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
    utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.

    ofn = shomedir()+'/psbl_de_32hz/2013_0607_0456_burst_data.tplot'
    tprobe='a'
    
    reload = 0
    
    
    pre0 = 'rbsp'+tprobe+'_'
    dr0 = 1d/4096
    tpad = 0.2  ; sec.
    utr = utr0+[-1,1]*tpad
    timespan, utr[0], utr[1]-utr[0], /second
    
    
    re = 6378d & re1 = 1d/re
    rgb = sgcolor(['red','green','blue'])
    c6s = sgcolor(['magenta','blue','cyan','green','yellow','red'])
    uvw = ['U','V','W']
    xyz = ['X','Y','Z']
    
    tplot_options, 'constant', 0
    tplot_options, 'labflag', -1
    tplot_options, 'num_lab_min', 5
    
    
    ; vars to be loaded and preprocessed.
    vars = ['mb1','vb1','eb1', $    ; raw data.
        'mb1_uvw','eb1_uvw', $
        ['u','v','w']+'_gsm','pos_gsm','b0_gsm','b0_uvw']
    var0s = pre0+vars
    
    
    ; vars on fac.
    fac = ['para','west','north']
    vars = [['mb1','eb1']+'_fac', $
        ['b','p','v']+'hat']
    var1s = pre0+vars
    
    ids = ['1','2','3']
    vars = ['eb1_fac_f'+ids,'mb1_fac_f'+ids]
    var2s = [pre0+vars,'morlet_info']
    allvars = [var0s,var1s,var2s]

    
    
    if file_test(ofn) eq 1 then tplot_restore, filename = ofn else reload = 1
    
    
    
    ; **** load data.
    load = 0
    foreach tvar, var0s do if tnames(tvar) eq '' then load = 1
    if reload eq 1 then load = 1
    if load then begin
        
        ; ** burst data.
        type = 'mscb1'
        tvar = pre0+'mb1'
        mb1 = sread_rbsp_efw_burst(utr, probes = tprobe, type = type)
        uts = sfmepoch(mb1.epoch,'unix')
        store_data, tvar, uts, mb1.mscb1, limits = $
            {ytitle:'(nT)', colors:rgb, labels:'dB '+uvw}
        rbsp_efw_cal_waveform, tvar, probe = tprobe, datatype = type, trange = utr
        
        type = 'vb1'
        tvar = pre0+'vb1'
        vb1 = sread_rbsp_efw_burst(utr, probes = tprobe, type = type)
        uts = sfmepoch(vb1.epoch,'unix')
        store_data, tvar, uts, vb1.vb1, limits = $
            {ytitle:'(V)', colors:c6s, labels:'V'+['1','2','3','4','5','6']}
        rbsp_efw_cal_waveform, tvar, probe = tprobe, datatype = type, trange = utr
        
        
        ; ** load UVW directions.
        defsysv,'!rbsp_spice', exists=flag
        if flag eq 0 then rbsp_load_spice_kernels, trange = utr
        
        uts = smkarthm(min(uts),max(uts),dr0*16,'dx')
        
        ; calc epoch times for spice (different than the epoch in cdfs).
        tmp = time_string(uts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        tets = tet0+uts-uts[0]
        
        scid = strupcase(pre0+'science')
        cspice_pxform, scid, 'GSM', tets, muvw2gsm
        
        ; muvw2gse[0,*,*] is u in GSM.
        ; muvw2gse[1,*,*] is v in GSM.
        ; muvw2gse[2,*,*] is w in GSM.
        store_data, pre0+'u_gsm', uts, transpose(reform(muvw2gsm[0,*,*])), $
            limits = {ytitle:'U', labels:'GSM '+xyz, colors:rgb}
        store_data, pre0+'v_gsm', uts, transpose(reform(muvw2gsm[1,*,*])), $
            limits = {ytitle:'V', labels:'GSM '+xyz, colors:rgb}
        store_data, pre0+'w_gsm', uts, transpose(reform(muvw2gsm[2,*,*])), $
            limits = {ytitle:'W', labels:'GSM '+xyz, colors:rgb}
        
        
        ; ** load B survey and sc position.
        rbsp_load_spice_state, probe = tprobe, coord = 'gsm', times = uts, /no_spice_load
        get_data, pre0+'state_pos_gsm', tmp, pos
        store_data, pre0+'state_*', /delete
        store_data, pre0+'pos_gsm', tmp, pos*re1, limits = $
            {yttile:'(Re)', colors:rgb, labels:'R GSM '+xyz}
        
        
        emfisis = sread_rbsp_emfisis_l3(utr, type = 'hires', probe = tprobe, coord = 'gsm')
        tuts = sfmepoch(emfisis.epoch,'unix')
        store_data, pre0+'b0_gsm', uts, sinterpol(emfisis.mag,tuts,uts), $
            limit = {ytitle:'(nT)', labels:'B0 GSM '+xyz, colors:rgb}
        
        ; convert B0 from gsm to uvw.
        get_data, pre0+'b0_gsm', uts, bgsm
        get_data, pre0+'u_gsm', tuts, ugsm & ugsm = sinterpol(ugsm,tuts,uts)
        get_data, pre0+'v_gsm', tuts, vgsm & vgsm = sinterpol(vgsm,tuts,uts)
        get_data, pre0+'w_gsm', tuts, wgsm & wgsm = sinterpol(wgsm,tuts,uts)
        buvw = [[sdot(bgsm,ugsm)],[sdot(bgsm,vgsm)],[sdot(bgsm,wgsm)]]
        store_data, pre0+'b0_uvw', uts, buvw, limits = $
            {ytitle:'B0 UVW!C(nT)', colors:rgb, labels:'B0 '+uvw}
    
    
        ; ** uniform time.
        ; uniform time for 'vb1','mb1','eb1', and all derived vars.
        ; rbspx_[mb1_uvw,eb1_uvw,vb1].
        tvar = pre0+'mb1'
        get_data, tvar, tuts, dat
        uts = smkarthm(min(tuts),max(tuts),dr0,'dx')
        nrec = n_elements(uts)
        dat = sinterpol(dat, tuts, uts)
        store_data, pre0+'mb1', uts, dat, limits = $
            {ytitle:'(nT)', colors:rgb, labels:'dB '+uvw}
        
        ; calc eb1 from vb1, 
        tvar = pre0+'vb1'
        get_data, tvar, tuts, dat
        dat = sinterpol(dat, tuts, uts)
        store_data, tvar, uts, dat
        dat = [(dat[*,0]-dat[*,1])/100d, (dat[*,2]-dat[*,3])/100d, $
            (dat[*,4]-dat[*,5])/12d]*1e3
        dat = reform(dat, [nrec,3])
        store_data, pre0+'eb1', uts, dat, limits = $
            {ytitle:'(mV/m)', colors:rgb, labels:'dE '+uvw}
            
            
            
        ; ** detrend, assume linear background.
        ; smoothing will reduce ~50% of the low frequency signal.
        ; trim edge > pad time.
        vars = pre0+['mb1','eb1']
        foreach tvar, vars do begin
            get_data, tvar, uts, dat, limits = lims
            for i = 0, 2 do begin
                tdat = dat[*,i]
                ; tdat = tdat-smooth(tdat, tpad/dr0, /edge_mirror)
                tdat = tdat-(tdat[0]+(tdat[nrec-1]-tdat[0])*findgen(nrec)/(nrec-1))
                dat[*,i] = tdat
            endfor
            store_data, tvar+'_uvw', uts, dat, limits = lims
        endforeach

        tplot_save, allvars, filename = ofn
    endif
    
    
    ; **** convert to FAC.
    load = 0
    foreach tvar, var1s do if tnames(tvar) eq '' then load = 1
    if reload eq 1 then load = 1
    if load then begin
    
        ; smooth B0 to remove wave.
        ; get direction for bhat, phat (west), vhat (north).
        b0p = 0.067     ; sec, the period in B0.
        b0w = 0.2       ; sec, width for smoothing.
        get_data, pre0+'vb1', uts
        get_data, pre0+'b0_uvw', tuts, dat
        dat = sinterpol(dat, tuts, uts)
        for i = 0, 2 do dat[*,i] = smooth(dat[*,i], b0w/dr0, /edge_mirror)
        bhat = sunitvec(dat)
        
        get_data, pre0+'vb1', uts
        get_data, pre0+'pos_gsm', tuts, dat
        dat = sinterpol(dat, tuts, uts)
        rhat = sunitvec(dat)
        phat = sunitvec(scross(rhat,bhat))
        vhat = scross(bhat,phat)
        
        store_data, pre0+'bhat', uts, bhat
        store_data, pre0+'phat', uts, phat
        store_data, pre0+'vhat', uts, vhat
        
        vars = pre0+['eb1','mb1']
        foreach tvar, vars do begin
            tmp = (tvar eq pre0+'eb1')? 'dE ': 'dB '
            get_data, tvar+'_uvw', uts, dat, limits = lim
            dat = [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]]
            store_data, tvar+'_fac', uts, dat, limits = $
                {ytitle:lim.ytitle, colors:rgb, labels:tmp+fac}
        endforeach
        
        tplot_save, allvars, filename = ofn
    endif
    
    
    
    ; **** filter use Morlet wavelet.
    load = 0
    foreach tvar, var2s do if tnames(tvar) eq '' then load = 1
    if reload eq 1 then load = 1
    if load then begin
        get_data, pre0+'vb1', uts
        nrec = n_elements(uts)
        
        
        w0 = 6
        s2p = 4*!dpi/(w0+sqrt(2+w0^2))
        
        ; coi.
        tmid = 0.5*(uts[0]+uts[nrec-1])
        coi = (tmid-uts[0]-abs(uts-tmid))/sqrt(2)
        
        ; prepare scales.
        s0 = 16*dr0
        dj = 0.125d
        nj = 40
        
        cdelta = 0.776d
        gamma = 2.32
        psi0 = !dpi^(-0.25)
        
        yfilts = 1d/[10,35,60,150]      ; in period, determined by looking at spectrogram.
        nfilt = n_elements(yfilts)-1
        
        vars = pre0+['mb1','eb1']+'_fac'
        foreach tvar, vars do begin     ; for each var.
            ytitle = (tvar eq pre0+'mb1_fac')? '(nT)': '(mV/m)'
            get_data, tvar, uts, dat
            
            fs = dblarr(nrec,3,nfilt)
            for i = 0, 2 do begin       ; for each component.
                tdat = dat[*,i]
                mor = wv_cwt(tdat, 'Morlet', w0, /pad, $
                    start_scale = s0/dr0, dscale = dj, nscale = nj, scale = recscls)
                timescls = recscls*dr0
                periods = timescls*s2p
                
                for j = 0, nfilt-1 do begin ; for each band.
                    filt = yfilts[j:j+1]
                    filt = filt[sort(filt)]
                    tys = dblarr(nrec)
                    idx = where(periods ge filt[0] and periods le filt[1], cnt)
                    for k = 0, cnt-1 do tys += real_part(mor[*,idx[k]]/sqrt(timescls[idx[k]]))
                    tys *= (dj*sqrt(dr0))/cdelta/psi0
                    fs[*,i,j] = tys
                endfor
            endfor
            
            for j = 0, nfilt-1 do begin
                id = string(j+1,format='(I0)')
                store_data, tvar+'_f'+id, uts, reform(fs[*,*,j]), limits = $
                    {ytitle:ytitle, colors:rgb, labels:fac}
            endfor
        endforeach
        
        info = {$
            s0:s0, $        ; min scale.
            dj:dj, $        ; scale step.
            nj:nj, $        ; # of scales.
            cdelta:cdelta, $;
            gamma:gamma, $  ;
            psi0:psi0, $    ; psi_0(0).
            filters:yfilts,$    ; filters in period in sec.
            scales:timescls,$   ; scales in sec.
            periods:periods}    ; periods in sec.
        store_data, 'morlet_info', 0, info

        tplot_save, allvars, filename = ofn
    endif
    
end
