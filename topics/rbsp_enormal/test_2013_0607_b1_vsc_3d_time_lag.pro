; calc cross correlation between the potentials of opposite probes
; save the cross correlation data


; settings.
utr = time_double(['2013-06-07/04:47','2013-06-07/05:05'])
autr = time_double(['2013-06-07/04:53:10','2013-06-07/04:53:55'])
probes = ['a','b']
boomlens = [100,100,12] ; m.


reload = 0
noplot = 1


tplot_options, 'labflag', -1


; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1
maxdv0 = 10     ; lower limit for bad spin axis potential difference.
badv6utr = time_double(['2013-06-07/04:57:44','2013-06-07/04:59:46'])

vlim = {colors:[1,2,3,4,5,6], labels:'V'+['1','2','3','4','5','6']}

vb1fn = shomedir()+'/psbl_de_32hz/2013_0607_vb1.tplot'
datfn = shomedir()+'/psbl_de_32hz/cross_correlation.tplot'

vars = ['vb1','vb1_low_res']                ; variable for burst data.
vb1vars = ['rbspa_'+vars,'rbspb_'+vars]
vars = 'corr_v'+['12','34','56']
datvars = ['rbspa_'+vars,'rbspb_'+vars]




; **** load VB1.
load = 0
if file_test(vb1fn) eq 1 then tplot_restore, filename = vb1fn
foreach tvar, vb1vars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        dr0 = (tprobe eq 'a')? 1d/4096: 1d/1024
        psbl_de_trim_b1_data, utr, tprobe
        
        tvar = pre0+'vb1'
        get_data, tvar, tuts, dat
        uts = smkarthm(utr[0], utr[1], dr0, 'dx')
        dat = sinterpol(dat, tuts, uts)
        store_data, tvar, uts, dat
        
        ; down sample to about 1 Hz.
        drec = 1/dr0
        uts = uts[0:*:drec]
        dat = dat[0:*:drec,*]
        store_data, tvar+'_low_res', uts, dat, limits = $
            {ytitle:'V 16Hz!C(V)', labels:'V'+['1','2','3','4','5','6'], colors:[1,2,3,4,5,6]}        
    endforeach
    
    options, 'rbspa_vb1'+['','_low_res'], 'yrange', [-35,-5]
    options, 'rbspb_vb1'+['','_low_res'], 'yrange', [-50,0]
    
    tplot_save, vb1vars, filename = vb1fn
endif






; **** calc cross correlation.

maxdv0 = 10
padt0 = 0.05
noplot = 0
nvb1 = 6
neb1 = 3

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    
    get_data, pre0+'vb1', uts, vb1, limits = lim
    nrec = n_elements(uts)
    scutr = (tprobe eq 'a')? autr: butr
    
    ; test on the half spin from 04:54:54.560 to 04:54:59.877
    ; truncate data to this time, remove high res spikes (<0.01 msec)
    spinutr0 = time_double(['2013-06-07/04:54:54.560','2013-06-07/04:54:59.877'])
    spinutr0 = time_double(['2013-06-07/04:53:15.000','2013-06-07/04:53:21.000'])
;    spinutr0 = time_double(['2013-06-07/04:53:15.000','2013-06-07/04:53:21.000'])+5.5
;    spinutr0 = time_double(['2013-06-07/04:49:38.550','2013-06-07/04:49:43.850'])   ; contain SDT, try despike.

    idx = where(uts ge spinutr0[0] and uts le spinutr0[1])
    uts = uts[idx]
    vb1 = vb1[idx,*]
    for i = 0, nvb1-1 do vb1[*,i] = vb1[*,i]-vb1[0,i]
    tvar = pre0+'vb1_current_spin'
    store_data, tvar, uts, vb1, limits = lim
    options, tvar, 'yrange', minmax(vb1)
    
    
    
    ; get the times for sectors.
    spindt0 = 5d ; use 5 sec of data out of ~5.5 sec.
    nsec = 40
    secdt = spindt0/nsec    ; 5/40 = 0.125 sec, 1/80 of spin, or 4.5 deg.
    secuts = mean(spinutr0)-0.5*spindt0+dindgen(nsec+1)*secdt
    
    
    
    ; get uvw2gsm.
    rbsp_load_spice_kernels, trange = spinutr0, probes = tprobe
    secmuts = (secuts[0:-2]+secuts[1:-1])*0.5
    tmp = time_string(secmuts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
    cspice_str2et, tmp, tet0
    tets = tet0+secmuts-secmuts[0]
    
    scid = strupcase(pre0+'science')
    cspice_pxform, scid, 'GSM', tets, muvw2gsm
    
    ; muvw2gse[3,3,nrec].
    ; muvw2gse[0,*,*] is u in GSM.
    ; muvw2gse[1,*,*] is v in GSM.
    ; muvw2gse[2,*,*] is w in GSM.
    store_data, pre0+'u_gsm', secmuts, transpose(reform(muvw2gsm[0,*,*])), $
        limits = {ytitle:'U GSM', labels:['x','y','z'], colors:[6,4,2]}
    store_data, pre0+'v_gsm', secmuts, transpose(reform(muvw2gsm[1,*,*])), $
        limits = {ytitle:'V GSM', labels:['x','y','z'], colors:[6,4,2]}
    store_data, pre0+'w_gsm', secmuts, transpose(reform(muvw2gsm[2,*,*])), $
        limits = {ytitle:'W GSM', labels:['x','y','z'], colors:[6,4,2]}

    
    
    
    ; settings for cross correlation.
    dr0 = uts[1]-uts[0]
    dtrg = [-1,1]*0.02      ; max time shift.
    dnrg = round(dtrg/dr0)
    dns = findgen(dnrg[1]*2+1)+dnrg[0]
    dts = dns*dr0
    ndt = n_elements(dts)
    
    
    ; remove high res signals, 0.01 sec, which corresponds to 100 m, 10 km/s.
    twidth = 0.01/dr0
    for i = 0, nvb1-1 do begin
        vb1[*,i] = smooth(vb1[*,i],twidth,/edge_truncate)
    endfor
    
    tvar = pre0+'vb1_current_spin_bg'
    store_data, tvar, uts, vb1, limits = lim
    options, tvar, 'yrange', minmax(vb1)

    
    
    
        
    corrs = dblarr(nsec,ndt,nvb1-1)
    corrdts = dblarr(nsec,nvb1-1)
    maxcorrs = dblarr(nsec,nvb1-1)
    maxcorr0 = 0.9
    
    vmags = dblarr(nsec)
    vangs = dblarr(nsec)
    yzangs = dblarr(nsec)
    spinangs = dblarr(nsec)
    uvwdts = dblarr(nsec,neb1)
    nhats = dblarr(nsec,3)
    

    for i = 0, nsec-1 do begin
        idx = where(uts ge secuts[i] and uts le secuts[i+1], tnrec)
        i0 = idx[0]
        i1 = idx[tnrec-1]
        
        tmp = vb1[i0:i1,0]
        v1 = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
        
        secvb1 = dblarr(tnrec,nvb1)
        secvb1[*,0] = v1
        
        for j = 1, nvb1-1 do begin
            for k = 0, ndt-1 do begin
                tmp = vb1[i0+dns[k]:i1+dns[k],j]
                v2 = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
                corrs[i,k,j-1] = c_correlate(v1,v2,0)
            endfor
        endfor
        
                
        for j = 0, nvb1-2 do begin
            maxcorrs[i,j] = max(corrs[i,*,j], idx)
            corrdts[i,j] = dts[idx]
            tmp = vb1[i0+dns[idx]:i1+dns[idx],j+1]
            secvb1[*,j+1] = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
        endfor
        store_data, pre0+'vb1_sec', uts[i0:i1], secvb1, limits = vlim
        
        
        
        for j = 0, nvb1-1 do begin
            tmp = vb1[i0:i1,j]
            secvb1[*,j] = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
        endfor
        store_data, pre0+'vb1_sec_orig', uts[i0:i1], secvb1, limits = vlim

        
        tmaxcorrs = reform(maxcorrs[i,*])
        tcorrdts = reform(corrdts[i,*])
        tcorrdns = tcorrdts/dr0
        idx = where(tmaxcorrs ge maxcorr0, cnt)
;        if cnt lt 4 then continue   ; does not contain enough good correlations.
        
        tplot, [10,11], trange = uts[[i0,i1]]
        print, tcorrdts/dr0, tmaxcorrs
        
        tmp = min(tmaxcorrs,idx)
        if idx[0] eq 0 then begin
            nmid = 0.8*(tcorrdns[1]+tcorrdns[2])+0.2*(tcorrdns[3]+tcorrdns[4])
            tcorrdns[0] = nmid
        endif
        if idx[0] eq 1 then begin
            nmid = 0.8*tcorrdns[0]+0.2*(tcorrdns[3]+tcorrdns[4])
            tcorrdns[1] = nmid-tcorrdns[2]
        endif
        if idx[0] eq 2 then begin
            nmid = 0.8*tcorrdns[0]+0.2*(tcorrdns[3]+tcorrdns[4])
            tcorrdns[2] = nmid-tcorrdns[1]
        endif
        if idx[0] eq 3 then begin
            nmid = 0.5*tcorrdns[0]+0.5*(tcorrdns[1]+tcorrdns[2])
            tcorrdns[3] = nmid-tcorrdns[4]
        endif
        if idx[0] eq 4 then begin
            nmid = 0.5*tcorrdns[0]+0.5*(tcorrdns[1]+tcorrdns[2])
            tcorrdns[4] = nmid-tcorrdns[3]
        endif
        
        
        ; try to use correlation along opposite probes.
        get_data, pre0+'vb1_sec', tmp, secvb1
        uvwmaxcorrs = dblarr(neb1)
        for j = 0, neb1-1 do begin
            uvwmaxcorrs[j] = c_correlate(secvb1[*,2*j],secvb1[*,2*j+1],0)
        endfor
        print, uvwmaxcorrs
        
        
        if uvwmaxcorrs[0] lt maxcorr0 then continue
        if uvwmaxcorrs[1] lt maxcorr0 then continue
        
        

        
        print, tcorrdns
        print, nmid
        tdns = [tcorrdns[0],tcorrdns[2]-tcorrdns[1],tcorrdns[4]-tcorrdns[3]]
        tdts = tdns*dr0
        uvwdts[i,*] = tdts
        
        corrdn1s = [tdns[0],[nmid+[[1,-1]*tdns[1],[1,-1]*tdns[2]]]*0.5]
        err = max(abs(corrdn1s-tcorrdns))
        print, corrdn1s
        
        nhat = tdns/boomlens
        nhat[2] = 0                     ; do not consider 56.
        nhat = sunitvec(nhat)           ; normal direction.
        spang = atan(nhat[1],nhat[0])*deg
        spang = (spang+360) mod 360
        
        vmag = abs(boomlens/tdts*nhat*1e-3)
        tmp = max(abs(nhat),idx)
        tmp = max(abs(nhat[0:1]),idx)
        vmag = vmag[idx]

        
        tut = mean(uts[[i0,i1]])
        spinang = (tut-uts[0])/11*360
        vuvw = vmag*nhat
        
        ; rotate from uvw to gsm.
        tmuvw2gsm = muvw2gsm[*,*,i]
        vgsm = vuvw # tmuvw2gsm
        yzang = atan(vgsm[2],vgsm[1])*deg
        
        vmags[i] = vmag
        spinangs[i] = spinang
        vangs[i] = spang
        yzangs[i] = yzang
        
        nhats[i,*] = nhat
        
        print, 'norm. dir in u-v plane:', spang, ' (deg)'
        print, 'norm. dir in y-z plane:', yzang, ' (deg)'
        print, 'spin angle:', spinang, ' (deg)'
        print, 'norm. dir in 3d:', nhat
        print, 'norm. vel in 3d:', vmag, ' (km/s)'
    endfor
    
    idx = where(vmags eq 0, cnt)
    vmags[idx] = !values.d_nan
    vangs[idx] = !values.d_nan
    spinangs[idx] = !values.d_nan
    yzangs[idx] = !values.d_nan
    nhats[idx,*] = !values.d_nan
    
;    vmags[where(vmags gt 200)] = !values.d_nan
    print, mean(vmags,/nan), stddev(vmags,/nan)
;    vangs[where(vangs ge 180)] -= 360
    
    
    plot, spinangs, vangs, psym = 1, /iso
    xr = !x.crange
    oplot, xr, -xr+mean(vangs+spinangs,/nan)
    stop
    
    
    
    
    
    ncorr = nsec-1
    corts = dblarr(ncorr)
    corrs = fltarr(ncorr,ndt,neb1)
    
    
    for i = 0, ncorr-1 do begin
        tutr = sects[i:i+1]
        
        ; this time 3-4 has no time lag.
        tutr = time_double('2013-06-07/04:53:26.125',tformat='YYYY-MM-DD/hh:mm:ss.fff')+[-1,1]*0.5*secdt
        tutr = time_double('2013-06-07/04:53:26.160',tformat='YYYY-MM-DD/hh:mm:ss.fff')+[-1,1]*0.5*secdt
        
        ; this time is within the sine wave period.
;        tutr = time_double('2013-06-07/04:52:11.625',tformat='YYYY-MM-DD/hh:mm:ss.fff')+[-1,1]*0.5*secdt
        
        
        corts[i] = (tutr[0]+tutr[1])*0.5
        
        ; slice out the data in the current sector.
        get_data, pre0+'vb1', uts, vb1
        idx = where(uts ge tutr[0] and uts le tutr[1], nrec)
        vb1 = vb1[idx,*]
        uts = uts[idx]
        
        ; original vb1.
        store_data, pre0+'vb1_sec0', uts, vb1, limits = $
            {labels: 'V'+string(findgen(nvb1)+1,format='(I0)'), $
            colors:findgen(nvb1)+1,ytitle:'VB1!C(V)'}
        
        
        ; remove background.
        for l = 0, nvb1-1 do begin
            tvb1 = vb1[*,l]
            vb1bg = tvb1[0]-(tvb1[0]-tvb1[nrec-1])/(uts[nrec-1]-uts[0])*(uts-uts[0])
            tvb1 = tvb1-vb1bg
            vb1[*,l] = tvb1
        endfor

        store_data, pre0+'vb1_sec1', uts, vb1, limits = $
            {labels: 'V'+string(findgen(nvb1)+1,format='(I0)'), $
            colors:findgen(nvb1)+1,ytitle:'VB1!C(V)'}
        
        
        ; use 3-4 to calc vsc as bg
;        vsc = (vb1[*,0]+vb1[*,1]+vb1[*,2]+vb1[*,3])*0.25
;        vb1[*,0] = vb1[*,0]-vsc
;        vb1[*,1] = vb1[*,1]-vsc
;        vb1[*,2] = vb1[*,2]-vsc
;        vb1[*,3] = vb1[*,3]-vsc
        
        store_data, pre0+'vb1_sec2', uts, vb1[*,0:3], limits = $
            {labels:'V'+['1','2','3','4'], $
            colors:[1,2,3,4],ytitle:'Vb1!C(V)'}
        
        
        tnrec = n_elements(uts)
        vb1bg = vb1
        for l = 0, nvb1-1 do begin
            vb1bg[*,l] = smooth(vb1[*,l],tnrec/40,/edge_truncate)
            vb1[*,l] = vb1[*,l]-vb1bg[*,l]
        endfor
        
        store_data, pre0+'vb1_sec3', uts, vb1, limits = $
            {labels: 'V'+string(findgen(nvb1)+1,format='(I0)'), $
            colors:findgen(nvb1)+1,ytitle:'VB1!C(V)'}
        store_data, pre0+'vb1_sec4', uts, vb1bg, limits = $
            {labels: 'V'+string(findgen(nvb1)+1,format='(I0)'), $
            colors:findgen(nvb1)+1,ytitle:'VB1!C(V)'}            
        
        tplot, pre0+'vb1_sec'+['1','3','4'], trange = tutr
        


        ; calc correlation for high freq part.
        tcorrs = dblarr(ndt,neb1)
        for j = 0, neb1-1 do begin
            v1 = vb1[*,j*2]
            v2 = vb1[*,j*2+1]
            corr = dblarr(ndt)

            for l = 0, ndt-1 do corr[l] = c_correlate(v1,v2,dns[l])
            tcorrs[*,j] = corr
        endfor
        


        
        ; calc correlation for each component.
        tbgcorrs = dblarr(ndt,nvb1-1)
        corrdts = dblarr(nvb1-1)
        maxcorrs = dblarr(nvb1-1)
        for j = 0, nvb1-2 do begin
            v1 = vb1bg[*,0]
            v2 = vb1bg[*,j+1]
            
            corr = dblarr(ndt)

            for l = 0, ndt-1 do corr[l] = c_correlate(v1,v2,dns[l])
            tbgcorrs[*,j] = corr
            maxcorrs[j] = max(corr, idx)
            corrdts[j] = dts[idx]
            plot, v1
            oplot, shift(v2,-dns[idx]), color = 6
            oplot, v2, color = 2
            
            print, dns[idx]
        endfor
                
        stop
    endfor
endforeach



stop




options, 'rbspa_corr_v'+['12','34','56'], 'ystyle', 1
options, 'rbspb_corr_v'+['12','34','56'], 'ystyle', 1



end
