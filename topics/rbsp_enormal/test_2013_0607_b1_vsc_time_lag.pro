

_2013_0607_load_burst_data

probes = ['a','b']


; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1


; **** make a low res version for plotting.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'vb1', uts, vb1, limits = lim
    store_data, pre0+'vb1_low', uts[0:*:1000], vb1[0:*:1000,*], limits = lim
    
    store_data, pre0+'vb1_12', uts, vb1[*,0:1], limits = $
        {ytitle:'(V)', labels:['V1','V2'], colors:[1,2]}
    store_data, pre0+'vb1_34', uts, vb1[*,2:3], limits = $
        {ytitle:'(V)', labels:['V3','V4'], colors:[3,4]}
    store_data, pre0+'vb1_56', uts, vb1[*,4:5], limits = $
        {ytitle:'(V)', labels:['V5','V6'], colors:[5,6]}
endforeach



fn = shomedir()+'/psbl_de_32hz/anti_corr.tplot'
tplot_restore, filename = fn
options, 'rbspa_corr_v'+['12','34','56'], 'ystyle', 1
options, 'rbspb_corr_v'+['12','34','56'], 'ystyle', 1

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    options, pre0+'corr_v*', 'ystyle', 1
    
    vars = ['_v12','_v34','_v56']
;    tutr = time_double(['2013-06-07/04:55:50','2013-06-07/04:56:30'])
;    tutr = time_double(['2013-06-07/04:53:15','2013-06-07/04:53:55'])
    foreach tvar, vars do begin
        get_data, pre0+'corr'+tvar, uts, corrs, dts
        dr0 = (tprobe eq 'a')? 1d/4096: 1d/1024
        tutr = minmax(uts)
 ;       idx = where(uts ge tutr[0] and uts le tutr[1])
 ;       uts = uts[idx]
 ;       corrs = corrs[idx,*]
        nrec = n_elements(uts)
        maxcorrdts = dblarr(nrec)
        for i = 0, nrec-1 do begin
            tcorrs = reform(corrs[i,*])
            maxcorr = max(tcorrs)
            mincorr = abs(min(tcorrs,idx))
;            if mincorr lt maxcorr then continue ; anti-correlate.
            if mincorr lt 0.5 then continue     ; weak anti-correlation.
            maxcorrdts[i] = dts[idx]
        endfor
        idx = where(abs(maxcorrdts) le dr0, cnt)
        if cnt ne 0 then maxcorrdts[idx] = !values.d_nan
        store_data, pre0+'corrdt'+tvar, uts, maxcorrdts
    endforeach
    
    poss = sgcalcpos(5)
    
    tplot, pre0+['corr_v12','corr_v34','corr_v56','u_gsm','vb1_low'], $
        position = poss, trange = tutr
    
    tpos = poss[*,0]
    get_data, pre0+'corrdt_v12', uts, maxcorrdts
    plot, minmax(uts), minmax(dts), /nodata, /noerase, position = tpos, $
        xstyle = 5, ystyle = 5
    oplot, uts, maxcorrdts, color = 255, psym = 1
    
    tpos = poss[*,1]
    get_data, pre0+'corrdt_v34', uts, maxcorrdts
    plot, minmax(uts), minmax(dts), /nodata, /noerase, position = tpos, $
        xstyle = 5, ystyle = 5
    oplot, uts, maxcorrdts, color = 255, psym = 1
    
    tpos = poss[*,2]
    get_data, pre0+'corrdt_v56', uts, maxcorrdts
    plot, minmax(uts), minmax(dts), /nodata, /noerase, position = tpos, $
        xstyle = 5, ystyle = 5
    oplot, uts, maxcorrdts, color = 255, psym = 1
    
    stop
endforeach

stop


; **** chop the data according to the spikes at half spin cadance.

maxdv0 = 10
padt0 = 0.05

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'vb1', uts, vb1
    nrec = n_elements(uts)
    scutr = minmax(uts)
    
    dvb1 = abs(vb1[*,4]-vb1[*,5])
    tmp = where(dvb1 gt maxdv0, cnt)
    if cnt eq 0 then message, 'no spin axis spikes ...'
    flags = [0,dvb1 gt maxdv0,0]
    tmp = flags[1:nrec+1]-flags[0:nrec]
    t1s = where(tmp eq 1)
    t2s = where(tmp eq -1)
    
    utrs = [[uts[t1s]],[uts[t2s]]]  ; trigger error when t1s/t2s mismatch.
    nutr = n_elements(utrs)/2
    for i = 0, nutr-1 do utrs[i,*]+= [-1,1]*padt0
    
    for i = 0, nutr-1 do begin
        idx = where(uts ge utrs[i,0] and uts le utrs[i,1])
        dvb1[idx] = !values.d_nan
    endfor
    
    ; half spin cadence.
    utrs = utrs[sort(utrs)]
    utrs = utrs[1:-2]
    nutr = n_elements(utrs)/2
    utrs = reform(utrs, 2, nutr)
    
    dr0 = sdatarate(uts)
    dtrg = [-1,1]*0.02      ; max time shift, 20 msec.
    dnrg = round(dtrg/dr0)
    dns = findgen(dnrg[1]*2+1)+dnrg[0]
    dts = dns*dr0
    ndt = n_elements(dts)
    
    nsec = 20           ; total length 5.5 sec, each sector is 9 deg.
    corrs = dblarr(nutr*nsec,ndt,3)
    corts = dblarr(nutr*nsec)

    for i = 0, nutr-1 do begin
        utr = utrs[*,i]     ; current time range.
        tutrs = smkarthm(utr[0],utr[1],nsec+1,'n')
        corts[i*nsec:i*nsec+nsec-1] = (tutrs[0:nsec-1]+tutrs[1:nsec])*0.5
        for j = 0, nsec-1 do begin
            tutr = tutrs[j:j+1]
            
            tvar = pre0+'vb1'
            get_data, tvar, uts, vb1
            
            idx = where(uts ge tutr[0] and uts le tutr[1], nrec)
            vb1 = vb1[idx,*]
            uts = uts[idx]
            if max(uts)-min(uts) lt dtrg[1] then continue ; section too short.
            
            ; remove background, which correlates and gives E.
            nvb1 = (size(vb1,/dimensions))[1]
            for l = 0, nvb1-1 do begin
                tvb1 = vb1[*,l]
                vb1bg = tvb1[0]-(tvb1[0]-tvb1[nrec-1])/(uts[nrec-1]-uts[0])*(uts-uts[0])
                tvb1 = tvb1-vb1bg
                vb1[*,l] = tvb1
                vb1[*,l] = vb1[*,l]-smooth(vb1[*,l],nrec/10)
                
                if l mod 2 eq 1 then vb1[*,l] *= -1    ; want anti-correlation.
            endfor
            
            store_data, tvar+'_sec', uts, vb1, limits = $
                {labels:'V'+string(findgen(nvb1)+1,format='(I0)'),colors:findgen(nvb1)+1,ytitle:'VB1!C(V)'}
                
                
            corr = dblarr(ndt)
            
            neb1 = nvb1/2
            poss = sgcalcpos(neb1+1, ypad = 4)
            ofn = shomedir()+'/corr2/'+pre0+'cross_correlation_'+time_string(tutr[0],tformat='YYYY_MMDD_hhmm_ss.ff')+'.pdf'
            ; ofn = 0
            sgopen, ofn, xsize = 11, ysize = 8.5, /inch
            device, decompose = 0
            loadct2, 43
            
            xchsz = double(!d.x_ch_size)/!d.x_size
            ychsz = double(!d.y_ch_size)/!d.y_size
            
            tpos = poss[*,0]
;            tplot, tvar+'_sec', position = tpos, /noerase, trange = tutr, $
            tplot, tvar+'_low', position = tpos, /noerase, trange = scutr, $
                /novtitle, $
                title = 'RBSP-'+strupcase(tprobe)+', '+$
                time_string(tutr[0],tformat='YYYY-MM-DD/hh:mm:')+ $
                strjoin(time_string(tutr,tformat='ss.fff'),' to ')
            timebar, tutr, color = 6
            
            for k = 0, neb1-1 do begin
                v1 = vb1[*,k*2]
                v2 = vb1[*,k*2+1]
                for l = 0, ndt-1 do corr[l] = c_correlate(v1,v2,dns[l])
                maxcorr = max(corr, idx)
                drec = dns[idx]
                tmaxcorr = dts[idx]*1e3
                corrs[i*nsec+j,*,k] = corr
                
                tpos = poss[*,k+1]
                tpos[2] = tpos[0]+(tpos[2]-tpos[0])*0.6
                tdat = [[v1],[shift(v2,-drec)]]
                tdat[0:abs(drec),*] = !values.d_nan
                tdat[nrec-1-abs(drec):nrec-1,*] = !values.d_nan
                tstr = ['V'+sgnum2str(k*2+1),'-V'+sgnum2str(k*2+2)]
                store_data, tvar+'_tmp', uts, tdat, limits = $
                    {labels:tstr,colors:k*2+[1,2],ytitle:strjoin(tstr,',')+'!Csmoothed,shfted!C(V)'}
                tplot, tvar+'_tmp', trange = tutr, position = tpos, /noerase, /novtitle
                
                
                xtitle = 'Time shift (msec)'
                tpos = poss[*,k+1]
                tpos[0] = tpos[0]+(tpos[2]-tpos[0])*0.75
                plot, dts*1e3, corr, /noerase, position = tpos, $
                    yrange = [-1,1], ystyle = 1, xstyle = 1, $
                    xtitle = xtitle, ytitle = 'Cross Correlation'
                plots, dtrg*1e3, [0,0], linestyle = 1
                plots, dtrg*1e3, 0.5+[0,0], linestyle = 1
                plots, dtrg*1e3,-0.5+[0,0], linestyle = 1
                
                plots, tmaxcorr, maxcorr, psym = 6, color = 6
                plots, tmaxcorr+[0,0], !y.crange, linestyle = 1
                xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1.2, /normal, $
                    'max='+string(maxcorr,format='(F4.1)')+', '+ $
                    sgnum2str(tmaxcorr,ndec=1)+' msec'
                
                
                mincorr = min(corr, idx)
                tmincorr = dts[idx]*1e3
                plots, tmincorr, mincorr, psym = 6, color = 6
                plots, tmincorr+[0,0], !y.crange, linestyle = 1
                xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*2.2, /normal, $
                    'min='+string(mincorr,format='(F4.1)')+', '+ $
                    sgnum2str(tmincorr,ndec=1)+' msec'
            endfor
            sgclose
        endfor
    endfor
    
    for k = 0, 2 do $
        store_data, pre0+'corr_v'+sgnum2str(2*k+1)+sgnum2str(2*k+2), $
        corts, reform(corrs[*,*,k]), dts, limits = $
        {spec:1, no_interp:1, $
        ytitle:'Time shift!C(msec)', zrange:[-1,1], $
        ztitle:'Correlation V'+sgnum2str(2*k+1)+',-V'+sgnum2str(2*k+2)}
endforeach



stop

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    
    ut1s = []
    ut2s = []
    
    
    get_data, pre0+'u_gsm', uts, dat
    ugsmy = dat[*,1]
    tmp = float([abs(ugsmy) ge 0.9,0])
    nrec = n_elements(tmp)
    idx = tmp[1:nrec-1]-tmp[0:nrec-2]
    
    ut1s = [ut1s,uts[where(idx eq  1)]]
    ut2s = [ut2s,uts[where(idx eq -1)]]
    
;    get_data, pre0+'v_gsm', uts, dat
;    vgsmy = dat[*,1]
;    tmp = float([abs(vgsmy) ge 0.7,0])
;    nrec = n_elements(tmp)
;    idx = tmp[1:nrec-1]-tmp[0:nrec-2]
;    
;    ut1s = [ut1s,uts[where(idx eq  1)]]
;    ut2s = [ut2s,uts[where(idx eq -1)]]
    
    
    nut = n_elements(ut1s)
    for i = 0, nut-1 do begin
        
        ; find where the maximum is.
        idx = where(uts ge ut1s[i] and uts le ut2s[i])
        tuts = uts[idx]
        tmp = max(abs(dat[idx]), idx)
        tutr = tuts[idx]+[-1,1]*0.05    ; 3 degree.
        
        if abs(tutr[0]-time_double('2013-06-07/04:54:54')) gt 11 then continue
        
        ; truncate to the time under focus, remove background.
        get_data, pre0+'vb1', uts, dat
        idx = where(uts ge tutr[0] and uts le tutr[1])
        tuts = uts[idx]
        tdat = dat[idx,*]
        tnrec = n_elements(tuts)
        for j = 0, 5 do begin
;            tdat[*,j] = tdat[*,j]-smooth(tdat[*,j],tnrec/8, /edge_truncate)
            tdat[*,j] = tdat[*,j]-tdat[0,j]
            tdat[*,j] = tdat[*,j]-findgen(tnrec)/(tnrec-1)*tdat[tnrec-1,j]
            tdat[*,j] = tdat[*,j]-smooth(tdat[*,j],tnrec/8, /edge_truncate)

            if j mod 2 eq 1 then tdat[*,j] *= -1    ; want anti-correlation.
        endfor
        
        dr0 = sdatarate(tuts)
        dtrg = [-1,1]*0.02   ; 20 msec.
        dnrg = ceil(dtrg/dr0)
        dns = smkarthm(dnrg[0],dnrg[1],1,'dx')  ; each shift in # of record.
        dts = dns*dr0
        ndt = n_elements(dts)
        corrs = fltarr(ndt)         ; the cross correlation at each shift.

        for k = 0, 1 do begin
            for j = 0, ndt-1 do corrs[j] = c_correlate(tdat[*,k*2], tdat[*,k*2+1], dns[j])
            maxcorr = max(corrs, idx)
            if maxcorr le 0.5 then continue ; no clear correlation.
            if idx eq 0 then continue       ; no time shift.
        endfor

stop

        ofn = shomedir()+'/test_vb1_vel/'+pre0+time_string(tutr[0],tformat='hhmm_ss')+'.pdf'
ofn = 0
        sgopen, ofn, xsize = 6, ysize = 4, /inch
        device, decomposed = 0
        loadct2, 43
       
        store_data, pre0+'vb1_12', tuts, tdat[*,0:1], limits = $
            {labels:['V1','-V2'], colors:[1,2]}
        store_data, pre0+'vb1_34', tuts, tdat[*,2:3], limits = $
            {labels:['V3','-V4'], colors:[1,2]}
        poss = sgcalcpos(3)
        tplot, pre0+['vb1_'+['12','34'],'u_gsm'], trange = tutr, position = poss, vlab_margin = 12
        sgclose
        stop
    endfor
endforeach



end
