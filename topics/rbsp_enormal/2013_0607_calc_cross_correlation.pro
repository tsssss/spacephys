; calc cross correlation between the potentials of opposite probes
; save the cross correlation data


; settings.
utr = time_double(['2013-06-07/04:47','2013-06-07/05:05'])
autr = time_double(['2013-06-07/04:52','2013-06-07/04:58:30'])
butr = time_double(['2013-06-07/04:54:30','2013-06-07/05:01'])
probes = ['a','b']

reload = 0
noplot = 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1
maxdv0 = 10     ; lower limit for bad spin axis potential difference.
badv6utr = time_double(['2013-06-07/04:57:44','2013-06-07/04:59:46'])


vb1fn = shomedir()+'/psbl_de_32hz/2013_0607_vb1.tplot'
datfn = shomedir()+'/psbl_de_32hz/cross_correlation.tplot'

vars = ['vb1','vb1_low_res']                ; variable for burst data.
vb1vars = ['rbspa_'+vars,'rbspb_'+vars]
vars = 'corr_v'+['12','34','56']
datvars = ['rbspa_'+vars,'rbspb_'+vars]




; **** load VB1.
load = 0
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

load = 0
foreach tvar, datvars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        get_data, pre0+'vb1', uts, vb1
        nrec = n_elements(uts)
        scutr = minmax(uts)
        
        
        ; **** truncate the data according to 
        secdt = 5d/20   ; 0.25 sec, 1/40 of spin, or 9 deg.
        nsec = floor((scutr[1]-scutr[0])/secdt)
        sects = scutr[0]+findgen(nsec)*secdt
        
        ; settings for cross correlation.
        dr0 = uts[1]-uts[0]
        dtrg = [-1,1]*0.015      ; max time shift.
        dnrg = round(dtrg/dr0)
        dns = findgen(dnrg[1]*2+1)+dnrg[0]
        dts = dns*dr0
        ndt = n_elements(dts)
        
        ncorr = nsec-1
        corts = dblarr(ncorr)
        corrs = fltarr(ncorr,ndt,neb1)
        
        
        for i = 0, ncorr-1 do begin
            tutr = sects[i:i+1]
            corts[i] = (tutr[0]+tutr[1])*0.5
            
            ; slice out the data in the current sector.
            get_data, pre0+'vb1', uts, vb1
            idx = where(uts ge tutr[0] and uts le tutr[1], nrec)
            vb1 = vb1[idx,*]
            uts = uts[idx]
            
            
            ; check if the sector contains bad spin axis data.
            badspin = 0
            if nvb1 ge 4 then begin
                ; rbspb v6 has a non-physical offset within badv6utr
                if tprobe eq 'b' then begin
                    if min(uts) ge min(badv6utr[0]) and max(uts) le max(badv6utr[1]) then begin
                        vb1[*,5] = -vb1[*,5]    ; seems to be the wrong sign as well.
                        vb1[*,5] = vb1[*,5]-vb1[0,5]+vb1[0,4]
                    endif
                endif
                tmp = abs(vb1[*,5]-vb1[*,4])
                if max(tmp) ge maxdv0 then badspin = 1
            endif
            
            
            ; remove background, save to rbspx_vb1_sec.
            for l = 0, nvb1-1 do begin
                tvb1 = vb1[*,l]
                vb1bg = tvb1[0]-(tvb1[0]-tvb1[nrec-1])/(uts[nrec-1]-uts[0])*(uts-uts[0])
                tvb1 = tvb1-vb1bg
                vb1[*,l] = tvb1
                vb1[*,l] = vb1[*,l]-smooth(vb1[*,l],nrec/10)
            endfor


            ; calc correlation for each component.
            for j = 0, neb1-1 do begin
                v1 = vb1[*,j*2]
                v2 = vb1[*,j*2+1]
                corr = dblarr(ndt)
                
                if j eq 2 and badspin then begin
                    corr[*] = !values.d_nan
                endif else begin
                    for l = 0, ndt-1 do corr[l] = c_correlate(v1,v2,dns[l])
                endelse
                corrs[i,*,j] = corr
            endfor
            
            
            ; make plots.
            if noplot eq 1 then continue
            
            tvar = pre0+'vb1'
            store_data, tvar+'_sec', uts, vb1, limits = $
                {labels: 'V'+string(findgen(nvb1)+1,format='(I0)'), $
                colors:findgen(nvb1)+1,ytitle:'VB1!C(V)'}
            
            poss = sgcalcpos(neb1+1, ypad = 4)
            ofn = shomedir()+'/corr/'+pre0+'cross_correlation_'+ $
                time_string(tutr[0],tformat='YYYY_MMDD_hhmm_ss.ff')+'.pdf'
            ; ofn = 0
            sgopen, ofn, xsize = 11, ysize = 8.5, /inch
            device, decompose = 0
            loadct2, 43
            
            xchsz = double(!d.x_ch_size)/!d.x_size
            ychsz = double(!d.y_ch_size)/!d.y_size
            
            tpos = poss[*,0]
            tplot, tvar+'_low_res', position = tpos, /noerase, trange = scutr, $
                /novtitle, $
                title = 'RBSP-'+strupcase(tprobe)+', '+$
                time_string(tutr[0],tformat='YYYY-MM-DD/hh:mm:')+ $
                strjoin(time_string(tutr,tformat='ss.fff'),' to ')
            timebar, tutr, color = 6
            
            for k = 0, neb1-1 do begin
                corr = reform(corrs[i,*,k])
                maxcorr = max(corr, idx)
                drec = dns[idx]
                tmaxcorr = dts[idx]*1e3
                
                tpos = poss[*,k+1]
                tpos[2] = tpos[0]+(tpos[2]-tpos[0])*0.6
                v1 = vb1[*,k*2]
                v2 = vb1[*,k*2+1]
                tdat = [[v1],[shift(v2,-drec)]]
                tdat[0:abs(drec),*] = !values.d_nan
                tdat[nrec-1-abs(drec):nrec-1,*] = !values.d_nan
                tstr = ['V'+sgnum2str(k*2+1),'V'+sgnum2str(k*2+2)]
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
        
        
        for j = 0, neb1-1 do $
            store_data, pre0+'corr_v'+sgnum2str(2*j+1)+sgnum2str(2*j+2), $
                corts, reform(corrs[*,*,j]), dts, limits = $
                {spec:1, no_interp:1, $
                ytitle:'Time shift!C(msec)', zrange:[-1,1], $
                ztitle:'Correlation V'+sgnum2str(2*j+1)+',-V'+sgnum2str(2*j+2)}
    endforeach
    tplot_save, datvars, filename = datfn
endif



stop




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
        
        nrec = n_elements(uts)
        maxcorrdts = dblarr(nrec)
        for i = 0, nrec-1 do begin
            tcorrs = reform(corrs[i,*])
            maxcorr = max(tcorrs,idx)
            mincorr = abs(min(tcorrs))
            if mincorr gt maxcorr then continue ; anti-correlate.
            if maxcorr lt 0.5 then continue     ; weak anti-correlation.
            flag = float(tcorrs ge 0.5)
            flag = flag-shift(flag,1)
            tmp = where(flag eq 1, cnt)
            if cnt gt 1 then continue           ; waves.
            maxcorrdts[i] = dts[idx]
        endfor
        idx = where(abs(maxcorrdts) le dr0, cnt)
        if cnt ne 0 then maxcorrdts[idx] = !values.d_nan
        store_data, pre0+'corrdt'+tvar, uts, maxcorrdts, limits = {psym:1}
    endforeach
    
    poss = sgcalcpos(5)
    
    tplot, pre0+['corr_v12','corr_v34','corr_v56','u_gsm','vb1_low'], $
        position = poss, trange = tutr
    
    sym = 8
    tmp = findgen(21)/20*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, color = 255, /fill
    
    tpos = poss[*,0]
    get_data, pre0+'corrdt_v12', uts, maxcorrdts
    plot, minmax(uts), minmax(dts), /nodata, /noerase, position = tpos, $
        xstyle = 5, ystyle = 5
    oplot, uts, maxcorrdts, color = 255, psym = sym
    
    tpos = poss[*,1]
    get_data, pre0+'corrdt_v34', uts, maxcorrdts
    plot, minmax(uts), minmax(dts), /nodata, /noerase, position = tpos, $
        xstyle = 5, ystyle = 5
    oplot, uts, maxcorrdts, color = 255, psym = sym
    
    tpos = poss[*,2]
    get_data, pre0+'corrdt_v56', uts, maxcorrdts
    plot, minmax(uts), minmax(dts), /nodata, /noerase, position = tpos, $
        xstyle = 5, ystyle = 5
    oplot, uts, maxcorrdts, color = 255, psym = sym
    
endforeach



; construct velocity from rbspx_corrdt_vxx
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    scdr = (tprobe eq 'a')? 1d/4096: 1d/1024

    
    vars = pre0+'corrdt_v'+['12','34','56']
    get_data, vars[0], uts
    nrec = n_elements(uts)
    vgsm = dblarr(nrec,3)
    vmag = dblarr(nrec)
    
    get_data, vars[0], uts, udts
    get_data, vars[1], uts, vdts
    
    stop

    idx = where(finite(udts))
    udts = interpol(udts[idx],uts[idx],uts)
    
    vsp = 0.1/sqrt(udts^2+vdts^2)
    idx = where(finite(vmag))
    vsp = interpol(vsp[idx],uts[idx],uts)
    
    get_data, vars[1], uts, wdts
    idx = where(wdts le scdr)
    wdts[idx] = !values.d_nan
    vsa = 0.012/wdts
    
    if tprobe eq 'a' then begin
        ut0 = time_double('2013-06-07/04:53')
        idx = where(uts ge ut0)
        print, '<v> spin plane', mean(vsp[idx],/nan)
        print, '<v> spin axis', mean(vsa[idx],/nan)
    endif else begin
        print, '<v> spin plane', mean(vsp,/nan)
        print, '<v> spin axis', mean(vsa,/nan)
    endelse
    stop
endforeach


end
