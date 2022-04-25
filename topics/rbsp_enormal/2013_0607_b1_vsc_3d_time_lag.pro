; calc cross correlation between the potentials of opposite probes
; save the cross correlation data



; **** constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1
nan = !values.d_nan




; **** settings.
noplot = 1
reload = 0


utr = time_double(['2013-06-07/04:47','2013-06-07/05:05'])
figdir = shomedir()+'/cc3d'

probes = ['a','b']
boomlens = [100,100,12] ; m.
maxdv0 = 10     ; lower limit for bad spin axis potential difference.
nvb1 = 6
neb1 = 3

badv6utr = time_double(['2013-06-07/04:57:44','2013-06-07/04:59:46'])

; settings for dividing half spin period and sub-sectors.
spinrate0 = 5.531*2       ; sec.
halfspinsecdt0 = 5d ; sec.
nsec0 = 40          ; 5/40 = 0.125 sec, 1/80 of a spin, or 4.5 deg.
secdt0 = halfspinsecdt0/nsec0
secdt0str = sgnum2str(secdt0)

; settings for cross correlation.
ccmaxdt0 = 0.02     ; sec, or 20 msec.
mincc0 = 0.9    ; must be > this value to be considered as correlated.
ccinfo0 = { $
    cc1x:dblarr(nvb1)+nan, $    ; the max cc b/w 1 and x.
    dn1x:dblarr(nvb1)+nan, $    ; the record lag b/w 1 and x.
    ccuvw:dblarr(neb1)+nan, $   ; the max cc b/w [12],[34],[56].
    dnuvw:dblarr(neb1)+nan, $   ; the record lag b/w opposite probes.
    dnmid:nan, $                ; center time converted to @ of record.
    dnerr:nan, $                ; error in determining center time.
    nhat:dblarr(3), $           ; normal direction of the plane.
    vmag:nan}                   ; normal velocity in km/s.


vb1fn = shomedir()+'/psbl_de_32hz/2013_0607_vb1.tplot'
datfn = shomedir()+'/psbl_de_32hz/cross_correlation.tplot'



tplot_options, 'labflag', -1
tplot_options, 'xticklen', 1
tplot_options, 'num_lab_min', 10
tplot_options, 'version', 0
tplot_options, 'constant', 0
vlim = {ytitle:'VB1!C(V)',colors:[1,2,3,4,5,6], labels:'V'+['1','2','3','4','5','6']}


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



; **** load spice kernels
defsysv,'!rbsp_spice', exists=flag
if flag eq 0 then rbsp_load_spice_kernels, trange = utr





; **** calc cross correlation.
foreach tprobe, probes do begin
    
    pre0 = 'rbsp'+tprobe+'_'
    scid = strupcase(pre0+'science')    ; used in spice kernal.

    get_data, pre0+'vb1', uts, vb1
    nrec = n_elements(uts)
    dr0 = uts[1]-uts[0]
    
    spikewd0 = 0.01/dr0     ; corresponds to 0.01 sec, or 100 m, 10 km/s
    
    ; prepare to do cross correlation.
    ccdtrg = [-1,1]*ccmaxdt0
    ccdnrg = round(ccdtrg/dr0)
    ccdns = findgen(ccdnrg[1]*2+1)+ccdnrg[0]
    ccdts = ccdns*dr0
    nccdt = n_elements(ccdts)
    
    
; **** chop data according to spikes at half-spin cadence.
; then loop through each half-spin period.
    
    dvb1 = abs(vb1[*,4]-vb1[*,5])
    tmp = where(dvb1 gt maxdv0, cnt)
    if cnt eq 0 then message, 'no spin axis spikes ...'
    flags = [0,dvb1 gt maxdv0,0]
    tmp = flags[1:nrec+1]-flags[0:nrec]
    t1s = uts[where(tmp eq 1)]
    t2s = uts[where(tmp eq -1)]
    nrec = n_elements(t1s)
    
    
    ; times where spin-axis voltage is bad.
    ; spinuts is the start times of each half-spin period.
    ; the duration is 5 sec.
    spinuts = (t1s+t2s)*0.5
    spindts = spinuts[1:*]-spinuts[0:-2]
    idx = where(spindts ge 5 and spindts le 6)
    halfspinrate = mean(spindts[idx])
    nspinut = ceil((spinuts[-1]-spinuts[0])/halfspinrate)
    spinuts = spinuts[0]+findgen(nspinut)*halfspinrate


    foreach tspinut, spinuts do begin
 
; if tspinut lt time_double('2013-06-07/04:53:15') then continue

        ; time range of the current half spin period.
        ; shrink to 5 sec to leave enough padding time.
        spinutr0 = tspinut+[0,0.5*spinrate0]
        spinutr1 = tspinut+0.25*spinrate0+[-1,1]*halfspinsecdt0*0.5
        
        ; truncate to the current time.
        get_data, pre0+'vb1', uts, vb1
        idx = where(uts ge spinutr0[0] and uts le spinutr0[1])
        uts = uts[idx]
        vb1 = vb1[idx,*]
        for i = 0, nvb1-1 do vb1[*,i] = vb1[*,i]-vb1[0,i]
        
        if tprobe eq 'b' then begin
            if min(uts) ge badv6utr[0] and max(uts) le badv6utr[1] then begin
                if nvb1 ge 4 then begin
                    vb1[*,5] = -vb1[*,5]
                endif
            endif
        endif
        store_data, pre0+'vb1_current_spin0', uts, vb1, limits = vlim
        
        
        ; remove high freq signals, most likely to be temporal.
        for i = 0, nvb1-1 do begin
            vb1[*,i] = smooth(vb1[*,i],spikewd0,/edge_truncate)
        endfor
        store_data, pre0+'vb1_current_spin1', uts, vb1, limits = vlim
        
        ; the high freq part.
        get_data, pre0+'vb1_current_spin0', uts, tmp
        store_data, pre0+'vb1_current_spin2', uts, tmp-vb1, limits = vlim
        
        ; get the times for sectors.
        ; secuts are the start times of each sector.
        secuts = spinutr1[0]+dindgen(nsec0)*secdt0
        
        ; get the rotation matrix at the center time of each sector.
        secmuts = secuts+0.5*secdt0 ; middle time of a sector.
        tmp = time_string(secmuts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        tets = tet0+secmuts-secmuts[0]
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
            
        
        
        ; variables to be saved for all the sectors.
        ccinfos = replicate(ccinfo0, nsec0)
        
        get_data, pre0+'vb1_current_spin1', uts, vb1

        for i = 0, nsec0-1 do begin
            tutr = secuts[i]+[0,secdt0]
            idx = where(uts ge tutr[0] and uts le tutr[1], tnrec)
            i0 = idx[0]
            i1 = idx[tnrec-1]
            
            secvb1 = dblarr(tnrec,nvb1) ; save the vb1 in the current sector.
            ccs = dblarr(nccdt,nvb1)    ; save the cc curves b/w 1 and x.

            ; V1, remove trend.
            tmp = vb1[i0:i1,0]
            v1 = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
            secvb1[*,0] = v1
            
            ; loop through V[23456] and all the time lags, remove trend.
            for j = 1, nvb1-1 do begin
                for k = 0, nccdt-1 do begin
                    tmp = vb1[i0+ccdns[k]:i1+ccdns[k],j]
                    v2 = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
                    ccs[k,j] = c_correlate(v1,v2,0)
                endfor
            endfor
            
            ; pick out the time lag at max cc, remove trend.
            ; save the vb1 after time shift.
            for j = 1, nvb1-1 do begin
                ccinfos[i].cc1x[j] = max(ccs[*,j], idx)
                ccinfos[i].dn1x[j] = ccdns[idx]
                tmp = vb1[i0+ccdns[idx]:i1+ccdns[idx],j]
                secvb1[*,j] = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
            endfor
            store_data, pre0+'vb1_sec1', uts[i0:i1], secvb1, limits = vlim


            ; save the original vb1.
            for j = 0, nvb1-1 do begin
                tmp = vb1[i0:i1,j]
                secvb1[*,j] = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
            endfor
            store_data, pre0+'vb1_sec0', uts[i0:i1], secvb1, limits = vlim

            get_data, pre0+'vb1_current_spin2', uts, tmp
            store_data, pre0+'vb1_sec2', uts[i0:i1], tmp[i0:i1,*], limits = vlim



        ; **** determine timing for the current sector.
        ; correlation between 4 spin plane probes must be **high**,
        ; do not count spin axis boom.
        
        ; use the cross correlations in all the sectors to
        ; calc the plane normal, velocity in the half spin.
        
            get_data, pre0+'vb1_sec1', tmp, secvb1
            for j = 0, neb1-1 do $
                ccinfos[i].ccuvw[j] = c_correlate(secvb1[*,2*j],secvb1[*,2*j+1],0)
            ccinfos[i].dnuvw[0] = ccinfos[i].dn1x[1]
            ccinfos[i].dnuvw[1] = ccinfos[i].dn1x[3]-ccinfos[i].dn1x[2]
            ccinfos[i].dnuvw[2] = ccinfos[i].dn1x[5]-ccinfos[i].dn1x[4]
            
            idx = where(ccinfos[i].ccuvw[0:1] lt mincc0, cnt)
            if cnt eq 0 then begin
                ccinfos[i].dnmid = 0.5*(total(ccinfos[i].dn1x[1:3]))
            endif
            
            
            tinfo = ccinfos[i]
            tdn1x = tinfo.dn1x
            tmp = min(tinfo.cc1x[1:*], idx)
            idx = idx+1
            case idx of
                1: begin
                    tinfo.dnmid = 0.8*(tdn1x[2]+tdn1x[3])+0.2*(dn1x[4]+dn1x[5])
                    tdn1x[1] = tinfo.dnmid & end
                2: begin
                    tinfo.dnmid = 0.8*tdn1x[1]+0.2*(tdn1x[4]+tdn1x[5])
                    tdn1x[2] = tinfo.dnmid-tdn1x[3] & end
                3: begin
                    tinfo.dnmid = 0.8*tdn1x[1]+0.2*(tdn1x[4]+tdn1x[5])
                    tdn1x[3] = tinfo.dnmid-tdn1x[2] & end
                4: begin
                    tinfo.dnmid = 0.5*tdn1x[1]+0.5*(tdn1x[2]+tdn1x[3])
                    tdn1x[4] = tinfo.dnmid-tdn1x[5] & end
                5: begin
                    tinfo.dnmid = 0.5*tdn1x[1]+0.5*(tdn1x[2]+tdn1x[3])
                    tdn1x[5] = tinfo.dnmid-tdn1x[4] & end
            endcase
            
            tinfo.dnuvw = [tdn1x[1],tdn1x[3]-tdn1x[2],tdn1x[5]-tdn1x[4]]
            tinfo.dnuvw = [tdn1x[1],tdn1x[3]-tdn1x[2],0]
            tinfo.dnerr = max(abs(tdn1x-ccinfos.dn1x),/nan)
            tinfo.nhat = sunitvec(tinfo.dnuvw/boomlens)
            tinfo.vmag = abs(boomlens[0]/tinfo.dnuvw[0]*tinfo.nhat[0]/dr0*1e-3)
            ccinfos[i] = tinfo


        ; **** save plot for the cc on the current subsector of the current spin.
            if noplot eq 1 then continue
            vars = pre0+['vb1_sec'+['1','0','2']]
            nvar = n_elements(vars)
            figlabs = ['a. VB1 time shifted','b. VB1, low freq comp.', 'c. VB1, high freq comp.']
            
            
            ofn = figdir+'/sector_overview/rbsp'+tprobe+'/'+ $
                'cc3d_'+pre0+'spin'+time_string(spinutr0[0],tformat='hhmmss')+ $
                '_sec'+time_string(tutr[0],tformat='ss.fff')+'.pdf'
; ofn = 0
            sgopen, ofn, xsize = 11, ysize = 8.5, /inch
            device, decomposed = 0
            loadct2, 43
            
            !x.gridstyle = 1
            !x.ticklen = 1
            
            poss = sgcalcpos(nvar,lmargin = 25, tmargin = 10)

            xchsz = double(!d.x_ch_size)/!d.x_size
            ychsz = double(!d.y_ch_size)/!d.y_size
            
            tplot, vars, position = poss, vlab_marg = 20
            
            
            tx = xchsz*5
            for j = 0, nvar-1 do begin
                ty = poss[3,j]-ychsz*1
                xyouts, tx, ty, /normal, figlabs[j], alignment = 0
            endfor
            
            xyouts, 0.5, 0.95, /normal, alignment = 0.5, charsize = 1.25, $
                'RBSP-'+strupcase(tprobe)+' VB1 of '+secdt0str+' sec from '+ $
                time_string(tutr[0], tformat='hh:mm:ss.fff')+', in the spin from '+ $
                time_string(spinutr0[0], tformat='hh:mm:ss')
            
            
            tx = poss[0,0]
            txs = tx+15*xchsz+findgen(nvb1-1)*5*xchsz

            ty = 0.95-ychsz*(2+1.2*0)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Cross corr. b/w V1 and'
            for j = 1, nvb1-1 do xyouts, txs[j-1], ty, /normal, alignment = 0, $
                'V'+sgnum2str(j+1)
                
            ty = 0.95-ychsz*(2+1.2*1)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Max cross corr.'
            for j = 1, nvb1-1 do xyouts, txs[j-1], ty, /normal, alignment = 0, $
                sgnum2str(ccinfos[i].cc1x[j],ndec=2)
            
            ty = 0.95-ychsz*(2+1.2*2)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Time lag (msec)'
            for j = 1, nvb1-1 do xyouts, txs[j-1], ty, /normal, alignment = 0, $
                sgnum2str(ccinfos[i].dn1x[j]*dr0*1e3,ndec=2)
            
            ty = 0.95-ychsz*(2+1.2*3)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Shift in # of record'
            for j = 1, nvb1-1 do xyouts, txs[j-1], ty, /normal, alignment = 0, $
                sgnum2str(ccinfos[i].dn1x[j])
            
            
            tx = poss[0,0]+(poss[2,0]-poss[0,0])*0.6
            txs = tx+15*xchsz+findgen(neb1)*5*xchsz

            ty = 0.95-ychsz*(2+1.2*0)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Cross corr. along '
            for j = 0, neb1-1 do xyouts, txs[j], ty, /normal, alignment = 0, $
                'V'+sgnum2str(j*2+1)+sgnum2str(j*2+2)

            ty = 0.95-ychsz*(2+1.2*1)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Max cross corr.'
            for j = 0, neb1-1 do xyouts, txs[j], ty, /normal, alignment = 0, $
                sgnum2str(ccinfos[i].ccuvw[j],ndec=2)
            
            ty = 0.95-ychsz*(2+1.2*2)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Time lag (msec)'
            for j = 0, neb1-1 do xyouts, txs[j], ty, /normal, alignment = 0, $
                sgnum2str(ccinfos[i].dnuvw[j]*dr0*1e3,ndec=2)
            
            ty = 0.95-ychsz*(2+1.2*3)
            xyouts, tx, ty, /normal, alignment = 0, $
                'Shift in # of record'
            for j = 0, neb1-1 do xyouts, txs[j], ty, /normal, alignment = 0, $
                sgnum2str(ccinfos[i].dnuvw[j])
            
            !x.gridstyle = 0
            !x.ticklen = 0
            
            sgclose
        endfor
        
        normangs = atan(ccinfos.nhat[1],ccinfos.nhat[0])*deg
        for j = 1, nsec0-1 do begin
            dangs = normangs[j]+[-1,0,1]*360
            tmp = min(dangs-normangs[j-1], idx, /absolute)
            normangs[j] = dangs[idx]
        endfor
        
        spinangs = findgen(nsec0)/(nsec0-1)*5/spinrate0*360
        dangs = normangs+spinangs
        idx = where(abs(dangs-mean(dangs,/nan)) le 0.8*stddev(dangs,/nan),cnt); remove points have large fluctuation.
        dangs = dangs[idx]
        normang0 = mean(dangs,/nan)
        normangerr = stddev(dangs,/nan)
        
        vmags = ccinfos.vmag
        vmags = vmags[where(finite(vmags))]
        idx = where(abs(vmags-mean(vmags)) le 0.8*stddev(vmags),cnt); remove points have large fluctuation.
        vmag0 = mean(vmags[idx],/nan)
        vmagerr = stddev(vmags[idx],/nan)

        ; convert the normal from uvw to gsm.        
        nuvw = [cos(normang0*rad),sin(normang0*rad),0]
        tsec0 = spinutr1[0]+halfspinsecdt0/nsec0*0.5
        tmp = time_string(tsec0,tformat='YYYY-MM-DDThh:mm:ss.ffffff')
        cspice_str2et, tmp, tet0
        
        scid = strupcase(pre0+'science')
        cspice_pxform, scid, 'GSM', tet0, muvw2gsm

        ngsm = nuvw # muvw2gsm
        
        
        ofn = figdir+'/spin_overview/rbsp'+tprobe+'/'+ $
                'cc3d_'+pre0+'spin'+time_string(spinutr0[0],tformat='hhmmss')+'.pdf'
; ofn = 0
        sgopen, ofn, xsize = 11, ysize = 8.5, /inch
        device, decomposed = 0
        loadct2, 43
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        poss = sgcalcpos(1,2, tmargin = 5.5, bmargin = 20, xpad = 15)
        poss[*,1]+= [-25*xchsz,0.30,0,0]
        xr = mean(spinangs)+[-1,1]*90
        yr = normang0+[-270,90]

        plot, spinangs, normangs, psym = 1, position = poss[*,0], /noerase, $
            xtitle = 'Spin angle!C(deg)', xrange = xr, xstyle = 1, $
            ytitle = 'Normal angle!C(deg)', yrange = yr, ystyle = 1
        oplot, xr, -xr+normang0, color = 6
        
        
        xr = [-1, nsec0]
        yr = [0,vmag0+vmagerr*2]
        plot, ccinfos.vmag, psym = 1, position = poss[*,1], /noerase, $
            xrange = xr, xstyle = 1, xtitle = 'Sector ID', $
            yrange = yr, ystyle = 1, ytitle = '|V|!Dnormal!N!C(km/s)'
        plots, xr, vmag0+[0,0], color = 6
        
        xyouts, 0.5, 0.95, /normal, alignment = 0.5, charsize = 1.25, $
            'RBSP-'+strupcase(tprobe)+' half spin from '+ $
            time_string(spinutr0[0], tformat='hh:mm:ss.fff')+ $
            ' to '+time_string(spinutr0[1], tformat='ss.fff')
        
        tx = poss[0,0]
        ty = 0.95-ychsz*(2+1.2*0)
        xyouts, tx, ty, /normal, alignment = 0, $
            'Despun normal angle at sector0 = '+sgnum2str(round(normang0))+' -/+ '+ $
            sgnum2str(round(normangerr))+' deg'
        
        tx = poss[0,1]
        ty = 0.95-ychsz*(2+1.2*0)
        xyouts, tx, ty, /normal, alignment = 0, $
            'Normal |V| = '+sgnum2str(vmag0,ndec=1)+' -/+ '+ $
            sgnum2str(vmagerr,ndec=1)+' km/s'
        

        tpos = poss[*,1] & tpos[1] = poss[1,0] & tpos[3] = poss[1,1]-ychsz*5
        tpos[2] = tpos[0]+(tpos[3]-tpos[1])*!d.y_size/!d.x_size
        plot, [-1,1],[-1,1], /noerase, position = tpos, /normal, /nodata, $
            xtitle = 'GSM Y', ytitle = 'GSM Z', xrange = [1,-1], xstyle = 1, $
            ystyle = 1
        arrow, 0,0, ngsm[1],ngsm[2], /data, /solid, hsize = !d.x_ch_size
        
        tx = tpos[2]+xchsz*3
        ty = tpos[3]-ychsz*(1+1.2*0)
        xyouts, tx, ty, /normal, alignment = 0, $
            'Normal direction in GSM: ('+strjoin(string(ngsm,format='(F5.2)'),',')+')'
        
        ty = tpos[3]-ychsz*(1+1.2*1)
        xyouts, tx, ty, /normal, alignment = 0, $
            'Normal direction in UVW: ('+strjoin(string(nuvw,format='(F5.2)'),',')+')'
        
        ty = tpos[3]-ychsz*(1+1.2*3)
        xyouts, tx, ty, /normal, alignment = 0, $
            'U in GSM: ('+strjoin(string(muvw2gsm[0,*],format='(F5.2)'),',')+')'
            
        ty = tpos[3]-ychsz*(1+1.2*4)
        xyouts, tx, ty, /normal, alignment = 0, $
            'V in GSM: ('+strjoin(string(muvw2gsm[1,*],format='(F5.2)'),',')+')'
        
        ty = tpos[3]-ychsz*(1+1.2*5)
        xyouts, tx, ty, /normal, alignment = 0, $
            'W in GSM: ('+strjoin(string(muvw2gsm[2,*],format='(F5.2)'),',')+')'        


        tpos = [poss[0,0],0.1,poss[2,1],poss[1,0]-ychsz*5]
        tplot, pre0+'vb1_low_res', trange = utr, position = tpos, /noerase
        timebar, spinutr0, color = 6
; stop
        sgclose
    endforeach
    
endforeach



end
