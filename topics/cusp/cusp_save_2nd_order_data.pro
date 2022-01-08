;+
; Read and save second order data, including
;   S, dE, dB and related FFT and wavelet results.
;-

function exponent, axis, index, number

    ; a special case.
    if number eq 0 then return, '0'
    
    ; assuming multiples of 10 with format.
    ex = string(number, format='(e8.0)')
    pt = strpos(ex, '.')
    
    first = strmid(ex, 0, pt)
    sign = strmid(ex, pt+2, 1)
    thisexponent = strmid(ex, pt+3)
    
    ; shave off leading zero in exponent
    while strmid(thisexponent, 0, 1) eq '0' do thisexponent = strmid(thisexponent, 1)
    
    ; fix for sign and missing zero problem.
    if (long(thisexponent) eq 0) then begin
        sign = ''
        thisexponent = '0'
    endif
    
    ; make the exponent a superscript.
    if sign eq '-' then begin
        return, sign + thisexponent
    endif else begin
        return, thisexponent
    endelse
    
end


pro cusp_save_2nd_order_data, eventid, test=test, $
    no_plot=noplot, save_data=save_data
    
    
    if n_elements(eventid) eq 0 then message, 'no event id ...'
    
;---constant.
    re = 6378d & re1 = 1d/re
    
;---settings.
    ; rootdir to save the data and plots.
    rootdir = sdiskdir('GoogleDrive')+'/My Drive/works'
    if file_test(rootdir) eq 0 then rootdir = shomedir()

    ; dir to save plots and data.
    figdir = rootdir+'/works/cusp/cusp list conjun'
    datdir = rootdir+'/data/cusp'
    if ~file_test(figdir,/directory) then file_mkdir, figdir
    if ~file_test(datdir,/directory) then file_mkdir, datdir

    ; to prevent overwriting useful data without confirm.
    if keyword_set(test) then begin
        figdir = shomedir()+'/cusp'
        datdir = shomedir()+'/cusp/data'
    endif


    deidx = 0   ; v: north-south.
    dbidx = 1   ; p: east-west.
    pfidx = 2   ; b: parallel.

    sats = ['po','fa']
    
    ; read basic info from list of conjunc events.
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event=eventid)
    if size(loginfo,/type) ne 8 then message, 'no id found ...'
    
    poinfo = loginfo.polar
    fainfo = loginfo.fast
    
    ; plot time range.
    poplotutr = poinfo.plot_time
    faplotutr = fainfo.plot_time
    if poplotutr[1] lt poplotutr[0] then poplotutr[1]+= 86400d
    if faplotutr[1] lt faplotutr[0] then faplotutr[1]+= 86400d
    
    ; cusp time range.
    pocusputr = poinfo.cusp_time
    facusputr = fainfo.cusp_time
    if pocusputr[1] lt pocusputr[0] then pocusputr[1]+= 86400d
    if facusputr[1] lt facusputr[0] then facusputr[1]+= 86400d
    
    poplotutr = mean(pocusputr)+[-1,1]*1.6*abs(pocusputr[1]-pocusputr[0])
    faplotutr = mean(facusputr)+[-1,1]*1.6*abs(facusputr[1]-facusputr[0])
    
    
    posu = [0,0.4,1,1]
    posd = [0,0,1,0.4]
    
    labfac = ['n-s','e-w','b']
    
    zchsz = 0.8
    
    red = 6
    blue = 2
    white = 255
    black = 0


;---plot settings and generate plots to disk.
    posl = [0.1,0.10,0.40,0.90]
    posr = [0.6,0.10,0.90,0.90]
    faclabs = ['v','p','b']
    falabs = 'fa_'+['ilat','mlt','dis']
    polabs = 'po_'+['ilat','mlt','dis']
    rgb = [6,4,2] & red = 6
    ct = 43
    charsz = 1
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'zcharsize', zchsz
    time_stamp, /off
    
    

;---load 1st order data.
    ifn = datdir+'/'+eventid+'_1st_order_data.tplot'
    if file_test(ifn) eq 0 then begin
        cusp_save_1st_order_data, eventid, /save_data
    endif
    store_data, '*', /delete
    tplot_restore, filename=ifn

    infofn = datdir+'/'+eventid+'_scinfo.tplot'
    tplot_restore, filename=infofn


;---calculate data and save to tplot.
; [po,fa]_[de,db,pf]_fac_fft_info
; [po,fa]_pf_fac
    get_data, 'scinfo', tutr, scinfo
    foreach tsat, sats do begin
        pre0 = tsat+'_'

        ; instataneous pflux.
        scaleinfo = {s0:0d,s1:0d,dj:1d/8,ns:0d}
        stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', $
            scaleinfo = scaleinfo
        
        if tsat eq 'po' then begin
            scinfo.polar.scaleinfo = scaleinfo
        endif else begin
            scinfo.fast.scaleinfo = scaleinfo
        endelse
        
        vars = pre0+['de','db','pf']
        idxs = [deidx,dbidx,pfidx]
        for i=0, 2 do begin
            tvar = vars[i]+'_mor'
            stplot_index, vars[i]+'_fac', idxs[i], newname=tvar
            stplot_mor, vars[i]+'_fac', newname=tvar
        endfor
        
        cusputr = (tsat eq 'po')? pocusputr: facusputr
        options, pre0+'pf_mor', 'constant', (cusputr[1]-cusputr[0])
    endforeach

    infofn = datdir+'/'+eventid+'_scinfo.tplot'
    store_data, 'scinfo', tutr, scinfo
    tplot_save, 'scinfo', filename = infofn

;    device, decomposed=0
;    loadct2, 43
;    erase
;    
;    poss = sgcalcpos(1,2)
;    
;    tpos = poss[*,0]
;    tplot, 'po_pf_mor', position=tpos, trange=poplotutr, /noerase
;    timebar, pocusputr, color=6
;    
;    tpos = poss[*,1]
;    tplot, 'fa_pf_mor', position=tpos, trange=faplotutr, /noerase
;    timebar, facusputr, color=6



;---save data to disk.
; save labels, de_fac, db_fac, ele_keflux, ion_keflux, in situ.
; save event info to (1) include other info and (2) info for calc
; quantities are not saved (e.g., pflux).
    if keyword_set(save_data) then begin
        vars = ['pf_fac',['pf','de','db']+'_mor_fft_info']
        vars = ['po_'+vars,'fa_'+vars]
        ofn = datdir+'/'+eventid+'_2nd_order_data.tplot'
        tplot_save, vars, filename = ofn
    endif
    while !d.window ne -1 do wdelete, !d.window



    if keyword_set(noplot) then return

;---plot 1 & 2: polar E, B, S.    
    foreach tsat, sats do begin
        
        pre0 = tsat+'_'
        plotutr = (tsat eq 'po')? poplotutr: faplotutr
        cusputr = (tsat eq 'po')? pocusputr: facusputr
        
        ofn = 0
        ofn = figdir+'/'+eventid+'/'+eventid+'_'+tsat+'_field.pdf'
        sgopen, ofn, xsize=6, ysize=8, /inch
        
        device, decomposed=0
        loadct2, 43
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        poss = sgcalcpos(3, position=[0,0.1,1,1], xpad=0, ypad=0)
        
        vars = pre0+['de','db','pf']
        figlabs = ['a','b','c']
        units = ['(mV/m)','(nT)','(mW/m!U2!N)']
        ytitl = ['dE!C(mV/m)','dB!C(nT)','S!C(mW/m!U2!N)']
        
        zlog = 1
        !x.ticklen = -0.03
        !y.ticklen = -0.01
        
        options, vars+'_mor', 'constant', cusputr[1]-cusputr[0]
        
        for i=0, 1 do begin
            tpos = sgcalcpos(1, region=poss[*,i], bmargin=0)
            tunit = units[i]
            
            xsize = tpos[2]-tpos[0]
            ysize = tpos[3]-tpos[1]
            
            ; spectrogram.
            pos1 = [tpos[0],tpos[3]-ysize*0.6,tpos[0]+xsize*0.75,tpos[3]]
            ; colorbar.
            pos2 = [pos1[0],pos1[3]+ychsz*0.2,pos1[2],pos1[3]+ychsz*0.7]
            ; spectrum.
            pos3 = [pos1[2]+xchsz*0.8,pos1[1],tpos[2],pos1[3]]
            ; waveform.
            pos4 = [pos1[0],tpos[1],pos1[2],pos1[1]-ychsz*0.2]
            
            get_data, vars[i]+'_mor_fft_info', tmp, info
            get_data, vars[i]+'_mor', uts, dat
            
            ; set options for spectrogram.
            zrng = 10d^double(ceil(alog10([median(dat),max(dat)])))
            tvar = vars[i]+'_mor'
            options, tvar, 'yrange', minmax(info.ps)
            options, tvar, 'ystyle', 1
            options, tvar, 'zposition', pos2
            options, tvar, 'zhorizontal', 1
            options, tvar, 'ztitle', tunit+'!U2'
            options, tvar, 'zlog', zlog
            options, tvar, 'zrange', zrng
            tplot, vars[i]+'_mor', position=pos1, /noerase, trange=plotutr, $
                /nouttick
            plot, uts, info.coi, color=white, /noerase, position=pos1, $
                xstyle=5, xrange=plotutr, $
                ystyle=5, yrange=minmax(info.ps), ylog=1
            timebar, cusputr, color=red
            xyouts, pos1[0]-xchsz*8, pos1[3]-ychsz*0.5, /normal, $
                figlabs[i]+'-1.'
            
            tvar = vars[i]+'_fac'
            if i eq 0 then begin
                get_data, tvar, uts, dat
                dwid = (tsat eq 'po')? 6d: 3
                drec = dwid/sdatarate(uts)
                for j=0, 2 do dat[*,j] = smooth(dat[*,j], drec, /nan)
                store_data, tvar, uts, dat
            endif
            options, tvar, 'yticks', 3
            options, tvar, 'labels', labfac
            options, tvar, 'ytitle', ytitl[i]
            tplot, tvar, position=pos4, /noerase, trange=plotutr, /nouttick
            if i eq 0 then xyouts, pos4[2]-xchsz*0.5, pos4[1]+ychsz*0.2, /normal, $
                alignment=1, sgnum2str(dwid)+' sec smooth', charsize=0.8
            timebar, cusputr, color=red
            xyouts, pos4[0]-xchsz*8, pos4[3]-ychsz*0.5, /normal, $
                figlabs[i]+'-2.'
            
            xrng = alog10(minmax(info.gws))
            xrng = (floor(xrng)+[0,1])
            xtickv = findgen(xrng[1]-xrng[0])+xrng[0]+1
            xtickv = 10d^(xtickv[0:*:2])
            xticks = n_elements(xtickv)-1
            xrng = 10d^xrng
            plot, info.gws, info.ps, xlog=1, ylog=1, position=pos3, /noerase, $
                ystyle=1, ytitle='', ytickformat='(A1)', yrange=minmax(info.ps), $
                xstyle=1, xtitle='', xtickformat='(A1)', xrange=xrng
            axis, xaxis=1, xtitle='Log '+tunit+'!U2', xstyle=1, xlog=1, xcharsize=zchsz, $
                xrange=xrng, xtickv=xtickv, xticks=xticks, xtickformat='exponent'
            tpos = pos3
            xyouts, tpos[2]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*0, /normal, $
                figlabs[i]+'-3.'
            
            tmp = linfit(alog(info.ps),alog(info.gws))
            xyouts, pos3[0]+xchsz*1, pos3[3]-ychsz*1.2, /normal, $
                'PS~f!U -'+sgnum2str(tmp[1],ndec=1), charsize=zchsz
        endfor
        
        
        ; pflux.
        i = 2
        tpos = sgcalcpos(1, region=poss[*,i], bmargin=0)
        tunit = 'Log (mW/m!U2!N)'
        
        xsize = tpos[2]-tpos[0]
        ysize = tpos[3]-tpos[1]
        
        ; spectrogram.
        pos1 = [tpos[0],tpos[3]-ysize*0.6,tpos[0]+xsize*0.75,tpos[3]]
        ; colorbar.
        pos2 = [pos1[0],pos1[3]+ychsz*0.2,pos1[2],pos1[3]+ychsz*0.7]
        ; spectrum.
        pos3 = [pos1[2]+xchsz*0.8,pos1[1],tpos[2],pos1[3]]
        ; waveform.
        pos4 = [pos1[0],tpos[1],pos1[2],pos1[1]-ychsz*0.2]
        
        get_data, pre0+'de_mor_fft_info', tmp, deinfo
        get_data, pre0+'db_mor_fft_info', tmp, dbinfo
        
        tvar = vars[i]+'_mor'
        get_data, tvar, uts, dat, val
        dat = sqrt(dat)
        store_data, tvar, uts, dat, val
        zrng = 10d^double(ceil(alog10([median(dat),max(dat)])))
        options, tvar, 'zposition', pos2
        options, tvar, 'zhorizontal', 1
        options, tvar, 'ztitle', '(mW/m!U2!N)'
        options, tvar, 'zlog', zlog
        options, tvar, 'zrange', zrng
        tplot, vars[i]+'_mor', position=pos1, /noerase, trange=plotutr, $
            /nouttick
        plot, uts, info.coi, color=white, /noerase, position=pos1, $
            xstyle=5, xrange=plotutr, $
            ystyle=5, yrange=minmax(info.ps), ylog=1
        timebar, cusputr, color=red
        xyouts, pos1[0]-xchsz*8, pos1[3]-ychsz*0.5, /normal, $
            figlabs[i]+'-1.'
        
        tvar = vars[i]+'_fac'
        options, tvar, 'yticks', 3
        options, tvar, 'labels', ['','','']
        options, tvar, 'ytitle', ytitl[i]
        tplot, tvar, position=pos4, /noerase, trange=plotutr
        timebar, cusputr, color=red
        xyouts, pos4[0]-xchsz*8, pos4[3]-ychsz*0.5, /normal, $
            figlabs[i]+'-2.'
        
        xrng = alog10(minmax(sqrt(info.gws)))
        xrng = (floor(xrng)+[0,1])
        xtickv = findgen(xrng[1]-xrng[0])+xrng[0]+1
        xtickv = 10d^(xtickv[0:*:2])
        xticks = n_elements(xtickv)-1
        xrng = 10d^xrng
        plot, sqrt(info.gws), info.ps, xlog=1, ylog=1, position=pos3, /noerase, $
            ystyle=1, ytitle='', ytickformat='(A1)', yrange=minmax(info.ps), $
            xstyle=5, xtitle='', xtickformat='(A1)', xrange=xrng
        axis, xaxis=1, xtitle=tunit, xstyle=1, xlog=1, xcharsize=zchsz, $
            xrange=xrng, xtickv=xtickv, xticks=xticks, xtickformat='exponent'
        
        
        tunit = 'Log (km/s)'
        
        get_data, pre0+'de_mor_fft_info', tmp, deinfo
        get_data, pre0+'db_mor_fft_info', tmp, dbinfo
        txs = sqrt(deinfo.gws/dbinfo.gws)*1e3
        xrng = alog10(minmax(txs))
        xrng = (floor(xrng)+[0,1])
        xtickv = findgen(xrng[1]-xrng[0])+xrng[0]+1
        xtickv = 10d^(xtickv[0:*:2])
        xticks = n_elements(xtickv)-1
        xrng = 10d^xrng
        
        plot, txs, deinfo.ps, $
            position=pos3, /noerase, color=6, $
            ystyle=5, ylog=1, ytitle='', ytickformat='(A1)', yrange=minmax(info.ps), $
            xstyle=5, xlog=1, xtitle='', xtickformat='(A1)', xrange=xrng
        axis, xaxis=0, xtitle=tunit, xstyle=1, xlog=1, xcharsize=zchsz, $
            xrange=xrng, xtickv=xtickv, xticks=xticks, xtickformat='exponent'
        
        ; model Va.
        get_data, 'scinfo', tmp, scinfo
        tinfo = (tsat eq 'po')? scinfo.polar: scinfo.fast
        va = tinfo.va
        va*= sqrt(2)   ; account for diff b/w val & max.
        oplot, va+[0,0], minmax(info.ps), linestyle=1, color=red
        oplot, va*0.25+[0,0], minmax(info.ps), linestyle=2, color=red
        
        tpos = pos3
        xyouts, tpos[2]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*0, /normal, $
            figlabs[i]+'-3.'
        xyouts, tpos[2]+xchsz*0.5, tpos[3]-ychsz*1-ychsz*1, /normal, $
            'PS(S)', charsize=zchsz
        xyouts, tpos[2]+xchsz*0.5, tpos[3]-ychsz*1-ychsz*2, /normal, $
            'R(E/B)', charsize=zchsz
        xyouts, tpos[2]+xchsz*0.5, tpos[3]-ychsz*1-ychsz*3, /normal, $
            'V!DA!N(O)', charsize=zchsz
        xyouts, tpos[2]+xchsz*0.5, tpos[3]-ychsz*1-ychsz*4, /normal, $
            'V!DA!N(H)', charsize=zchsz
            
        plots, tpos[2]+xchsz*[4,6], tpos[3]-ychsz*0.8-ychsz*1, /normal, $
            color=black, linestyle=0
        plots, tpos[2]+xchsz*[4,6], tpos[3]-ychsz*0.8-ychsz*2, /normal, $
            color=red, linestyle=0
        plots, tpos[2]+xchsz*[4,6], tpos[3]-ychsz*0.8-ychsz*3, /normal, $
            color=red, linestyle=2
        plots, tpos[2]+xchsz*[4,6], tpos[3]-ychsz*0.8-ychsz*4, /normal, $
            color=red, linestyle=1
        sgclose
    endforeach
    
end

;id = '1998_0925_05'
;id = '1998_1001_02'
;cusp_save_2nd_order_data, id, /save_data

;ids = cusp_id('all')
;foreach id, ids do begin
;    print, '----processing '+id+' ...'
;    cusp_save_2nd_order_data, id, /save_data
;endforeach
;end
;

id = '1998_0925_05'
id = '1999_1010_15'
cusp_save_2nd_order_data, id
end