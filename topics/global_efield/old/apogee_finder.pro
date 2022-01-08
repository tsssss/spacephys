;+
; Find apogee position for given time
; posvar is any 3-D vector whose magnitude gives distance, or the distance directly.
; more_vars can contain the mlt in hr.
;-

pro apogee_finder, posvar, period=period, apogee=apogee, more_vars=mvars

    if n_elements(posvar) eq 0 then message, 'no pos info ...'
    get_data, posvar, uts, pos, limits=lim
    nrec = n_elements(uts)
    sz = size(pos,/dimensions)
    if n_elements(sz) eq 1 then dis = pos else dis = snorm(pos)

    ; find apogee using derivative.
    df = dis[1:nrec-1]-dis[0:nrec-2]
    idx = where(df le 0, cnt)
    if cnt eq 0 then begin
        message, 'no apogee found ...', /continue
        return
    endif

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
    tvar = posvar+'_apogee'
    store_data, tvar, uts[idx1], dis[idx1], limits=lim
    
    period = sdatarate(uts[idx1])
    apogee = mean(dis[idx1])
    
    foreach tvar, mvars do begin
        get_data, tvar, tuts, tdat, limits=lim
        store_data, tvar+'_apogee', tuts[idx1], tdat[idx1], limits=lim
    endforeach

end

;---use on themis.
    tplot_options, 'ynozero', 1
    tplot_options, 'lab_num_min', 6
    tplot_options, 'labflag', -1
    tplot_options, 'yticklen', -0.01
    tplot_options, 'xticklen', -0.02

    utr0 = time_string(['2008-01-01','2015-12-31'])
    posvar0 = ['Epoch','RADIUS','SM_LCT_T']
    posvar1 = ['epoch','radius','mlt']
    probes = ['a','d','e','b','c']
    probes = ['b','c']
    foreach tprobe, probes do begin
        pre0 = 'th'+tprobe+'_'
        if tnames(pre0+'dis') eq '' then begin
            pos = sread_thm_orbit(utr0, vars=posvar0, newname=posvar1, probe=tprobe)
            uts = sfmepoch(pos.epoch[0:*:10],'unix')
            tvar = pre0+'dis'
            store_data, tvar, uts, pos.radius[0:*:10], limits={ytitle:'(Re)',labels:'Apogee'}
            tvar = pre0+'mlt'
            store_data, tvar, uts, pos.mlt[0:*:10], limits={ytitle:'(hr)',labels:'MLT',yrange:[0,24],$
                ystyle:1,constant:[6,12,18],yticks:4,yminor:6}
        endif
        apogee_finder, pre0+'dis', more_vars=pre0+['mlt']
        
        ofn = shomedir()+'/fig_th'+tprobe+'_apogee_history.pdf'
        sgopen, ofn, xsize=5, ysize=3, /inch
        pos0 = [0.15,0.15,0.85,0.85]
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        vars = pre0+['dis','mlt']+'_apogee'
        nvar = n_elements(vars)
        poss = sgcalcpos(nvar, position=pos0)
        tplot, vars, trange=utr0, position=poss
        titl = 'TH-'+strupcase(tprobe)+' apogee history'
        xyouts, (pos0[0]+pos0[2])*0.5, pos0[3]+ychsz*0.5, /normal, titl, alignment=0.5
        sgclose
    endforeach
    
    ofn = shomedir()+'/thm_apogee_history.sav'
    vars = tnames('th?_'+['dis','mlt']+'_apogee')
    tplot_save, vars, filename=ofn


stop

;---test on fake data.
    nrec = 1000
    txs = findgen(nrec)
    tys = sin(txs/(nrec-1)*!dpi*10)

    tvar = 'test_pos'
    store_data, tvar, txs, tys
    apogee_finder, tvar

end
