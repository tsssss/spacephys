;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    secofday = 86400d


;---Settings.
    utr1 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    utr4 = time_double(['2014-12-22/04:00','2014-12-22/06:30'])
    pres = ['rba','rbb','tha','the','thd']+'_'
    sites = ['whit','chbg']
    
    
    ;---Use this block to determine the mlat to calculate ewogram.
;        dval0 = 150
;        ; RB.
;        ; low latitude structure 62.5-63.0 deg.
;        ; high latitude structure 63.5-64.0 deg.
;        ; TH.
;        ; 64.0-64.5 deg.
;        idx = where(uts ge time_double('2014-12-22/05:00') and uts le time_double('2014-12-22/05:40'))
;        for i=idx[0],idx[-1],5 do begin
;            ;tv, bytscl(imgs[*,*,i+1]-imgs[*,*,i], min=-dval0, max=dval0)
;            timg = imgs[*,*,i]
;            idx = where(mlats ge 64.0 and mlats le 64.5, complement=idx2)
;            timg[idx2] = 0
;            tv, bytscl(timg, min=10, max=500, top=254)
;            wait, 0.2
;        endfor
;        stop


    ;_2014_1222_load_data
    

;---Calc auroral motion.
; step down: for aurora moves eastward,
; step up: for aurora moves westward.
; max: for the brightest aurora.
    infos = [$
        {var:'rbx_ewo_63.5', method:'stepdown', val0:300d, $
        utr:time_double(['2014-12-22/05:20','2014-12-22/05:40']),mltrng:[0.3,1.7]}, $
        {var:'thx_ewo_64.0', method:'max', val0:250d, $
        ;utr:time_double(['2014-12-22/05:21','2014-12-22/05:26']),mltrng:[-5.3,-5.0]}]
        utr:time_double(['2014-12-22/05:17','2014-12-22/05:26']),mltrng:[-5.3,-4.5]}, $
        {var:'thx_ewo_64.0', method:'stepup', val0:250d, $
        utr:time_double(['2014-12-22/05:11','2014-12-22/05:15']),mltrng:[-6.0,-5.2]}]

    ofn = 0
    ofn = shomedir()+'/fig_asi_motion.pdf'
    sgopen, ofn, xsize=6, ysize=4, /inch
    
    device, decomposed=0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    poss = sgcalcpos(2, tmargin=1, lmargin=13, rmargin=9, bmargin=4)
    thick = (size(ofn,/type) eq 7)? 4: 2
    
    foreach tinfo, infos, k do begin 
        tvar = tinfo.var
        tutr = tinfo.utr
        vrng = tinfo.mltrng
        tmethod = tinfo.method
        val0 = tinfo.val0
        
        options, tvar, 'ztitle', 'Photon Count'
        options, tvar, 'ytitle', 'MLT (hr)!C'+strmid(tvar,3,/reverse)+' deg'
        ttitle = 'Eowgram!C    '
        if k eq 0 then begin
            tpos = poss[*,0]
            tplot, tvar, position=tpos, /noerase, /novtitle, /nouttick
            xyouts, tpos[0]-xchsz*12, tpos[3]-ychsz*1, /normal, 'a. '+ttitle+'post-midn'
        endif else begin
            tpos = poss[*,1]
            if k eq 1 then tplot, tvar, position=tpos, /noerase, /novtitle
            xyouts, tpos[0]-xchsz*12, tpos[3]-ychsz*1, /normal, 'b. '+ttitle+'pre-midn'
        endelse
        
        
        get_data, tvar, uts, dat, val, limits=lim
        idx = where(uts ge tutr[0] and uts le tutr[1],nrec)
        uts = uts[idx]
        dat = dat[idx,*]
        idx = where(val ge vrng[0] and val le vrng[1])
        val = val[idx]
        dat = dat[*,idx]
        mlts = dblarr(nrec)
        for i=0, nrec-1 do begin
            tdat = dat[i,*]
            if tmethod eq 'stepdown' then begin
                ; remove noise.
                idx = where(val ge 0.0032*i+0.5, cnt)
                if cnt ne 0 then tdat[idx] = 0
                ; detect boundary.
                idx = where(tdat ge val0, cnt)
                mlts[i] = (cnt eq 0)? !values.d_nan: val[idx[cnt-1]]
            endif else if tmethod eq 'stepup' then begin
                ; detect boundary
                idx = where(tdat ge val0, cnt)
                tmp = max(tdat[idx[0]:*])
                idx = where(tdat eq tmp)
                mlts[i] = (cnt eq 0)? !values.d_nan: val[idx[0]]
            endif else if tmethod eq 'max' then begin
                ; remove noise
                idx = where(val le -0.004*i-4.8, cnt)
                if cnt ne 0 then tdat[idx] = 0
                ; detect boundary
                idx = where(tdat ge val0, cnt)
                tmp = max(tdat[idx[0]:*])
                idx = where(tdat eq tmp)
                mlts[i] = (cnt eq 0)? !values.d_nan: val[idx[0]]
            endif
        endfor
        idx = uniq(mlts)
        mlts = mlts[idx]
        uts = uts[idx]
        idx = where(finite(mlts))
        mlts = mlts[idx]
        uts = uts[idx]
        ;plot, uts, mlts, psym=1
        res = linfit(uts, mlts)
        
        tx = tpos[0]+xchsz*1
        ty = tpos[3]-ychsz*1
        if k eq 2 then ty -= ychsz*1
        xyouts, tx, ty, /normal, 'V'+string(k,format='(I0)')+'='+$
            string((res[1]*15),format='(F6.3)')+' deg/sec', color=255
        
        plot, utr4, lim.yrange, /noerase, /nodata, position=tpos, xstyle=5, ystyle=5
        plots, uts, mlts, psym=1, symsize=0.5, color=255
        oplot, utr4, res[1]*utr4+res[0], color=255, thick=thick
    endforeach
    
    
    sgclose

end
