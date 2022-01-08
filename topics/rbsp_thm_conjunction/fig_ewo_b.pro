;+
; Caclulate the auroral azimuthal velocity from single camera.
;-

;---Constant.
    secofday = 86400d
    top = 254
    

;---Settings.
    utr0 = time_double(['2014-08-27/08:40','2014-08-27/10:10'])
    utr0 = time_double(['2014-08-27/08:50','2014-08-27/10:10'])
    sites = ['atha','pina']
    rootdir = shomedir()+'/Google Drive/works/works/rbsp_thm_conjunction'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir

    minelev = 5d
    height = 110d       ; km.
    hgtidx = where([90,110,150] eq height)
    dmlt = 0.05         ; hr. resolution in mlt.

    
    mlatrng = [62d,63]  ; gen ewgram for this mlat range.
    ;mlatrng = [63d,64]  ; gen ewgram for this mlat range.
    mlatrng = [64d,65]  ; gen ewgram for this mlat range.
    mlatrng = [65d,66]  ; gen ewgram for this mlat range.
    datfn = rootdir+'/ewo_'+time_string(utr0[0],tformat='YYYY_MMDD')+'_'+sgnum2str(mlatrng[0])+'-'+sgnum2str(mlatrng[1])+'deg.tplot'
    figfn = rootdir+'/ewo_'+time_string(utr0[0],tformat='YYYY_MMDD')+'_'+sgnum2str(mlatrng[0])+'-'+sgnum2str(mlatrng[1])+'deg.pdf'

    
    load_ewo = 0
    
    
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'xticklen', -0.02
    tplot_options, 'yticklen', -0.01
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 6
    

;---Read asf, and mlat/mlon/midn.
    ; based on linear fit for all sites.
    ; cmlon = -15.5179*cmidn+72.5124
    
    if file_test(datfn) ne 0 then tplot_restore, filename=datfn
    load = 0
    foreach site, sites do if tnames(site+'_ewo') eq '' then load = 1
    if keyword_set(load_ewo) then load = 1
    if load then begin
        foreach site, sites do begin
            asf = sread_thg_asi(utr0, site, type='asf')
            asc = sread_thg_asc(0, site, type='asf', vars=['mlat','mlon','elev'])
            
            uts = asf.utsec
            nrec = n_elements(uts)
            mlats = reform(asc.(0).mlat[hgtidx,*,*])
            mlons = reform(asc.(0).mlon[hgtidx,*,*])
            elevs = asc.(0).elev
            midns = -15.5179*(uts*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.
            
            mlats = (mlats[1:*,1:*]+mlats[1:*,0:255]+mlats[0:255,1:*]+mlats[0:255,0:255])*0.25
            mlons = (mlons[1:*,1:*]+mlons[1:*,0:255]+mlons[0:255,1:*]+mlons[0:255,0:255])*0.25
            
            tmp = minmax(mlons) & tmp = [tmp-min(midns),tmp-max(midns)]/15
            mltrng = minmax(tmp) & mltrng = mltrng-(mltrng mod dmlt)+[0,dmlt]
            mlts = smkarthm(mltrng[0], mltrng[1], dmlt, 'dx')
            nmlt = n_elements(mlts)
            
            ewos = dblarr(nrec,nmlt)
            
            
            for i=0, nrec-1 do begin
                img = double(reform(asf.img[i,*,*]))
                ; remove edge.
                edge = where(elevs lt minelev or ~finite(elevs))
                img[edge] = mean(img[0:10,0:10], /nan)
                ; crude corner average subtraction. use thg_asf_site_offset?
                img = img - img[0] > 0
                ; scale luminosity, adopted from thm_asi_merge_mosaic.
                img *= 64d/(median(img) > 1)
                tmlts = (mlons-midns[i])/15
                
                ; map to ewogram.
                for j=0, nmlt-1 do begin
                    idx = where(mlats ge mlatrng[0] and mlats le mlatrng[1] $
                        and tmlts ge mlts[j]-dmlt*0.5 and tmlts lt mlts[j]+dmlt*0.5, cnt, complement=idx2)
                    if cnt ne 0 then ewos[i,j] = mean(img[idx],/nan)
                    ;            sgopen, 0, xsize=5, ysize=5, /inch
                    ;            tmp = img & tmp[idx2] = 0
                    ;            sgtv, bytscl(tmp), ct=1, position=[0,0,1,1]
                    ;            stop
                endfor
            endfor
            
            tvar = site[0]+'_ewo'
            idx = where(ewos eq 0, cnt)
            if cnt ne 0 then ewos[idx] = !values.d_nan
            store_data, tvar, uts, ewos, mlts, limits={ytitle:'(hr)', spec:1, yrange:mltrng, no_interp:1}
        endforeach
        
        tplot_save, sites+'_ewo', filename=datfn
    endif
    
    
    
    mltrng = []
    foreach site, sites do begin
        get_data, site+'_ewo', uts, dat, val
        mltrng = [mltrng,val[where(total(dat,1,/nan) ne 0)]]
    endforeach
    nrec = n_elements(uts)
    mltrng = minmax(mltrng)
    mlts = smkarthm(mltrng[0], mltrng[1], dmlt, 'dx')
    nmlt = n_elements(mlts)
    ewos = dblarr(nrec,nmlt)
    cnts = dblarr(nrec,nmlt)
    
    
;---load atha for low latitude and pina for high latitude.
    datfn = rootdir+'/ewo_'+time_string(utr0[0],tformat='YYYY_MMDD')+'_62-63deg.tplot'
    tplot_restore, filename=datfn
    get_data, 'atha_ewo', uts, atha, val1
    
    datfn = rootdir+'/ewo_'+time_string(utr0[0],tformat='YYYY_MMDD')+'_65-66deg.tplot'
    tplot_restore, filename=datfn
    get_data, 'pina_ewo', uts, pina, val2
    
    store_data, 'atha_ewo', uts, atha, val1 & atha = 0
    store_data, 'pina_ewo', uts, pina, val2 & pina = 0

    siteweight=[1,1]    ; adjust site luminosity.
    foreach site, sites, j do begin
        get_data, site+'_ewo', uts, dat, val
        for i=0, n_elements(val)-1 do begin
            idx = where(val[i] eq mlts, cnt)
            if cnt eq 0 then continue
            tmp = where(finite(dat[*,i]), cnt)
            if cnt eq 0 then continue
            ewos[tmp,idx[0]] += dat[tmp,i]*siteweight[j]
            cnts[tmp,idx[0]] += 1
        endfor
    endforeach
    idx = where(cnts ne 0, cnt)
    if cnt ne 0 then ewos[idx] /= cnts[idx]
    idx = where(ewos eq 0, cnt)
    if cnt ne 0 then ewos[idx] = !values.d_nan
    store_data, 'ewo', uts, ewos, mlts, limits={ytitle:'(hr)', spec:1, yrange:mltrng, no_interp:1, ystyle:1}
    

;---Load other data.
    _2014_0827_load_data
    
    
    pres = 'th'+['d','e','a']+'_'
    foreach pre0, pres do begin
        get_data, pre0+'b_gsm', uts, dat
        store_data, pre0+'bmag', uts, snorm(dat)
    endforeach
    
    stop
    
    tvar = 'thd_b_gsm'
    bidx = 2
    comp = 'z'
;    bidx = 0
;    comp = 'x'
    stplot_index, tvar, bidx, newname=tvar+comp
    tvar = 'thd_b_gsm'+comp
    if comp eq 'x' then options, tvar, 'yrange', [40,130] else options, tvar, 'yrange', [5,55]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5

    tvar = 'the_b_gsm'
    stplot_index, tvar, bidx, newname=tvar+comp
    tvar = 'the_b_gsm'+comp
    if comp eq 'x' then options, tvar, 'yrange', [40,110] else options, tvar, 'yrange', [-15,35]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5

    tvar = 'tha_b_gsm'
    stplot_index, tvar, bidx, newname=tvar+comp
    tvar = 'tha_b_gsm'+comp
    if comp eq 'x' then options, tvar, 'yrange', [30,80] else options, tvar, 'yrange', [-15,25]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5

    tvar = 'rba_b_gsm'
    stplot_index, tvar, bidx, newname=tvar+comp
    tvar = 'rba_b_gsm'+comp
    if comp eq 'x' then options, tvar, 'yrange', [0,120] else options, tvar, 'yrange', [20,90]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5


    vars = ['thd_b_gsm','the_b_gsm','tha_b_gsm','rba_b_gsm']+comp
    options, vars, 'ytitle', '(nT)'
    options, vars, 'labels', ['']
    options, vars, 'ystyle', 1
    
    
    
    
    
    figfn = shomedir()+'/fig_ewo_b'+comp+'.pdf'
    ;figfn = 0
    sgopen, figfn, xsiz=6, ysize=6, /inch

    
    device, decomposed=0
    loadct2, 43
    
    pos1 = [0.2,0.7,0.85,0.95]
    pos2 = [0.2,0.1,0.85,0.65]
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    if size(figfn,/type) eq 7 then hsize = 100 else hsize = 5
    if size(figfn,/type) eq 7 then thick = 4 else thick = 2

    yrng = [-0.5,5]
    
    
    tvar = 'ewo'
    options, tvar, 'yrange', yrng
    options, tvar, 'zrange', [0,350]
    options, tvar, 'ztitle', 'Photon Count'
    options, tvar, 'ytitle', 'MLT (hr)'
    options, tvar, 'ystyle', 1
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [0,2,4]
    options, tvar, 'yminor', 5
    options, tvar, 'constant', -100
    

    tplot, tvar, trange=utr0, position=pos1, /nouttick
    
    
    dfuts = time_double('2014-08-27/'+['09:28','09:36','09:42','09:47'])
    fmlts = [0.11d,1.13,1.90,4.64]
    fdmlts = [0.01d,0.09,0.12,0.20]
    fmlats = [63.4d,66.6,68.4,64.9]
    fdmlats = [2.3d,3.4,3.6,1.2]
    sccolors = [6,4,2,1]
    
    plot, utr0, yrng, /nodata, /noerase, xstyle=5, ystyle=5, position=pos1
    for i=0,3 do plots, dfuts[i], fmlts[i], /data, psym=6, color=sccolors[i], thick=thick
    
    scposs = dblarr(4,2)
    for i=0,3 do begin
        tmp = convert_coord(dfuts[i],fmlts[i], /data, /to_normal)
        scposs[i,*] = tmp[0:1]
    endfor

    
    xyouts, pos1[0]-xchsz*12, pos1[3]-ychsz*0.8, 'a. Ewogram', /normal
    xyouts, pos1[0]+xchsz*1, pos1[3]-ychsz*1.2, /normal, '1. PINA!C 65 deg', alignment=0
    xyouts, pos1[2]-xchsz*1, pos1[1]+ychsz*1.4, /normal, '2. ATHA!C62 deg', alignment=1
    plots, utr0, [1.6,3], /data, thick=thick*2, color=0

    
    xyouts, pos2[0]-xchsz*12, pos2[3]-ychsz*0.8, 'b. B'+comp+' GSM', /normal

    vars = ['thd_b_gsm','the_b_gsm','tha_b_gsm','rba_b_gsm']+comp
    labs = ['THD','THE','THA','RBA']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, lmargin=15, rmargin=10, position=pos2)
    tplot, vars, trange=utr0, position=poss, /novtitle, /noerase
    for i=0, nvar-1 do xyouts, poss[2,i]+xchsz*1, (poss[3,i]+poss[1,i])*0.5-ychsz*0.3, labs[i], /normal, color=sccolors[i]
    for i=0, nvar-1 do begin
        ;tx = poss[0,i]+(poss[2,i]-poss[0,i])*(dfuts[i]-utr0[0])/(utr0[1]-utr0[0])
        ;plots, tx+[0,0], poss[[1,3],i], color=6, /normal
        get_data, vars[i], uts, dat, limits=lim
        tdat = interpol(dat,uts,dfuts[i])
        ty = poss[1,i]+(poss[3,i]-poss[1,i])*(tdat-lim.yrange[0])/(lim.yrange[1]-lim.yrange[0])
        ;plots, scposs[i,0]+[0,0], [ty, scposs[i,1]-ychsz*0.3], color=sccolors[i], /normal, linestyle=1
        arrow, scposs[i,0], scposs[i,1]-ychsz*0.3, scposs[i,0], ty+ychsz*0.5, color=sccolors[i], /normal, $
            /solid, hsize=hsize
    endfor

    
    
    
    sgclose
    
    stop

end
