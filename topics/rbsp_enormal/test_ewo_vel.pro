;+
; Show EWO at specific latitude range.
;-

;---Constant.
    secofday = 86400d
    top = 254
    

;---Settings.
    utr0 = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    utr1 = time_double(['2013-06-07/04:54','2013-06-07/04:57'])
    tsite = 'pina'
;    rootdir = shomedir()+'/Google Drive/works/works/xxx'
    rootdir = shomedir()
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir

    minelev = 10d
    height = 110d       ; km.
    hgtidx = where([90,110,150] eq height)
    dmlt = 0.05         ; hr. resolution in mlt.

    
    
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'xticklen', -0.02
    tplot_options, 'yticklen', -0.01
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 6
    

;---Load other data.
    
    mlatrng = [64d,64.5]    ; mlat range.
    mlatrng = [63.5d,64]    ; mlat range.
    ;mlatrng = [64d,65]    ; mlat range.
    tutr = utr1 ; common time when all sites have data.
    
    suf0 = '_'+string(mlatrng[0],format='(I2)')

    asf = sread_thg_asi(tutr, tsite, type='asf')
    asc = sread_thg_asc(0, tsite, type='asf', vars=['mlat','mlon','elev'])
    
    uts = asf.utsec
    nrec = n_elements(uts)
    mlats = reform(asc.(0).mlat[hgtidx,*,*])
    mlons = reform(asc.(0).mlon[hgtidx,*,*])
    elevs = asc.(0).elev
    midns = -15.5179*(uts*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.
    
    ; converts corner position to center position.
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
        endfor
    endfor
    
    tvar = tsite+'_ewo'+suf0
    idx = where(ewos eq 0, cnt, complement=idx2)
    if cnt ne 0 then ewos[idx] = !values.d_nan
    store_data, tvar, uts, ewos, mlts, limits={ytitle:'(hr)', spec:1, yrange:[-2.2,-1.3], no_interp:1}
    

    tplot, tvar, trange=tutr
    
    stop
    
    
    
    figfn = shomedir()+'/fig_ewo.pdf'
    figfn = 0
    sgopen, figfn, xsiz=12, ysize=12, /inch

    
    device, decomposed=0
    loadct2, 43
    
    pos1 = [0.2,0.1,0.85,0.95]
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size


    yrng = [-0.5,5]
    
    
    lab0s = ['a','b','c','d']+'. '
    labs = string(reverse(mlats),format='(I2)')
    suf1 = '_'+labs
    vars = 'ewo'+suf1
    nvar = n_elements(vars)

    options, vars, 'yrange', yrng
    options, vars, 'ztitle', 'Photon Count'
    options, vars, 'ytitle', 'MLT (hr)'
    options, vars, 'ystyle', 1
    options, vars, 'yticks', 2
    options, vars, 'ytickv', [0,2,4]
    options, vars, 'yminor', 5
    options, vars, 'constant', -100
    
    zrngs = [[1,500],[1,700],[1,900],[1,1000]]
    ;zrngs = [[1,500],[1,500],[1,500],[1,500]]
    for i=0, nvar-1 do options, vars[i], 'zrange', zrngs[*,i]
    
    poss = sgcalcpos(nvar, position=pos1)
    tplot, vars, position=poss, /novtitle
    for i=0, nvar-1 do begin
        xyouts, poss[0,i]-xchsz*12, poss[3,i]-ychsz*0.8, /normal, lab0s[i]+labs[i]+' deg'
        xyouts, poss[0,i]+xchsz*1, poss[3,i]-ychsz*1.2, /normal, '1. PINA', alignment=0
        xyouts, poss[2,i]-xchsz*1, poss[1,i]+ychsz*0.4, /normal, '2. ATHA x '+sgnum2str(siteweights[0]), alignment=1
    endfor
    
    sgclose
    
end
