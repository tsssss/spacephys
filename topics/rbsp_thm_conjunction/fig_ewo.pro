;+
; Show EWO at different latitude.
;-

;---Constant.
    secofday = 86400d
    top = 254
    

;---Settings.
    utr0 = time_double(['2014-08-27/08:40','2014-08-27/10:10'])
    utr1 = time_double(['2014-08-27/08:50','2014-08-27/10:10'])
    sites = ['atha','pina']
    rootdir = shomedir()+'/Google Drive/works/works/rbsp_thm_conjunction'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir

    minelev = 5d
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
    _2014_0827_load_data
    
    
    mlats = smkarthm(62,65,1,'dx')
    sites = ['atha','pina']
    siteweights = [2,1]
    tutr = utr1 ; common time when all sites have data.
    
    foreach tmlat, mlats do begin
        
        suf0 = '_'+string(tmlat,format='(I2)')

        ; get common y-axis.
        mltrng = []
        foreach tsite, sites do begin
            get_data, tsite+'_ewo', uts, dat, val
            mltrng = [mltrng,val[where(total(dat,1,/nan) ne 0)]]
        endforeach
        idx = where(uts ge tutr[0] and uts le tutr[1])
        ut0s = uts[idx]
        nrec = n_elements(ut0s)
        mltrng = minmax(mltrng)
        mlts = smkarthm(mltrng[0], mltrng[1], dmlt, 'dx')
        nmlt = n_elements(mlts)
        ewos = dblarr(nrec,nmlt)
        cnts = dblarr(nrec,nmlt)

        ; fill total ewo with each site.
        foreach tsite, sites, j do begin
            get_data, tsite+'_ewo'+suf0, uts, dat, val
            idx = where(uts ge tutr[0] and uts le tutr[1])
            uts = uts[idx]
            dat = dat[idx,*]
            for i=0, n_elements(val)-1 do begin
                idx = where(val[i] eq mlts, cnt)
                if cnt eq 0 then continue
                tmp = where(finite(dat[*,i]), cnt)
                if cnt eq 0 then continue
                ewos[tmp,idx[0]] += dat[tmp,i]*siteweights[j]
                cnts[tmp,idx[0]] += 1
            endfor
        endforeach
        
        idx = where(cnts ne 0, cnt)
        if cnt ne 0 then ewos[idx] /= cnts[idx]
        idx = where(ewos eq 0, cnt)
        if cnt ne 0 then ewos[idx] = !values.d_nan
        store_data, 'ewo'+suf0, ut0s, ewos, mlts, limits=$
            {ytitle:'(hr)', spec:1, yrange:mltrng, no_interp:1, ystyle:1}
    endforeach
    
    
    
    
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
