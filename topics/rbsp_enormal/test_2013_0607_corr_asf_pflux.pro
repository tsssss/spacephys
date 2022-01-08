;+
; correlation map between each pixel and pflux.
;-


;; constants.
;deg = 180d/!dpi
;rad = !dpi/180
;re = 6378d & re1 = 1d/re
;r0 = 100d/re+1
;dr0 = 3
;
;top = 254
;maxcnt = 800
;maxcnt = 500
;mincnt = 50
;minlat = 55
;
;; range of correlation for bytscl.
;mincor0 = 0.2
;maxcor0 = 0.7
;
;xr = [-90,0]
;;xr = [-90,0]
;yr = [minlat,90]
;yticks = [55,65,75]
;
;
;tutr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
;
;probes = ['a','b']
;
;
;
;; **** load asf data.
;_2013_0607_load_data
;
;tplot_options, 'ygridstyle', 0
;tplot_options, 'yticklen', 0.01
;tplot_options, 'num_lab_min', 5
;tplot_options, 'ystyle', 1
;
;
;
;get_data, 'asf_mos', uts, mos, pxidx
;get_data, 'asf_info', tmp, info
;imgsz = info.imgsz
;minlat = info.minlat
;
;
;; mlt, mlats at each pixel.
;txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
;tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
;txs = txs/imgsz[0]*2
;tys = tys/imgsz[0]*2
;mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
;mlts = atan(tys, txs)*deg+90    ; in deg, 0 deg at midnight.
;mlts = mlts/15          ; in hr.
;idx = where(mlts gt 12, cnt)    ; make mlt in [-12,12]
;if cnt ne 0 then mlts[idx] = mlts[idx]-12
;
;
;
;
;idx = where(uts ge tutr[0] and uts le tutr[1], idx)
;uts = uts[idx]
;mos = mos[idx,*]
;ets = stoepoch(uts,'unix')
;nrec = n_elements(uts)
;
;asfs = fltarr(nrec,imgsz[0],imgsz[1])
;for i = 0, nrec-1 do begin
;    timg = fltarr(imgsz)
;    timg[pxidx] = mos[i,*]
;    asfs[i,*,*] = timg
;endfor
;
;ofn = shomedir()+'/correlation_map.pdf'
;;ofn = 0
;sgopen, ofn, xsize = 10, ysize = 5, /inch
;
;
;pos1 = [0,0,0.5,1]
;pos2 = [0.5,0,1,1]
;
;
;foreach tprobe, probes do begin
;    pre0 = 'rbsp'+tprobe+'_'
;    tvar = pre0+'pf_fac'
;    get_data, tvar, tuts, dat
;    spara = interpol(dat[*,0], tuts, uts)
;    
;    ; set the range to calculate correlation.
;    pidx = where(mlats ge 58 and mlats le 65 and mlts ge -2.5 and mlts le -1)
;    pidx = where(mlats ge 58 and mlts ge -2.5 and mlts le -1)
;
;    cormap = fltarr(imgsz)
;    for i = 0, n_elements(pidx)-1 do begin
;        tpos = array_indices(imgsz, pidx[i], /dimensions)
;        cormap[tpos[0],tpos[1]] = c_correlate(asfs[*,tpos[0],tpos[1]],spara,0)
;    endfor
;    cormap0 = cormap    ; save a copy of the original one.
;    cormap = cormap[0:imgsz[0]/2,*] ; only keep the left half.
;    
;    
;    xchsz = double(!d.x_ch_size)/!d.x_size
;    ychsz = double(!d.y_ch_size)/!d.y_size
;    
;    
;    pos = (tprobe eq 'a')? pos1: pos2
;    device, decomposed = 0
;    loadct, 1
;    sgtv, bytscl(cormap, min=mincor0, max=maxcor0, top = 254), position = pos
;    loadct2, 43
;
;    get_data, pre0+'fpt_mlt', tuts, fptmlts, limits = lim
;    get_data, pre0+'fpt_mlat', tuts, fptmlats
;    
;    idx = where(tuts ge tutr[0] and tuts le tutr[1])
;    cmlt = mean(minmax(fptmlts[idx,*])) & if cmlt ge 12 then cmlt = cmlt-24
;    dmlt = 0.1
;    pidx = where(mlts le cmlt-dmlt or mlts ge cmlt+dmlt)
;    tmp = cormap0
;    tmp[pidx] = 0
;    mcval = max(tmp, mcidx)
;    mcmlt = mlts[mcidx]+24  ; the mlt of the max correlection point, in hr.
;    mcmlat = mlats[mcidx]   ; the mlat in deg.
;    print, 'MLT, MLat:', mcmlt, mcmlat
;    print, 'Mac Correlation:', mcval
;    sgset_map, position = pos, color = 255, xrange = xr, $
;        yrange = yr, ytickv = yticks
;    plots, (mcmlt-24)*15, mcmlat, /data, psym = 1, thick = 4, color = 254
;
;    
;    get_data, pre0+'map_coef', tuts, coef
;    mapcoef = interpol(coef[*,0], tuts, uts)
;    store_data, pre0+'smap', uts, spara*mapcoef, limits = $
;        {ytitle:'(mW/m!U2!N)', labels:'S!D||!N@100 km', ylog:1, yrange:[1,200]}
;    tpos = array_indices(imgsz, mcidx, /dimensions)
;    store_data, pre0+'pcnt', uts, asfs[*,tpos[0],tpos[1]], limits = $
;        {ytitle:'Photon Count', labels:'MLT:'+sgnum2str(mlts[mcidx]+24,ndec=1)+' hr'+$
;        '!C  MLat:'+sgnum2str(mlats[mcidx],ndec=1)+' deg'}
;
;    
;    models = lim.labels
;    nmodel = n_elements(models)
;    for i = 0, nmodel-1 do begin
;        tmlats = sinterpol(fptmlats[*,i], tuts, uts, /quadratic)
;        tmlts = sinterpol(fptmlts[*,i], tuts, uts, /quadratic)
;        
;        tcolor = 6-i
;        tidx = n_elements(tmlts)/2
;        plots, tmlts[tidx]*15, tmlats[tidx], /data, color = tcolor, psym = 6, symsize=0.7
;        xyouts, pos[0]+xchsz, pos[1]+ychsz*(i+2), /normal, $
;            strupcase(models[i]), color = tcolor
;    endfor
;    xyouts, pos[0]+xchsz, pos[1]+ychsz*1, /normal, 'RBSP-'+strupcase(tprobe), color=255
;
;endforeach
;
;sgclose
;
;
;options, 'rbspa_pcnt', 'yrange', [100,400]
;options, 'rbspa_pcnt', 'yticks', 2
;options, 'rbspa_pcnt', 'yminor', 5
;options, 'rbspb_pcnt', 'yrange', [100,300]
;options, 'rbspb_pcnt', 'yticks', 2
;options, 'rbspb_pcnt', 'yminor', 5
;
;ofn = 0
;ofn = shomedir()+'/s_vs_asf_best_corr.pdf'
;sgopen, ofn, xsize = 5, ysize = 5, /inch
;
;xchsz = double(!d.x_ch_size)/!d.x_size
;ychsz = double(!d.y_ch_size)/!d.y_size
;
;vars = ['rbspa_'+['smap','pcnt'],'rbspb_'+['smap','pcnt']]
;nvar = n_elements(vars)
;poss = sgcalcpos(nvar)
;
;tplot, vars, trange = tutr, position = poss
;figlabs = ['a','b','c','d']+'.'
;for i = 0, nvar-1 do xyouts, poss[0,i]-xchsz*8, poss[3,i]-ychsz*1, /normal, figlabs[i]
;
;sgclose
;
;


vars = ['rbspa_'+['smap','pcnt'],'rbspb_'+['smap','pcnt']]
;foreach tvar, vars do stplot_mor, tvar


ofn = shomedir()+'/s_vs_asf_psd.pdf'
; ofn = 0
sgopen, ofn, xsize=5, ysize=6, /inch

poss = sgcalcpos(2, ypad=8)
xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


for i = 0, 1 do begin
    
    tprobe = (i eq 0)? 'a': 'b'
    
    get_data, 'rbsp'+tprobe+'_smap_mor_fft_info', tmp, fftinfo0
    get_data, 'rbsp'+tprobe+'_pcnt_mor_fft_info', tmp, fftinfo1
    
    tpos = poss[*,i]
    xr = [1,1e3]
    yr = [1e-2,1e2]
    title = (i eq 0)? 'RBSP-A': 'RBSP-B'
    
    plot, xr, yr, /nodata, /noerase, position=tpos, $
        xrange=xr, xstyle=1, xlog=1, xtitle='Period (sec)', $
        yrange=yr, ystyle=1, ylog=1, ytitle='Normalized PS', $
        title=title
    tx = fftinfo0.ps
    ty = fftinfo0.gws/fftinfo0.sigma2
    oplot, tx, ty, color=sgcolor('blue')
    tx = fftinfo1.ps
    ty = fftinfo1.gws/fftinfo1.sigma2
    oplot, tx, ty, color=sgcolor('red')
    
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1
    plots, tx+[0,2]*xchsz, ty+[0,0], /normal, color=sgcolor('blue')
    xyouts, tx+3*xchsz, ty-ychsz*0.2, /normal, 'Poynting flux'
    
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*2
    plots, tx+[0,2]*xchsz, ty+[0,0], /normal, color=sgcolor('red')
    xyouts, tx+3*xchsz, ty-ychsz*0.2, /normal, 'Aurora Brightness'
    
endfor

sgclose


end
