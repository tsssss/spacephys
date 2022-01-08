
pro neupert_survey, utr

    tstr = time_string(utr[0],tformat='YYYY_MMDD')
    rootdir = shomedir()+'/neupert_image'
    
    ; load ae.
    omni = sread_omni(utr, vars = ['Epoch','AE_INDEX'])
    store_data, 'ae', sfmepoch(omni.epoch, 'unix'), omni.ae_index, $
        limits = {ytitle:'AE (nT)', constant:500, yrange:[0,1500]}

    ; load image aurora.
    minv = 500
    maxv = 5000
    top = 254
    ct = 39
    picpos = [0,0,1,1]
    xr = [-90,90]   ; longitude range.
    xr = [-60,60]   ; longitude range.
    yr = [50,80]    ; latitude range.
    xticks = [-60,0,60]
    yticks = [60,70,80]
    ytickpos = -60  ; latitude tick position.
    xtickname = ['20','00','04']
    white = sgcolor('white')
    xsz = 5
    ysz = 2.5
    
    deg = 180d/!dpi
    rad = !dpi/180d
    
    ; **** IMAGE/WIC.
    wic = sread_image_fuv(utr, /half)
    nrec = n_elements(wic.epoch)
    cnts = dblarr(nrec)
    
    timg = reform(wic.mltimg[0,*,*])
    tsz = size(timg,/dimensions)
    txs = indgen(tsz[0]) # (intarr(tsz[1])+1) & txs = txs-tsz[0]/2
    tys = (intarr(tsz[0])+1) # indgen(tsz[1]) & tys = tys-tsz[1]
    
    trs = sqrt(txs^2+tys^2)  & trs = 90-trs*abs(90-50)/tsz[1]
    tas = atan(tys, txs)*deg & tas = tas+90
    
    for i = 0, nrec-1 do begin
        timg = reform(wic.mltimg[i,*,*])
        
        ; remove points out of xyrange.                
        idx = where(trs lt yr[0] or trs gt yr[1] or tas lt xr[0] or tas gt xr[1], cnt)
        timg[idx] = 0
        cnts[i] = total(timg)

        timg = bytscl(timg, min = minv, max = maxv, /nan)<top
        idx = where(timg ne 0, cnt)
        if cnt le 0.1*tsz[0]*tsz[1] then cnts[i] = !values.d_nan
        
        ofn = 1
        ofn = rootdir+'/'+tstr+'/image_wic/image_wic_'+$
            sfmepoch(wic.epoch[i],'YYYY_MMDD_hhmm_ss')+'.png'
        sgopen, ofn, xsize = xsz, ysize = ysz, /inch
        sgtv, timg, position = picpos, ct = ct
        sgset_map, xrange = xr, yrange = [50,90], pos = picpos, $
            ytickpos = ytickpos, xtickname = xtickname, color = white, $
            xtickv = xticks, ytickv = yticks
        sgclose
    endfor
    
    store_data, 'img_wic_cnts', sfmepoch(wic.epoch,'unix'), cnts, $
        limits = {ytitle: 'IMAGE/WIC!CCounts!Cin FOV'}
    
    ; **** Polar/UVI.
    uvi = sread_polar_uvi(utr, /half)
    nrec = n_elements(uvi.epoch)
    cnts = dblarr(nrec)
    
    timg = reform(uvi.mltimg[0,*,*])
    tsz = size(timg,/dimensions)
    txs = indgen(tsz[0]) # (intarr(tsz[1])+1) & txs = txs-tsz[0]/2
    tys = (intarr(tsz[0])+1) # indgen(tsz[1]) & tys = tys-tsz[1]
    
    trs = sqrt(txs^2+tys^2)  & trs = 90-trs*abs(90-50)/tsz[1]
    tas = atan(tys, txs)*deg & tas = tas+90
    
    for i = 0, nrec-1 do begin
        timg = reform(uvi.mltimg[i,*,*])
        
        ; remove points out of xyrange.
        idx = where(trs lt yr[0] or trs gt yr[1] or tas lt xr[0] or tas gt xr[1], cnt)
        timg[idx] = 0
        cnts[i] = total(timg)
        
        timg = timg<top
        idx = where(timg ne 0, cnt)
        if cnt le 0.1*tsz[0]*tsz[1] then cnts[i] = !values.d_nan
        
        ofn = 1
        ofn = rootdir+'/'+tstr+'/polar_uvi/polar_uvi_'+$
            sfmepoch(uvi.epoch[i],'YYYY_MMDD_hhmm_ss')+'.png'
        sgopen, ofn, xsize = xsz, ysize = ysz, /inch
        sgtv, timg, position = picpos, ct = ct
        sgset_map, xrange = xr, yrange = [50,90], pos = picpos, $
            ytickpos = ytickpos, xtickname = xtickname, color = white, $
            xtickv = xticks, ytickv = yticks
        sgclose
    endfor
    
    store_data, 'po_uvi_cnts', sfmepoch(uvi.epoch,'unix'), cnts, $
        limits = {ytitle: 'Polar/UVI!CCounts!Cin FOV'}
    
    
    ; **** polar hydra.
    hyd = sread_polar_hydra(utr)
    uts = sfmepoch(hyd.epoch, 'unix')
    store_data, 'po_hyd_ele', uts, hyd.jee, hyd.ene, limits = $
        {ytitle:'Polar!CHydra!CElectron!Energy (eV)', spec:1, $
        ylog:1, yrange:hyd.ene[[0,-1]], zlog:1, zrange:[1e4,5e8], $
        ztitle:'1/(cm!U2!N-s-sr)', ystyle:1}
    store_data, 'po_hyd_ion', uts, hyd.jei, hyd.eni, limits = $
        {ytitle:'Polar!CHydra!CIon!Energy (eV)', spec:1, $
        ylog:1, yrange:hyd.eni[[0,-1]], zlog:1, zrange:[1e4,1e8], $
        ztitle:'1/(cm!U2!N-s-sr)', ystyle:1}
    
    ofn = rootdir+'/fig/neupert_cnt_'+tstr+'.pdf'
    sgopen, ofn, xsize = 6, ysize = 8, /inch
    device, decomposed = 0
    loadct2, 43
    vars = ['ae','img_wic_cnts','po_uvi_cnts','po_hyd_ele','po_hyd_ion']
    tplot, vars, trange = utr
    sgclose
    
    ofn = rootdir+'/data/neupert_cnt_'+tstr+'.tplot'
    tplot_save, vars, filename = ofn
    
end

utrs = []
;utr = time_double(['2001-04-23/01:00','2001-04-23/04:00'])  ; no image at the start of substorm.
;utrs = [[utrs],[utr]]

utr = time_double(['2001-10-21/16:00','2001-10-21/19:00'])
utrs = [[utrs],[utr]]
utr = time_double(['2001-10-22/10:00','2001-10-22/13:00'])
utrs = [[utrs],[utr]]
utr = time_double(['2002-10-27/13:00','2002-10-27/16:00'])
utrs = [[utrs],[utr]]
utr = time_double(['2002-10-29/09:00','2002-10-29/12:00'])
utrs = [[utrs],[utr]]
utr = time_double(['2002-11-03/15:30','2002-11-03/18:30'])
utrs = [[utrs],[utr]]
;utr = time_double(['2003-05-23/02:30','2003-05-23/05:30'])  ; air glow polution.
;utr = time_double(['2003-07-16/06:00','2003-07-16/09:00'])  ; need to fix spinphase.

nutr = n_elements(utrs)/2
for i = 0, nutr-1 do neupert_survey, utrs[*,i]
;for i = 0, 0 do neupert_survey, utrs[*,i]
end