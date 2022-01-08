;+
; generate auroral movie to be used in camshift.
;-

    rad = !dpi/180d
    deg = 180d/!dpi
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    
    
    mosmincnt = 10d
    mosmaxcnt = 400d
    top = 254
    white = 255
    ct = 1
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    
    models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    
    rootdir = shomedir()+'/Google Drive/works/works/rbsp_thm_conjunction'
    
    utr0 = time_double(['2014-12-22/05:16','2014-12-22/05:36'])
    utr0 = time_double(['2014-12-22/05:26','2014-12-22/05:36'])
    sites = ['whit','chbg','fykn','nrsq']
    

    utr0 = time_double(['2014-12-22/02:00','2014-12-22/02:30'])
    ;utr0 = time_double(['2014-12-22/01:30','2014-12-22/02:00'])
    sites = ['fykn','whit','snkq','chbg','nrsq']
    sites = ['fykn','snkq','chbg','nrsq']
    
    
    utr0 = time_double(['2015-04-17/03:00','2015-04-17/03:30'])
    sites = ['snkq','kuuj','pina']
    
        
    type = 'asf'
    
    asifn = rootdir+'/thg_asf_mosaic_'+stodate(utr0[0],'%Y_%m%d_%H%M')+'.cdf'
    if file_test(asifn) eq 0 then begin
        tfn = sread_thg_mosaic(utr0, sites, type=type, /weight, notop=1)
        file_copy, tfn, asifn
        file_delete, tfn
    endif
    
    
    if tnames('asf_mos') eq '' then begin
        cdfs = scdfread(asifn)

        uts = (*cdfs[0].value)
        mos = (*cdfs[1].value)
        midn= (*cdfs[2].value)
        mlt = (*cdfs[3].value)
        mlat= (*cdfs[4].value)
        imgsz = (*cdfs[5].value)
        pxidx = (*cdfs[6].value)
        minlat = (*cdfs[7].value)[0]

        img = bytarr(imgsz)
        nrec = n_elements(uts)

        txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
        tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
        txs = txs/imgsz[0]*2
        tys = tys/imgsz[0]*2
        mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
        mlts = atan(tys, txs)*deg+90    ; in deg.

        store_data, 'asf_mos', uts, mos, pxidx
        store_data, 'asf_info', 0, {imgsz:imgsz, minlat:minlat, mlats:mlats, mlts:mlts}
    endif
    
    get_data, 'asf_mos', uts, mos, pxidx
    get_data, 'asf_info', 0, asfinfo
    
    
    ; plot all images in a folder.
    tpos = [0d,0,1,1]
    mosmltrng = [-1,1]*90
    moslatrng = [minlat,90]
    moslatticks = [55,65,75,85]>minlat

    odir = rootdir+'/movie_'+time_string(uts[0],tformat='YYYY_MMDD_hhmm')
    if file_test(odir,/directory) eq 0 then file_mkdir, odir
    
    
    foreach tut, uts, i do begin
        ofn = odir+'/thg_asi_'+stodate(tut,'%Y_%m%d_%H%M_%S')+'.png'
        ;ofn = 0
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        timg = bytscl(timg, min=mosmincnt, max=mosmaxcnt, top=top)
        
        sgopen, ofn, xsize=10, ysize=5, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        device, decomposed = 0
        loadct, ct
        sgtv, timg, position=tpos
        
        loadct2, 43
        sgset_map, position=tpos, color=white, xrange=mosmltrng, $
            yrange=moslatrng, ytickv=moslatticks
        xyouts, tpos[0]+xchsz, tpos[1]+ychsz, /normal, $
            time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color=white

        sgclose
    endforeach


    moviefn = rootdir+'/thg_asi_'+stodate(utr0[0],'%Y_%m%d_%H%M')+'.mp4'
;    ;fns = strarr(n_elements(uts))
;    ;for i=0, n_elements(fns)-1 do fns[i] = 'thg_asi_'+stodate(uts[i],'%Y_%m%d_%H%M_%S')+'.png'
    ext = 'png'
    spic2movie, odir, moviefn, ext;, filenames=odir+fns

end
