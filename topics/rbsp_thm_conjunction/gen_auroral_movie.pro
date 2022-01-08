;+
; generate auroral movie with s/c on
;-

    rad = !dpi/180d
    deg = 180d/!dpi
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    
    
    mosmincnt = 20d
    mosmaxcnt = 200d
    top = 254
    white = 255
    modelcolors = [6,5,4,2]
    scsym = [1,4,5,6]
    
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    
    models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    
    rootdir = shomedir()+'/Google Drive/works/works/rbsp_thm_conjunction'
    
    utr0 = time_double(['2014-08-27/09:00','2014-08-27/10:30'])
    ;utr0 = ['2014-08-27/09:34:00','2014-08-27/09:34:03']
    sites = ['mcgr','whit','atha','gill','pina']
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
    
    
    utr0 = minmax(uts)
    probes = ['rba','tha','thd','the']
    nprobe = n_elements(probes)
    
    
    if tnames(probes[0]+'_fpt_mlon') eq '' then begin
        foreach tmodel, models do sgeopack_par, utr0, tmodel
        dir = -1
        
        foreach tprobe, probes do begin
            pre0 = tprobe+'_'
            
            ; load spacecraft position.
            ; tuts, rgsm.
            if tprobe eq 'rba' then begin
                spice = sread_rbsp_spice_product(utr0, probe='a')
                tuts = spice.ut_pos
                rgsm = spice.pos_gsm
                store_data, pre0+'pos_gsm', tuts, rgsm, limits=$
                    {ytitle:'(Re)', colors:rgb, labels:'GSM '+xyz}
            endif else begin
                prb = strmid(tprobe,2)
                dat = sread_thm_orbit(utr0, probes=prb)
                tuts = sfmepoch(dat.epoch,'unix')
                idx = where(tuts ge utr0[0] and tuts le utr0[1])
                tuts = tuts[idx]
                rgsm = dat.xyz_gsm[idx,*]
                store_data, pre0+'pos_gsm', tuts, rgsm, limits=$
                    {ytitle:'(Re)', colors:rgb, labels:'GSM '+xyz}
            endelse
            
            tnrec = n_elements(tuts)
            hfmlats = dblarr(tnrec,nmodel)
            hfmlons = dblarr(tnrec,nmodel)
            hfmlts = dblarr(tnrec,nmodel)
            
            for j=0, nmodel-1 do begin
                tmodel = models[j]
                t89 = (tmodel eq 't89')? 1: 0
                t96 = (tmodel eq 't96')? 1: 0
                t01 = (tmodel eq 't01')? 1: 0
                t04s = (tmodel eq 't04s')? 1: 0
                storm = (tmodel eq 't04s')? 1: 0
                
                get_data, pre0+'pos_gsm', tuts, posgsm
                tets = stoepoch(tuts,'unix')
                get_data, tmodel+'_par', tmp, pars
                pars = sinterpol(pars, tmp, tuts)
                tnrec = n_elements(tuts)
                
                for i=0, tnrec-1 do begin
                    tet = tets[i]
                    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
                    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
                    
                    xp = posgsm[i,0]
                    yp = posgsm[i,1]
                    zp = posgsm[i,2]
                    ; map to the ionosphere.
                    ; [xyz]f is the footpoint position.
                    par = reform(pars[i,*])
                    geopack_trace, xp,yp,zp, dir, par, xf,yf,zf, r0=r0, $
                        /refine, /ionosphere, $
                        t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm
                    ; convert footpoint from gsm to mag.
                    geopack_conv_coord, xf,yf,zf, /from_gsm, $
                        txf,tyf,tzf, /to_mag
                    hfmlats[i,j] = asin(tzf/r0)*deg
                    hfmlons[i,j] = atan(tyf,txf)*deg
                    hfmlts[i,j] = slon2lt(hfmlons[i,j], tet, /mag, /deg)/15
                endfor
            endfor
            store_data, pre0+'fpt_mlat', tuts, hfmlats
            store_data, pre0+'fpt_mlon', tuts, hfmlons
            store_data, pre0+'fpt_mlt', tuts, hfmlts
        endforeach
    endif

    
    
    stop
    
    
    ; plot all images in a folder.
    tpos = [0d,0,1,1]
    mosmltrng = [-90,90]
    moslatrng = [minlat,90]
    moslatticks = [55,65,75,85]>minlat

    odir = rootdir+'/movie'
    if file_test(odir,/directory) eq 0 then file_mkdir, odir
    
    
    foreach tprobe, probes do begin
        get_data, tprobe+'_fpt_mlat', tuts, mlat
        get_data, tprobe+'_fpt_mlon', tuts, mlon
        mlat = reform(sinterpol(mlat, tuts, uts))
        mlon = reform(sinterpol(mlon, tuts, uts))
        store_data, tprobe+'_fpt_mlat2', uts, mlat
        store_data, tprobe+'_fpt_mlon2', uts, mlon
    endforeach
    
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
        loadct, 1
        sgtv, timg, position=tpos
        
        loadct2, 43
        sgset_map, position=tpos, color=white, xrange=mosmltrng, $
            yrange=moslatrng, ytickv=moslatticks
        xyouts, tpos[0]+xchsz, tpos[1]+ychsz, /normal, $
            time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color=white
        
        foreach tprobe, probes, k do begin
            get_data, tprobe+'_fpt_mlat2', uts, mlat
            get_data, tprobe+'_fpt_mlon2', uts, mlon
            for j=0, nmodel-1 do $
                plots, mlon[i,j]-midn[i], mlat[i,j], psym=scsym[k], color=modelcolors[j]
            for j=0, nmodel-1 do $
                xyouts, tpos[0]+xchsz*3, tpos[3]-ychsz*(j+1)*1.5, /normal, strupcase(models[j]), color=modelcolors[j]
            for j=0, nprobe-1 do begin
                xyouts, tpos[0]+xchsz*10, tpos[3]-ychsz*(j+1)*1.5, /normal, strupcase(probes[j]), color=white
                plots, tpos[0]+xchsz*14, tpos[3]-ychsz*(j+0.8)*1.5, /normal, psym=scsym[j], color=white
            endfor
        endforeach
        
        sgclose
    endforeach

stop
    moviefn = rootdir+'/thg_asi_'+stodate(utr0[0],'%Y_%m%d')+'.mp4'
    ext = 'png'
    spic2movie, odir, moviefn, ext

end
