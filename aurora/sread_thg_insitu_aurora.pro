;---constant.
    rad = !dpi/180d
    deg = 180d/!dpi
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    
;---parameters for 2015-04-16.
    utr0 = time_double(['2015-04-16/07:45','2015-04-16/08:20'])
    sites = ['mcgr','whit']
    probes = ['a','b']

    dis0 = 5.5d
    mltrng = [19,22]+0.2
    zzzrng = [-2.5,0]
    dddrng = [4.5,6.5]
    mlts = smkarthm(mltrng[0], mltrng[1], 0.1, 'dx')*15*rad
    zzzs = smkarthm(zzzrng[0], zzzrng[1], 0.4, 'dx')
    ddds = smkarthm(dddrng[0], dddrng[1], 0.2, 'dx')
    
    
    mosmincnt = 20d
    mosmaxcnt = 250d

;---parameters for 2013-05-01.
    utr0 = time_double(['2013-05-01/07:25','2013-05-01/07:50'])
    sites = ['atha','tpas']
    probes = ['b']
    
    dis0 = 5.5d
    mltrng = [-2,1]
    zzzrng = [0,2.5]
    dddrng = [4,6]
    mlts = smkarthm(mltrng[0], mltrng[1], 0.1, 'dx')*15*rad
    zzzs = smkarthm(zzzrng[0], zzzrng[1], 0.4, 'dx')
    ddds = smkarthm(dddrng[0], dddrng[1], 0.2, 'dx')


    ; settings.
    load_pos = 0    ; load rbsp pos, fpt, and grid pos, fpt.
    load_cdf = 0
    plot_mlt = 0
    to_insitu = 1
    notilt = 0
    
    modidx = 0      ; t89 is usually the best model...
    
    tet = stoepoch(utr0[0],'unix')
    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date, tilt=tilt


    vgridinfo = {yvalue:zzzs, ytitle:'GSM Z (Re)', $
        xvalue:mlts/15/rad, xtitle:'MLT (hr)', title:'Vertical plane'}
    hgridinfo = {yvalue:ddds, ytitle:'R (Re)', $
        xvalue:mlts/15/rad, xtitle:'MLT (hr)', title:'Equatorial plane'}


;---the in-situ grids of the vertical and horizontal planes.
; vgrids, hgrids.
    lgrid = n_elements(ddds)
    mgrid = n_elements(mlts)
    ngrid = n_elements(zzzs)
    tilt = tilt*rad
    
    vgrids = dblarr(mgrid,ngrid,3)
    for i=0, mgrid-1 do begin
        for j=0, ngrid-1 do begin
            vgrids[i,j,0] = -dis0*cos(mlts[i])
            vgrids[i,j,1] = -dis0*sin(mlts[i])
            vgrids[i,j,2] = zzzs[j]
            tmp = transpose(reform(vgrids[i,j,*]))
            if keyword_set(notilt) then continue
            srotate, tmp, tilt, 1
            vgrids[i,j,*] = tmp
        endfor
    endfor
    
    hgrids = dblarr(mgrid,lgrid,3)
    for i=0, mgrid-1 do begin
        for j=0, lgrid-1 do begin
            hgrids[i,j,0] = -ddds[j]*cos(mlts[i])
            hgrids[i,j,1] = -ddds[j]*sin(mlts[i])
            hgrids[i,j,2] = 0
            tmp = transpose(reform(hgrids[i,j,*]))
            if keyword_set(notilt) then continue
            srotate, tmp, tilt, 1
            hgrids[i,j,*] = tmp
        endfor
    endfor
    




;---input.
    ; sites, must provide.
    if n_elements(sites) eq 0 then message, 'no site is specified ...'
    ; utr0, time range in UT sec.
    if n_elements(utr0) ne 2 then message, 'need start and end time ...'
    ; vgrids, in [m,n,3], the x,y,z coord in gsm for a in-situ surface.
    if n_elements(vgrids) eq 0 then message, 'need vgrids in gsm ...'
    if n_elements(hgrids) eq 0 then message, 'need hgrids in gsm ...'


    ; model, default to use all models.
    if n_elements(models) eq 0 then models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    modcolors = [6,5,4,3,1]
    modcolors = modcolors[0:nmodel-1]
    
    ; root directory for holding the cdf file and aurora images and tplot data.
    datestr = time_string(utr0[0],tformat='YYYY_MMDD')
    rootdir = shomedir()+'/'+datestr
    if file_test(rootdir) eq 0 then file_mkdir, rootdir
    



;---generate the MLT specialized cdf file.
    if n_elements(minlat) eq 0 then minlat = 55
    if n_elements(type) eq 0 then type = 'asf'
    asifn = rootdir+'/thg_asf_mosaic_'+datestr+'.cdf'

    if load_cdf then begin
        tfn = sread_thg_mosaic(utr0, sites, type=type, minlat=minlat, $
            dark=0, notop=1)
        file_copy, tfn, asifn, /allow_same, /overwrite
        file_delete, tfn
        vars = ['asf_mos','asf_info']
        store_data, vars, /delete
    endif


;---generate the MLT images in ionosphere.
    mltfn = rootdir+'/thg_asf_mosaic_'+datestr+'.tplot'
    if tnames('asf_mos') eq '' then load_mlt = 1 else load_mlt = 0
    get_data, 'asf_mos', uts
    if uts[0] ne utr0[0] then load_mlt = 1

    if load_mlt then begin
        ; read mosaic images.
        tmp = scdfread(asifn)
        uts = (*tmp[0].value)
        mos = (*tmp[1].value)
        midn= (*tmp[2].value)
        mlt = (*tmp[3].value)
        mlat= (*tmp[4].value)
        imgsz = (*tmp[5].value)
        pxidx = (*tmp[6].value)
        minlat = (*tmp[7].value)[0]
        
        img = bytarr(imgsz)
        nrec = n_elements(uts)
        
        ; mosaic grid in the ionosphere.
        ixs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & ixs = ixs-imgsz[0]/2
        iys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & iys = iys-imgsz[0]/2
        ixs = ixs/imgsz[0]*2
        iys = iys/imgsz[0]*2
        imlats = 90-sqrt(ixs^2+iys^2)*(90-minlat)
        imlts = atan(iys, ixs)*deg+90    ; in deg.
        
        store_data, 'asf_mos', uts, mos, pxidx
        store_data, 'asf_info', 0, {imgsz:imgsz, minlat:minlat, $
            mlats:imlats, mlts:imlts, xxxs:ixs, yyys:iys, midn:midn}
        
        ; do not save it, since the file size is huge.
        ; tplot file is 2GB while the cdf file is only 50MB.
        ; ; update data file.
        ; vars = ['asf_mos','asf_info']
        ; tplot_save, vars, filename = mltfn
    endif
    
    
    
;---load s/c location in gsm, map grid and s/c pos to 100km using all models.
    nprobe = n_elements(probes)
    scsyms = [1,6]

    gridfn = rootdir+'/'+datestr+'_mapping_grid.tplot'
    dir = -1    ; trace to northern hemisphere.
    
    if load_pos then begin
        ; load spice kernel, model parameters.
        defsysv,'!rbsp_spice', exists=flag
        if flag eq 0 then rbsp_load_spice_kernels, trange = utr0
        get_data, 'asf_mos', uts
        foreach tmodel, models do sgeopack_par, utr0, tmodel
        
    ;---load pos and map it to the ionosphere using all models.
    ; rbspx_[pos_gsm,fpt_[mlat,mlon,mlt]].
    ; [model]_par.
        foreach tprobe, probes do begin
            pre0 = 'rbsp'+tprobe+'_'
            ; load pos in gsm.
            tuts = smkarthm(utr0[0],utr0[1],60,'dx')
            rbsp_load_spice_state, probe=tprobe, coord='gsm', times=tuts, /no_spice_load
            get_data, pre0+'state_pos_gsm', tuts, posgsm
            posgsm = posgsm*re1
            tdat = posgsm
            store_data, pre0+'pos_gsm', tuts, posgsm
            store_data, pre0+'state_*', /delete

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
            store_data, pre0+'fpt_mlat', uts, sinterpol(hfmlats,tuts,uts)
            store_data, pre0+'fpt_mlon', uts, sinterpol(hfmlons,tuts,uts)
            store_data, pre0+'fpt_mlt', uts, sinterpol(hfmlts,tuts,uts)
        endforeach



    ;---map in-situ grid to ionosphere.
    ; [model]_[vfmlat,vfmlon,hfmlat,hfmlon].
    ; insitu_[vgrid,hgrid].
        tuts = smkarthm(utr0[0], utr0[1], 3, 'n')   ; will try other ways.
        ntut = n_elements(tuts)
        tets = stoepoch(tuts, 'unix')
        foreach tmodel, models do begin
            t89 = (tmodel eq 't89')? 1: 0
            t96 = (tmodel eq 't96')? 1: 0
            t01 = (tmodel eq 't01')? 1: 0
            t04s = (tmodel eq 't04s')? 1: 0
            storm = (tmodel eq 't04s')? 1: 0
            
            ; map the vgrids for each time.
            vfmlats = dblarr(ntut,mgrid,ngrid) ; fpt mlat, in deg.
            vfmlons = dblarr(ntut,mgrid,ngrid) ; fpt mlon, in deg.

            hfmlats = dblarr(ntut,mgrid,lgrid) ; fpt mlat, in deg.
            hfmlons = dblarr(ntut,mgrid,lgrid) ; fpt mlon, in deg.
            
            get_data, tmodel+'_par', tmp, pars
            pars = sinterpol(pars, tmp, uts)

            for i=0, ntut-1 do begin
                tet = tets[i]
                geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
                geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date

                for j=0, mgrid-1 do begin
                    for k=0, ngrid-1 do begin
                        ; in-situ position.
                        xp = vgrids[j,k,0]
                        yp = vgrids[j,k,1]
                        zp = vgrids[j,k,2]
                        ; map to the ionosphere.
                        ; [xyz]f is the footpoint position.
                        par = reform(pars[i,*])
                        geopack_trace, xp,yp,zp, dir, par, xf,yf,zf, r0=r0, $
                            /refine, /ionosphere, $
                            t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm
                        ; convert footpoint from gsm to mag.
                        geopack_conv_coord, xf,yf,zf, /from_gsm, $
                            txf,tyf,tzf, /to_mag
                        vfmlats[i,j,k] = asin(tzf/r0)*deg
                        vfmlons[i,j,k] = atan(tyf,txf)*deg
                    endfor
                endfor
                
                for j=0, mgrid-1 do begin
                    for k=0, lgrid-1 do begin
                        ; in-situ position.
                        xp = hgrids[j,k,0]
                        yp = hgrids[j,k,1]
                        zp = hgrids[j,k,2]
                        ; map to the ionosphere.
                        ; [xyz]f is the footpoint position.
                        par = reform(pars[i,*])
                        geopack_trace, xp,yp,zp, dir, par, xf,yf,zf, r0=r0, $
                            /refine, /ionosphere, $
                            t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm
                        ; convert footpoint from gsm to mag.
                        geopack_conv_coord, xf,yf,zf, /from_gsm, $
                            txf,tyf,tzf, /to_mag
                        hfmlats[i,j,k] = asin(tzf/r0)*deg
                        hfmlons[i,j,k] = atan(tyf,txf)*deg
                    endfor
                endfor
            endfor
            tvar = tmodel+'_vfmlat'
            store_data, tvar, tuts, vfmlats, vgridinfo
            tvar = tmodel+'_vfmlon'
            store_data, tvar, tuts, vfmlons, vgridinfo
            
            tvar = tmodel+'_hfmlat'
            store_data, tvar, tuts, hfmlats, hgridinfo
            tvar = tmodel+'_hfmlon'
            store_data, tvar, tuts, hfmlons, hgridinfo
        endforeach
        
        store_data, 'insitu_vgrid', utr0, vgrids, vgridinfo
        store_data, 'insitu_hgrid', utr0, hgrids, hgridinfo
        vars = ['insitu_vgrid', 'insitu_hgrid', models+'_par', $
            models+'_hfmlat', models+'_hfmlon', $
            models+'_vfmlat', models+'_vfmlon', $
            'rbsp'+probes+'_fpt_mlat','rbsp'+probes+'_fpt_mlon']
        tplot_save, vars, filename=gridfn
    endif else begin
        vars = ['insitu_vgrid', 'insitu_hgrid', models+'_par', $
            models+'_hfmlat', models+'_hfmlon', $
            models+'_vfmlat', models+'_vfmlon', $
            'rbsp'+probes+'_fpt_mlat','rbsp'+probes+'_fpt_mlon']
        store_data, vars, /delete
        tplot_restore, filename=gridfn
    endelse



;---export mosaic images to png.
    mosdir = rootdir+'/asf_mlt_img'
    tpos = [0,0,1,1]
    white = 255
    top = 254
    symsz = 0.8
    
    mosxsz = 10d
    mosysz = 5d
    mosmltrng = [-90,90]
    moslatrng = [minlat,90]
    moslatticks = [55,65,75,85]
    

    if plot_mlt then begin     
        get_data, 'asf_mos', uts, mos, pxidx
        get_data, 'asf_info', 0, asfinfo
        
        nrec = n_elements(uts)
        for i=0, nrec-1 do begin
            ; prepare time and image.
            timg = fltarr(imgsz)
            timg[pxidx] = mos[i,*]
            timg = bytscl(timg, min=mosmincnt, max=mosmaxcnt, top=top)
            tut = uts[i]

            tfn = rootdir+'/thg_asf_'+$
                time_string(tut,tformat='YYYY_MMDD_hhmm_ss')+'.png'
            ;tfn = 0
            sgopen, tfn, xsize=mosxsz, ysize=mosysz, /inch

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
            
            for j=0, nprobe-1 do begin
                pre0 = 'rbsp'+probes[j]+'_'
                get_data, pre0+'fpt_mlat', tuts, scmlats
                get_data, pre0+'fpt_mlon', tuts, scmlons
                for k=0, nmodel-1 do begin
                    tmlt = scmlons[i,k]-midn[i]
                    tlat = scmlats[i,k]
                    plots, tmlt, tlat, /data, color=modcolors[k], $
                        psym=scsyms[j], symsize=symsz
                endfor
                xyouts, tpos[0]+xchsz, tpos[1]+ychsz*(2+j), /normal, $
                    'RBSP-'+strupcase(probes[j]), color=white
                plots, tpos[0]+xchsz*8, tpos[1]+ychsz*(2+0.3+j), /normal, $
                    psym=scsyms[j], color=white, symsize=symsz
            endfor
            
            for j=0, nmodel-1 do begin
                xyouts, tpos[0]+xchsz, tpos[1]+ychsz*(4+j), /normal, $
                    strupcase(models[j]), color=modcolors[j]
            endfor
            sgclose
        endfor
    endif


;---map the MLT images to in situ grids.
    if to_insitu then begin
        ; plot in-situ vgrids to aurora for test purpose.
        get_data, 'asf_mos', uts, mos, pxidx
        get_data, 'asf_info', 0, asfinfo
        
        midn = asfinfo.midn
        imgsz = asfinfo.imgsz

        get_data, models[0]+'_vfmlat', tuts, vfmlats
        get_data, models[0]+'_vfmlon', tuts, vfmlons

        get_data, models[0]+'_hfmlat', tuts, hfmlats
        get_data, models[0]+'_hfmlon', tuts, hfmlons

        timg = fltarr(imgsz)
        timg[pxidx] = mos[0,*]
        timg = bytscl(timg, min=mosmincnt, max=mosmaxcnt, top=top)
        tut = uts[0]

        if keyword_set(test_grid) then begin
            tfn = 0
            sgopen, tfn, xsize=mosxsz, ysize=mosysz, /inch
    
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
            ; plot s/c footpoint and labels.
            for j=0, nprobe-1 do begin
                pre0 = 'rbsp'+probes[j]+'_'
                get_data, pre0+'fpt_mlat', tuts, scmlats
                get_data, pre0+'fpt_mlon', tuts, scmlons
                tmlt = scmlons[i,0]-midn[i]
                tlat = scmlats[i,0]
                plots, tmlt, tlat, /data, color=modcolors[0], $
                    psym=scsyms[j], symsize=symsz
                xyouts, tpos[0]+xchsz, tpos[1]+ychsz*(2+j), /normal, $
                    'RBSP-'+strupcase(probes[j]), color=white
                plots, tpos[0]+xchsz*8, tpos[1]+ychsz*(2+0.3+j), /normal, $
                    psym=scsyms[j], color=white, symsize=symsz
            endfor
            xyouts, tpos[0]+xchsz, tpos[1]+ychsz*(4+0), /normal, $
                strupcase(models[0]), color=modcolors[0]
            ; plot the vertical and horizontal grids on aurora.
            plots, reform(vfmlons[0,*,*])-midn[0], reform(vfmlats[0,*,*]), /data, $
                color=white, psym=3            
            plots, reform(hfmlons[0,*,*])-midn[0], reform(hfmlats[0,*,*]), /data, $
                color=6, psym=3
        endif
        
        
    ;---loop through each grid.
        vars = models[modidx]+['_v','_h']
        foreach tvar, vars do begin
        ;---load the grid to work with.
            get_data, tvar+'fmlat', tuts, fmlats, gridinfo
            get_data, tvar+'fmlon', tuts, fmlons, gridinfo
            
            picdir = rootdir+'/insitu_'+tvar            
            nrec = n_elements(uts)
            gridsz0 = size(fmlats,/dimensions)  ; orig spatial size.
            gridsz0 = gridsz0[1:*]
            gridsz1 = gridsz0*8                 ; full spatial size.

            
        ;---loop through each time to interpolate the grid at that time.
        ; saving the grids for all times cause memeory problem...
            for i=0, nrec-1 do begin
            ;---interpolate grid in time and space.
                tfmlats = dblarr(gridsz0)
                tfmlons = dblarr(gridsz0)
                for j=0, gridsz0[0]-1 do begin
                    for k=0, gridsz0[1]-1 do begin
                        tfmlats[j,k] = interpol(fmlats[*,j,k],tuts,uts[i],/quadratic)
                        tfmlons[j,k] = interpol(fmlons[*,j,k],tuts,uts[i],/quadratic)
                    endfor
                endfor
                ; interpolate to get the full spatial resolution.
                tfmlats = congrid(reform(tfmlats[*,*]), gridsz1[0],gridsz1[1],/interp)
                tfmlons = congrid(reform(tfmlons[*,*]), gridsz1[0],gridsz1[1],/interp)

            ;---convert mlat/mlon to image pixel position.
                tfrs = (90-tfmlats)/(90-minlat)
                tfts = (tfmlons-midn[i])*rad
                tfxs =  tfrs*sin(tfts)
                tfys = -tfrs*cos(tfts)
                txs = (tfxs+1)*imgsz[0]*0.5
                tys = tfys*imgsz[1]+imgsz[1]
                ; map using pixel position.
                timg0 = fltarr(imgsz)
                timg0[pxidx] = mos[i,*]
                timg1 = fltarr(gridsz1)
                for j=0, gridsz1[0]-1 do begin
                    for k=0, gridsz1[1]-1 do begin
                        tj = txs[j,k]
                        tk = tys[j,k]
                        timg1[j,k] = timg0[tj,tk]
                    endfor
                endfor
                ; export image.
                ofn = picdir+'/asf_insitu_'+time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
                ;ofn = 0
                xsz = gridsz1[0]*2
                ysz = gridsz1[1]*2*2
                sgopen, ofn, xsize=xsz, ysize=ysz
                erase, sgcolor('white')
                tpos = [0.15,0.35,0.95,0.75]
                device, decomposed=0
                loadct2, 1
                timg = bytscl(timg1, min=mosmincnt, max=mosmaxcnt, top=top)
                sgtv, timg, position=tpos, ct=1
                device, decomposed=1
                
                
                xchsz = double(!d.x_ch_size)/!d.x_size
                ychsz = double(!d.y_ch_size)/!d.y_size
                
                xr = minmax(gridinfo.xvalue)
                yr = minmax(gridinfo.yvalue)
                plot, xr, yr, /nodata, /noerase, position=tpos, title=gridinfo.title, $
                    xticklen=-0.02, xtitle=gridinfo.xtitle, xrange=xr, xstyle=1, $
                    yticklen=-0.01, ytitle=gridinfo.ytitle, yrange=yr, ystyle=1, yticks=2, yminor=5
                xyouts, xchsz, ychsz, /normal, time_string(uts[i])
                
                foreach tprobe, probes do begin
                    pre0 = 'rbsp'+tprobe+'_'
                    tcolor = (tprobe eq 'a')? sgcolor('red'): sgcolor('green')
                    tvar = pre0+'pos_gsm'
                    get_data, tvar, tmp, rgsm
                    rgsm = sinterpol(rgsm, tmp, uts[i])
                    tzzz = rgsm[2]
                    tddd = sqrt(rgsm[0]^2+rgsm[1]^2)
                    tmlt = -atan(rgsm[1],-rgsm[0])*deg/15
                    if strmid(gridinfo.title,0,1) eq 'V' then begin
                        plots, tmlt, tzzz, psym=1, color=tcolor
                    endif else begin
                        plots, tddd, tzzz, psym=1, color=tcolor
                    endelse
                endforeach
                xyouts, xchsz, 1-ychsz*1.5, /normal, 'RBSP-A', color=sgcolor('red')
                xyouts, xchsz, 1-ychsz*2.5, /normal, 'RBSP-B', color=sgcolor('green')
                sgclose
            endfor
        endforeach
    endif

end
