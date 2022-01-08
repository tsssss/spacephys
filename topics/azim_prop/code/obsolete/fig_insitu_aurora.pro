
;---Settings.
    id = '2014_0828_10'
    rootdir = sparentdir(srootdir())
    datadir = rootdir+'/data' & if file_test(datadir,/directory) eq 0 then file_mkdir, datadir
    plotdir = rootdir+'/plot' & if file_test(plotdir,/directory) eq 0 then file_mkdir, plotdir
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    deg = 180d/!dpi
    rad = !dpi/180d
    cns = -1
    top = 254

    
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero'
    tplot_options, 'ystyle', 1
    tplot_options, 'yticklen', -0.01
    tplot_options, 'xticklen', -0.03

    
    utr2 = time_double(['2014-08-28/10:05','2014-08-28/10:20']) ; time to map aurora to in-situ.

    

;---Set up the in-situ grid.
    nrec = n_elements(uts)
    ; a vertical grid at D = 6Re.
    mltgridrng = [-1,2]     ; in hr.
    mltgrids = smkarthm(mltgridrng[0],mltgridrng[1], 0.05, 'dx')
    nmltgrid = n_elements(mltgrids)
    zzzgridrng = [0,4.4]      ; in Re.
    zzzgrids = smkarthm(zzzgridrng[0],zzzgridrng[1], 0.10, 'dx')
    nzzzgrid = n_elements(zzzgrids)
    vgridd = 5d

;---Map in-situ grid to ionosphere.
    re = 6378d  ; km.
    r0 = 110d/re+1  ; in Re.
    dir = -1    ; to northern hemisphere.
    par = 2d
    tut = time_double('2014-08-28/10:13')
    tet = stoepoch(tut,'unix')
    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date, tilt=tilt
    grid_fmlat = fltarr(nmltgrid, nzzzgrid) ; in deg.
    grid_fmlon = fltarr(nmltgrid, nzzzgrid) ; in deg.
    for i=0, nmltgrid-1 do begin
        for j=0, nzzzgrid-1 do begin
            tang = mltgrids[i]*15*rad    ; in rad.
            ; xyz in sm.
            x0 = -vgridd*cos(tang)
            y0 = -vgridd*sin(tang)
            z0 = zzzgrids[j]
            geopack_conv_coord, x0,y0,z0, /from_sm, $
                x1,y1,z1, /to_gsm
            geopack_trace, x1,y1,z1, dir, par, x2,y2,z2, r0=r0, $
                /refine, /ionosphere, /t89
            geopack_conv_coord, x2,y2,z2, /from_gsm, $
                x3,y3,z3, /to_mag
            grid_fmlat[i,j] = asin(z3/r0)*deg
            grid_fmlon[i,j] = atan(y3,x3)*deg
        endfor
    endfor
    grid_fmlt = slon2lt(grid_fmlon, tet, /mag, /deg)


    
;---Load data, need mos.
    _2014_0828_load_data
    get_data, 'asf_mos', uts, mos, pxidx
    get_data, 'asf_info', 0, asinfo
    mlts = asinfo.mlts
    mlats = asinfo.mlats
    imgsz = asinfo.imgsz
    minlat = asinfo.minlat
    idx = where(uts ge utr2[0] and uts le utr2[1])
    uts = uts[idx]
    mos = mos[idx,*]
    

;---Test grid on aurora.
    test_grid_on_aurora = 0
    if keyword_set(test_grid_on_aurora) then begin
        timg = fltarr(imgsz)
        timg[pxidx] = mos[where(uts eq tut),*]
        timg = bytscl(timg, max=200, top=top)
        tpos = [0,0,1,1]
        red = sgcolor('red')
        white = sgcolor('white')
        
        sgopen, 0, xsize=imgsz[0], ysize=imgsz[1]
        sgtv, timg, ct=1, position=tpos
        sgset_map, position=tpos, color=white, xrange=[-1,1]*90, yrange=[minlat,90]
        
        tmp = findgen(11)/10*2*!dpi
        txs = cos(tmp)
        tys = sin(tmp)
        usersym, txs, tys, /fill
        plots, grid_fmlt, grid_fmlat, color=red, psym=3
        xyouts, 10,10, /device, time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color=white
    endif
    
    
    
;---map aurora to in-situ grid.
    gridsz1 = [nmltgrid, nzzzgrid]
    
    tfrs = (90-grid_fmlat)/(90-minlat)
    tfts = grid_fmlt*rad
    tfxs =  tfrs*sin(tfts)
    tfys = -tfrs*cos(tfts)
    txs = (tfxs+1)*imgsz[0]*0.5
    tys = tfys*imgsz[1]+imgsz[1]
    
    tpos = [0.1,0.1,0.95,0.95]
    azimgridrng = mltgridrng*15*rad*vgridd
    scale = (zzzgridrng[1]-zzzgridrng[0])/(azimgridrng[1]-azimgridrng[0])
    figxsz = 500  ; pixel, roughly 5 inch.
    figysz = figxsz*scale
    
    fps = 10d
    vfn = plotdir+'/asf_img_insitu.mp4'
    ovid = idlffvideowrite(vfn)
    vidstream = ovid.AddVideoStream(figxsz, figysz, fps)
    
    ; map aurora to in-situ grid using pixel position.
    ofn = rootdir+'/fig_insitu_tmp.png'
    foreach tut, uts, i do begin
        printf, cns, time_string(tut)
        
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
        
        sgopen, ofn, xsize=figxsz, ysize=figysz
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        device, decomposed=1
        white = sgcolor('white')
        
        timg = bytscl(timg1, max=200, top=top)
        
        sgtv, timg, position=tpos, ct=1, /resize
        plot, azimgridrng, zzzgridrng, /nodata, /noerase, position=tpos, /isotropic, $
            xstyle=1, xticklen=-0.01, xtitle='Azim (Re)', $
            ystyle=1, yticklen=-0.01*scale, ytitle='Z (Re)', $
            color=white
        xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz*1.2, /normal, $
            time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color=white
        sgclose
        
        read_png, ofn, timg
        time = ovid.put(vidstream, timg)
    endforeach
    
    

    ovid.Cleanup
    file_delete, ofn

end
