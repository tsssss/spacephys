

;---Settings.
    pres = ['rbb','rba','tha','the','thd']+'_'
    colors = [1,3,4,5,6]
    
    models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    modelsyms = [1,6,5,7]
    
    ;_2014_1222_load_data


;---Plot mosaic.
    tut = time_double('2014-12-22/05:18')
    asifn = sread_thg_mosaic(tut+[-3,0], sites, type='asf', minlat=minlat, dark=0, /notop)
    cdfs = scdfread(asifn)

    uts = (*cdfs[0].value)
    mos = (*cdfs[1].value)
    midn= (*cdfs[2].value)
    mlt = (*cdfs[3].value)
    mlat= (*cdfs[4].value)
    imgsz = (*cdfs[5].value)
    pxidx = (*cdfs[6].value)
    minlat = (*cdfs[7].value)[0]

    img = fltarr(imgsz)
    nrec = n_elements(uts)
    img[pxidx] = mos[1,*]

    ofn = shomedir()+'/fig_mosaic.pdf'
    ;ofn = 0
    sgopen, ofn, xsize=imgsz[0], ysize=imgsz[1]
    polyfill, [0,1,1,0,0], [0,0,1,1,0], /normal, color=sgcolor('black')
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tpos = [0.05,0.05,0.95,0.95]
    
    device, decomposed=0
    loadct, 1
    sgtv, bytscl(img, min=10, max=400, top=254), position=tpos

    loadct2, 43
    xyouts, xchsz*1, ychsz*0.2, /normal, time_string(tut), color=255
    foreach pre0, pres, i do $
        xyouts, xchsz*1, ychsz*0.2+(i+1)*ychsz, /normal, strupcase(strmid(pre0,0,3)), color=colors[i]
    foreach tmodel, models, i do begin
        xyouts, xchsz*8, ychsz*0.2+(i+1)*ychsz, /normal, strupcase(tmodel), color=255
        plots, xchsz*8+xchsz*5, ychsz*0.2+(i+1)*ychsz+0.4*ychsz, /normal, psym=modelsyms[i], color=255
    endforeach

    
    sgset_map, position=tpos, xrange=[-1,1]*90, yrange=[minlat,90], /noxtick, $
        color=255, ytickv=smkarthm(minlat,90,5,'dx'), xtickv=smkarthm(-90,90,30,'dx')
    for i=0,2 do xyouts, -30, 60+i*10, /data, string(60+i*10,format='(I02)'), alignment=0.5, color=255
    for i=0,2 do xyouts, i*90-90, 59, /data, string((18+i*6) mod 24,format='(I02)'), alignment=0.5, color=255
    

    mlts = dblarr(nmodel+1,n_elements(pres))
    foreach pre0, pres, i do begin
        get_data, pre0+'fpt_mlat', uts, mlat
        mlat = sinterpol(mlat, uts, tut)
        get_data, pre0+'fpt_mlon', uts, mlon
        mlon = sinterpol(mlon, uts, tut)
        tx = mlon-midn[1]
        ty = mlat
        mlts[i,1:*] = tx/15
        get_data, pre0+'mlt', uts, mlt
        mlts[i,0] = interpol(mlt,uts,tut)
        for j=0, nmodel-1 do plots, tx[j], ty[j], psym=modelsyms[j], color=colors[i]
        ;xyouts, tx, ty, /data, alignment=0, '  '+strupcase(strmid(pre0,0,3)), color=colors[i]
    endforeach
    
    sgclose
    
    print, 'RBB-RBA', string(mean(mlts[1,*]-mlts[0,*]),format='(F5.2)'), string(stddev(mlts[1,*]-mlts[0,*]),format='(F5.2)')
    print, 'THA-THE', string(mean(mlts[2,*]-mlts[3,*]),format='(F5.2)'), string(stddev(mlts[2,*]-mlts[3,*]),format='(F5.2)')
    print, 'THE-THD', string(mean(mlts[3,*]-mlts[4,*]),format='(F5.2)'), string(stddev(mlts[3,*]-mlts[4,*]),format='(F5.2)')
    
    
    
end
