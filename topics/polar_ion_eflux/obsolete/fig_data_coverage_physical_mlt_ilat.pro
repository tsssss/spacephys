;+
; slicing according to distance.
;-


deg = 180d/!dpi
rad = !dpi/180


rootdir = shomedir()+'/polar_ion_eflux'
fns = rootdir+'/polar_ion_eflux_bin_info_physical_'+['1996','1997','1998']+'.sav'

;---Load bin info, contains binxs, binys, binzs, bininfos.
    restore, fns[0]
    binsize = size(bininfos,/dimensions)
    bincnts = dblarr(binsize)


    foreach tfn, fns do begin
        restore, tfn
        for i=0,binsize[0]-1 do $
            for j=0,binsize[1]-1 do $
                for k=0, binsize[2]-1 do begin
                    if ~ptr_valid(bininfos[i,j,k].uts) then continue
                    bincnts[i,j,k] += n_elements(*bininfos[i,j,k].uts)
                endfor

    endforeach



;---Settings.
    ct0 = 40
    top = 255
    thick = 2

    plotxrng = [-1,1]   ; polar coord.
    plotyrng = [-1,1]
    plotzrng = [2,9]    ; range in dis.
    minlat = 55d        ; in deg.

    zrng = [0,4]
    ztickv = smkarthm(zrng[0],zrng[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')
    ztitl = 'Count per bin!C '

    ; slice according to y.
    slicerngs = smkarthm(plotzrng[0],plotzrng[1],1,'dx')
    slicerngs = [2,3]
    nslicerng = n_elements(slicerngs)
    

;---Plot.
    for k=0, nslicerng-2 do begin
        kidx = where(bindiss ge slicerngs[k] and bindiss lt slicerngs[k+1], cnt)
        if cnt eq 0 then continue
        
        
        ofn = shomedir()+'/fig_data_coverage_r_'+$
            sgnum2str(slicerngs[k])+'_'+sgnum2str(slicerngs[k+1])+'_re.pdf'
        ;ofn = 0
        ;stop

        sgopen, ofn, xsize=6, ysize=6, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        pos0 = [0.1,0.05,0.9,0.85]
        pos1 = pos0 & pos1[1] = pos0[3]+ychsz*2 & pos1[3] = pos1[1]+ychsz*0.5
        
        
        plot, plotxrng, plotyrng, position=pos0, $
            xstyle=5, xrange=plotxrng, $
            ystyle=5, yrange=plotyrng, $
            /nodata, /noerase, /isotropic
        
        
        device, decomposed=0
        loadct, ct0
        tbincnts = total(bincnts[*,*,kidx],3)
        tdat = tbincnts
        tdat = alog10(tdat)
        tdat = bytscl(tdat, min=zrng[0], max=zrng[1], top=top)
        
        for i=0,binsize[0]-1 do $
            for j=0,binsize[1]-1 do begin
                tc = tdat[i,j]
                if tc le zrng[0] then continue

                mltrng = bininfos[i,j,kidx[0]].mlt_rng
                latrng = bininfos[i,j,kidx[0]].lat_rng
                rrng = abs(90-latrng)/(90-minlat)
                trng = mltrng*15*rad
                txs = [rrng[0]*cos(trng[0]),rrng[1]*cos(trng[0]),rrng[1]*cos(trng[1]),rrng[0]*cos(trng[1]),rrng[0]*cos(trng[0])]
                tys = [rrng[0]*sin(trng[0]),rrng[1]*sin(trng[0]),rrng[1]*sin(trng[1]),rrng[0]*sin(trng[1]),rrng[0]*sin(trng[0])]
                polyfill, txs, tys, /data, color=tc
            endfor

        
        ; draw axis.    
        lats = [minlat,smkarthm(65,85,10,'dx')]
        foreach tlat, lats do begin
            tmp = findgen(51)/50*2*!dpi
            tr = (90-tlat)/(90-minlat)
            txs = tr*cos(tmp)
            tys = tr*sin(tmp)
            plots, txs, tys, color=white, linestyle=1
            ti = 50*60d/360
            xyouts, txs[ti],tys[ti], /data, sgnum2str(tlat), color=white, alignment=0.5
        endforeach
            
        lons = smkarthm(0,360,30,'dx')
        foreach tlon, lons do begin
            tmp = smkarthm(plotxrng[0],0,50,'n')
            txs = tmp*cos(tlon*rad)
            tys = tmp*sin(tlon*rad)
            plots, txs, tys, color=white, linestyle=1
        endforeach
        lons = smkarthm(0,90,4,'x0')
        foreach tlon, lons do $
            xyouts, 1.03*cos(tlon*rad), 1.03*sin(tlon*rad)-0.02, /data, alignment=0.5, string(tlon/15,format='(I02)'), color=white
            
        ; draw colorbar.
        sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
            zrange=zrng, ztitle='R in '+sgnum2str(slicerngs[k])+' to '+sgnum2str(slicerngs[k+1])+' Re!C'+ztitl, $
            ztickv=ztickv, ztickname=ztickn
            
        sgclose

    endfor

end
