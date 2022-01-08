;+
; slicing according to mlt.
;-


deg = 180d/!dpi
rad = !dpi/180


rootdir = shomedir()+'/polar_ion_eflux'
fns = rootdir+'/polar_ion_eflux_bin_info_spherical_'+['1996','1997','1998']+'.sav'

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

    plotxrng = [-6,6]  ; in deg.
    plotyrng = [0,9]    ; in dis in re.
    plotzrng = [11,13]  ; range in mlt.
    minlat = 55d        ; in deg.

    zrng = [0,4]
    ztickv = smkarthm(zrng[0],zrng[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')
    ztitl = 'Count per bin!C '

    ; slice according to mlt.
    slicerngs = smkarthm(plotzrng[0],plotzrng[1],1,'dx')
    slicerngs = [-1,1]*2+12
    nslicerng = n_elements(slicerngs)
    

;---Plot.
    for k=0, nslicerng-2 do begin
        kidx = where(binmlts ge slicerngs[k] and binmlts lt slicerngs[k+1], cnt)
        if cnt eq 0 then continue
        tbincnts = total(bincnts[kidx,*,*],1)
        
        
        ofn = shomedir()+'/fig_data_coverage_mlt_'+$
            sgnum2str(slicerngs[k])+'_'+sgnum2str(slicerngs[k+1])+'_hr.pdf'
;        ofn = 0
;        stop

        sgopen, ofn, xsize=6, ysize=5, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        pos0 = [0.1,0.1,0.9,0.85]
        pos1 = pos0 & pos1[1] = pos0[3]+ychsz*1 & pos1[3] = pos1[1]+ychsz*0.5
        
        
        plot, plotxrng, plotyrng, position=pos0, $
            xstyle=1, xrange=plotxrng, xtitle='X SM (Re)', $
            ystyle=1, yrange=plotyrng, ytitle='Y SM (Re)', $
            /nodata, /noerase, /isotropic
        
        
        device, decomposed=0
        loadct, ct0
        tdat = tbincnts
        tdat = alog10(tdat)
        tdat = bytscl(tdat, min=zrng[0], max=zrng[1], top=top)
        
        for i=0,binsize[1]-1 do $
            for j=0,binsize[2]-1 do begin
                tc = tdat[i,j]
                if tc le zrng[0] then continue

                latrng = bininfos[kidx[0],i,j].lat_rng*rad
                disrng = bininfos[kidx[0],i,j].dis_rng
                txs = [disrng[0]*cos(latrng[0]),disrng[1]*cos(latrng[0]),disrng[1]*cos(latrng[1]),disrng[0]*cos(latrng[1]),disrng[0]*cos(latrng[0])]
                tys = [disrng[0]*sin(latrng[0]),disrng[1]*sin(latrng[0]),disrng[1]*sin(latrng[1]),disrng[0]*sin(latrng[1]),disrng[0]*sin(latrng[0])]
                polyfill, txs, tys, /data, color=tc
            endfor


        ; the other half.
        idx = where(binmlts lt 12) & binmlts[idx]+= 24
        kidx = where(binmlts ge slicerngs[k]+12 and binmlts lt slicerngs[k+1]+12, cnt)
        if cnt eq 0 then continue
        tbincnts = total(bincnts[kidx,*,*],1)
        
        tdat = tbincnts
        tdat = alog10(tdat)
        tdat = bytscl(tdat, min=zrng[0], max=zrng[1], top=top)
        
        for i=0,binsize[1]-1 do $
            for j=0,binsize[2]-1 do begin
                tc = tdat[i,j]
                if tc le zrng[0] then continue

                latrng = bininfos[kidx[0],i,j].lat_rng*rad
                disrng = bininfos[kidx[0],i,j].dis_rng
                txs =-[disrng[0]*cos(latrng[0]),disrng[1]*cos(latrng[0]),disrng[1]*cos(latrng[1]),disrng[0]*cos(latrng[1]),disrng[0]*cos(latrng[0])]
                tys = [disrng[0]*sin(latrng[0]),disrng[1]*sin(latrng[0]),disrng[1]*sin(latrng[1]),disrng[0]*sin(latrng[1]),disrng[0]*sin(latrng[0])]
                polyfill, txs, tys, /data, color=tc
            endfor
        

        
        ; draw earth.
        device, decompose=1
        white = sgcolor('white')
        black = sgcolor('black')
        tmp = findgen(51)/50*!dpi
        txs = cos(tmp)
        tys = sin(tmp)
        ;idx = where(txs ge 0)
        ;polyfill, txs[idx], tys[idx], color=sgcolor('white')
        idx = where(txs lt 0)
        ;polyfill, txs[idx], tys[idx], color=sgcolor('silver')
        polyfill, [txs[idx],0], [tys[idx],0], /line_fill, color=black, orientation=45
        polyfill, [txs[idx],0], [tys[idx],0], /line_fill, color=black, orientation=-45
        plots, txs, tys, /data, color=black, thick=thick

        for i=1, 9 do oplot, txs*i, tys*i, linestyle=1

            
        ; draw colorbar.
        sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
            zrange=zrng, ztitle='MLT in '+sgnum2str(slicerngs[k])+' to '+sgnum2str(slicerngs[k+1])+' hr!C'+ztitl, $
            ztickv=ztickv, ztickname=ztickn
            
        sgclose

    endfor

end
