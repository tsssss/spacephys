
deg = 180d/!dpi
rad = !dpi/180


rootdir = shomedir()+'/polar_ion_eflux'
fns = rootdir+'/polar_ion_eflux_bin_info_cartesian_'+['1996','1997','1998']+'.sav'

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

    plotxrng = [-7.5,6.5]
    plotyrng = [-7.5,6.5]
    plotzrng = [-3,9]

    zrng = [0,5]
    ztickv = smkarthm(zrng[0],zrng[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')
    ztitl = 'Count per bin!C '

    ; slice according to y.
    slicerngs = smkarthm(plotyrng[0],plotyrng[1],1,'dx')
    slicerngs = [-1,1]*1.5
    nslicerng = n_elements(slicerngs)


;---Plot.
    for k=0, nslicerng-2 do begin
        kidx = where(binys ge slicerngs[k] and binys le slicerngs[k+1], cnt)
        if cnt eq 0 then continue
        
        
        ofn = shomedir()+'/fig_data_coverage_cartesian_y_'+$
            string(slicerngs[k],format='(F4.1)')+'_'+string(slicerngs[k+1],format='(F4.1)')+'_re.pdf'
        ;ofn = 0
        ;stop

        sgopen, ofn, xsize=6, ysize=5, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        pos0 = [0.15,0.1,0.9,0.85]
        pos1 = pos0 & pos1[1] = pos0[3]+ychsz*1 & pos1[3] = pos1[1]+ychsz*0.5
        
        
        plot, plotxrng, plotzrng, position=pos0, $
            xstyle=1, xrange=plotxrng, xtitle='X GSM (Re)', $
            ystyle=1, yrange=plotzrng, ytitle='Z GSM (Re)', $
            /nodata, /noerase, /isotropic
        
        
        device, decomposed=0
        loadct, ct0
        tbincnts = total(bincnts[*,kidx,*],2)
        tdat = tbincnts
        tdat = alog10(tdat)
        tdat = bytscl(tdat, min=zrng[0], max=zrng[1], top=top)
        
        for i=0,binsize[0]-1 do $
            for j=0,binsize[2]-1 do begin
                tc = tdat[i,j]
                if tc le zrng[0] then continue
            
                txxrng = bininfos[i,kidx[0],j].x_rng
                tyyrng = bininfos[i,kidx[0],j].z_rng
                txs = [txxrng[0],txxrng[0],txxrng[1],txxrng[1],txxrng[0]]
                tys = [tyyrng[0],tyyrng[1],tyyrng[1],tyyrng[0],tyyrng[0]]
                polyfill, txs, tys, /data, color=tc
            endfor
        
        
        




        ; draw earth.
        device, decompose=1
        white = sgcolor('white')
        black = sgcolor('black')
        tmp = findgen(51)/50*2*!dpi
        txs = cos(tmp)
        tys = sin(tmp)
        ;idx = where(txs ge 0)
        ;polyfill, txs[idx], tys[idx], color=sgcolor('white')
        idx = where(txs lt 0)
        ;polyfill, txs[idx], tys[idx], color=sgcolor('silver')
        polyfill, txs[idx], tys[idx], /line_fill, color=black, orientation=45
        polyfill, txs[idx], tys[idx], /line_fill, color=black, orientation=-45
        plots, txs, tys, /data, color=black, thick=thick









        ; draw colorbar.
        sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
            zrange=zrng, ztitle='Y GSM in '+sgnum2str(slicerngs[k])+' to '+sgnum2str(slicerngs[k+1])+' Re!C'+ztitl, $
            ztickv=ztickv, ztickname=ztickn
            
        sgclose

    endfor

end
