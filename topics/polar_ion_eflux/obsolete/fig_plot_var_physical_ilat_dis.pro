;+
; slicing according to mlt.
; plot physical quantities.
;-


deg = 180d/!dpi
rad = !dpi/180


yrs = ['1996','1997','1998']
nyr = n_elements(yrs)
var0 = 'h_density'
var0 = 'o_density'

rootdir = shomedir()+'/polar_ion_eflux'
binfns = rootdir+'/polar_ion_eflux_bin_info_physical_'+yrs+'.sav'
datfns = rootdir+'/polar_ion_eflux_'+var0+'_'+yrs+'.sav'


;---Load bin info, contains binxs, binys, binzs, bininfos.
    restore, binfns[0]
    binsize = size(bininfos,/dimensions)
    bincnts = dblarr(binsize)
    bindata = dblarr(binsize)


    for ii=0, nyr-1 do begin
        restore, binfns[ii]
        restore, datfns[ii]
        for i=0,binsize[0]-1 do for j=0,binsize[1]-1 do for k=0, binsize[2]-1 do begin
            if ~ptr_valid(bininfos[i,j,k].uts) then continue
            bincnts[i,j,k] += n_elements(*bininfos[i,j,k].uts)
            bindata[i,j,k] += total(dat0[*bininfos[i,j,k].idx])
        endfor
    endfor



;---Settings.
    ct0 = 40
    top = 255
    thick = 2

    plotxrng = [55,90]  ; in deg.
    plotyrng = [1,9]    ; in dis in re.
    plotzrng = [11,13]  ; range in mlt.
    minlat = 55d        ; in deg.

    zrng = [0,2]
    ztickv = smkarthm(zrng[0],zrng[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')
    ztitl = 'Density (cm!U-3!N)!C '

    ; slice according to mlt.
    slicerngs = smkarthm(plotzrng[0],plotzrng[1],1,'dx')
    slicerngs = [-1,1]*1.5+12
    nslicerng = n_elements(slicerngs)
    

;---Plot.
    for k=0, nslicerng-2 do begin
        kidx = where(binmlts ge slicerngs[k] and binmlts lt slicerngs[k+1], cnt)
        if cnt eq 0 then continue
        
        
        ofn = shomedir()+'/fig_data_coverage_mlt_'+$
            sgnum2str(slicerngs[k])+'_'+sgnum2str(slicerngs[k+1])+'_hr.pdf'
        ofn = 0
        stop

        sgopen, ofn, xsize=6, ysize=6, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        pos0 = [0.1,0.05,0.9,0.85]
        pos1 = pos0 & pos1[1] = pos0[3]+ychsz*2 & pos1[3] = pos1[1]+ychsz*0.5
        
        
        plot, plotxrng, plotyrng, position=pos0, $
            xstyle=1, xrange=plotxrng, xtitle='ILat (deg)', $
            ystyle=1, yrange=plotyrng, ytitle='R (Re)', $
            /nodata, /noerase
        
        
        device, decomposed=0
        loadct, ct0
        tbincnts = total(bincnts[kidx,*,*],1)
        tbindata = total(bindata[kidx,*,*],1)
        tdat = tbindata/tbincnts
        tdat = alog10(tdat)
        tdat = bytscl(tdat, min=zrng[0], max=zrng[1], top=top)
        
        for i=0,binsize[1]-1 do $
            for j=0,binsize[2]-1 do begin
                tc = tdat[i,j]
                if tc le zrng[0] then continue

                latrng = bininfos[kidx[0],i,j].lat_rng
                disrng = bininfos[kidx[0],i,j].dis_rng
                txs = [latrng[0],latrng[1],latrng[1],latrng[0],latrng[0]]
                tys = [disrng[0],disrng[0],disrng[1],disrng[1],disrng[0]]
                polyfill, txs, tys, /data, color=tc
            endfor

        
        ; draw earth.    

            
        ; draw colorbar.
        sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
            zrange=zrng, ztitle='MLT in '+sgnum2str(slicerngs[k])+' to '+sgnum2str(slicerngs[k+1])+' hr!C'+ztitl, $
            ztickv=ztickv, ztickname=ztickn
            
        sgclose

    endfor

end
