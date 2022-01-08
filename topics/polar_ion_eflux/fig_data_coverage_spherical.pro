

pro fig_data_coverage_spherical, bininfos=bininfos, tdat=tdat, bincnts=bincnts, noplot=noplot, slicerng=slicerng, _extra=extra

deg = 180d/!dpi
rad = !dpi/180


;---Settings.
    bintype = 'mlt_mlat_dis'

    ; slice on mlt.
    sliceidx = 0
    if n_elements(slicerng) eq 0 then slicerng = [-1,1]*1.5+12

    ; plot settings.
    ct0 = 40
    top = 255
    thick = 2

    case sliceidx of
        0: begin
            xtitl = 'X SM (Re)'
            ytitl = 'X SM (Re)'
            ztitl = 'MLT in ('+sgnum2str(mean(slicerng))+','+sgnum2str(mean(slicerng)-12)+')+/-1.5 hr!CCount per bin'
            xrang = [-1d,1]*7
            yrang = [0d,1]*9
            zrang = [0d,5]
            minlat = 45d
            end
    endcase

    ztickv = smkarthm(zrang[0],zrang[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')

    xsz = 6d ; inch
    ysz = abs((yrang[1]-yrang[0])/(xrang[1]-xrang[0])*xsz)

    ; filenames.
    stryrs = string(smkarthm(1996,1998,1,'dx'),format='(I4)')
    rootdir = shomedir()+'/polar_ion_eflux'
    datadir = shomedir()+'/polar_ion_eflux/data'
    datafns = datadir+'/'+stryrs+'/polar_ion_eflux_bin_info_'+bintype+'_'+stryrs+'.sav'

    figfn = shomedir()+'/fig_data_coverage_'+bintype+'.pdf'


;---Load bin info, contains bin1s, bin2s, bin3s, bininfos.
    restore, datafns[0]
    binsz = size(bininfos,/dimensions)
    nbin1 = binsz[0]
    nbin2 = binsz[1]
    nbin3 = binsz[2]
    bincnts = dblarr(binsz)


    foreach tfn, datafns do begin
        restore, tfn
        for i=0,nbin1-1 do for j=0,nbin2-1 do for k=0, nbin3-1 do begin
            if ~ptr_valid(bininfos[i,j,k].uts) then continue
            bincnts[i,j,k] += n_elements(*bininfos[i,j,k].uts)
        endfor
    endforeach



;---reduce bin dimension.
    case sliceidx of
        0: begin
            tbins = bin1s
            tbins[where(tbins ge 18)]-= 24
            end
        1: tbins = bin2s
        2: tbins = bin3s
    endcase
    kidx = where(tbins ge slicerng[0] and tbins le slicerng[1], cnt)
    if cnt eq 0 then return

    case sliceidx of
        0: bincnts = bincnts[kidx,*,*]
        1: bincnts = bincnts[*,kidx,*]
        2: bincnts = bincnts[*,*,kidx]
    endcase
    tdat = total(bincnts, sliceidx+1)
    tdat = alog10(tdat)
    tdat = bytscl(tdat, min=zrang[0], max=zrang[1], top=top)

    case sliceidx of
        0: bininfos = reform(bininfos[kidx[0],*,*])
        1: bininfos = reform(bininfos[*,kidx[0],*])
        2: bininfos = reform(bininfos[*,*,kidx[0]])
    endcase

    tsz = size(bininfos,/dimensions)
    if keyword_set(noplot) then return
    
    
;---Plot.
;    figfn = 0
;    stop

    sgopen, figfn, xsize=xsz, ysize=ysz, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    pos0 = [0.15,0.1,0.9,0.85]
    pos1 = pos0 & pos1[1] = pos0[3]+ychsz*1 & pos1[3] = pos1[1]+ychsz*0.5
    
    
    plot, xrang, yrang, position=pos0, $
        xstyle=1, xrange=xrang, xtitle=xtitl, $
        ystyle=1, yrange=yrang, ytitle=ytitl, $
        /nodata, /noerase, /isotropic
    
    
    device, decomposed=0
    loadct, ct0
    

    for i=0,tsz[0]-1 do for j=0,tsz[1]-1 do begin
        tc = tdat[i,j]
        if tc le zrang[0] then continue
    
        case sliceidx of
            0: begin
                tyyrng = bininfos[i,j].range2
                txxrng = bininfos[i,j].range3
                end
            1: begin

                end
            2: begin

                end
        endcase
        tyyrng = tyyrng*rad
        txs = [txxrng[0]*cos(tyyrng[0]),txxrng[0]*cos(tyyrng[1]),txxrng[1]*cos(tyyrng[1]),txxrng[1]*cos(tyyrng[0]),txxrng[0]*cos(tyyrng[0])]
        tys = [txxrng[0]*sin(tyyrng[0]),txxrng[0]*sin(tyyrng[1]),txxrng[1]*sin(tyyrng[1]),txxrng[1]*sin(tyyrng[0]),txxrng[0]*sin(tyyrng[0])]
        polyfill, txs, tys, /data, color=tc
    endfor
    
    fig_data_coverage_spherical, bininfos=bininfos, tdat=tdat, bincnts=bincnts, /noplot, slicerng=slicerng-12
    for i=0,tsz[0]-1 do for j=0,tsz[1]-1 do begin
        tc = tdat[i,j]
        if tc le zrang[0] then continue

        case sliceidx of
            0: begin
                tyyrng = bininfos[i,j].range2
                txxrng = bininfos[i,j].range3
            end
            1: begin

            end
            2: begin

            end
        endcase
        tyyrng = tyyrng*rad
        txs =-[txxrng[0]*cos(tyyrng[0]),txxrng[0]*cos(tyyrng[1]),txxrng[1]*cos(tyyrng[1]),txxrng[1]*cos(tyyrng[0]),txxrng[0]*cos(tyyrng[0])]
        tys = [txxrng[0]*sin(tyyrng[0]),txxrng[0]*sin(tyyrng[1]),txxrng[1]*sin(tyyrng[1]),txxrng[1]*sin(tyyrng[0]),txxrng[0]*sin(tyyrng[0])]
        polyfill, txs, tys, /data, color=tc
    endfor
    

    ; draw earth.
    device, decompose=1
    tcolor = sgcolor('black')
    tmp = findgen(51)/50*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    idx = where(txs lt 0)
    polyfill, [txs[idx],0], [tys[idx],0], /line_fill, color=tcolor, orientation=45
    polyfill, [txs[idx],0], [tys[idx],0], /line_fill, color=tcolor, orientation=-45
    plots, txs, tys, /data, color=tcolor, thick=thick

    plot, xrang, yrang, position=pos0, $
        xstyle=1, xrange=xrang, xtitle=xtitl, $
        ystyle=1, yrange=yrang, ytitle=ytitl, $
        /nodata, /noerase, /isotropic
    
    
    ; draw axis.
    for i=2,8,2 do oplot, txs*i, tys*i, linestyle=1, color=tcolor
    for i=0,180,15 do oplot, [1,15]*cos(i*rad), [1,15]*sin(i*rad), linestyle=1, color=tcolor


    ; draw colorbar.
    sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
        zrange=zrang, ztitle=ztitl, $
        ztickv=ztickv, ztickname=ztickn
    
    sgclose


end
