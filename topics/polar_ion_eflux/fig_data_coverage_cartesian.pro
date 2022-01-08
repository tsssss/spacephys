
deg = 180d/!dpi
rad = !dpi/180


;---Settings.
    bintype = 'gsm_xyz'

    ; slice on gsm_y.
    sliceidx = 1
    slicerng = [-1,1]*1.5

    ; plot settings.
    ct0 = 40
    top = 255
    thick = 2

    case sliceidx of
        1: begin
            xtitl = 'X GSM (Re)'
            ytitl = 'Z GSM (Re)'
            ztitl = 'Y GSM in ('+sgnum2str(slicerng[0])+','+sgnum2str(slicerng[1])+') Re!CCount per bin'
            xrang = [-7.5,6.5]
            yrang = [-3,9]
            xrang = [-1,1]*10
            yrang = [-1,1]*10
            zrang = [0,5]
            end
    endcase

    ztickv = smkarthm(zrang[0],zrang[1],1,'dx')
    ztickn = '10!U'+string(ztickv,format='(I0)')

    xsz = 6 ; inch
    ysz = abs((yrang[1]-yrang[0])/(xrang[1]-xrang[0])*xsz)

    ; filenames.
    stryrs = string(smkarthm(1996,2008,1,'dx'),format='(I4)')
    rootdir = shomedir()+'/polar_ion_eflux'
    datadir = shomedir()+'/polar_ion_eflux/data'
    datafns = datadir+'/'+stryrs+'/polar_ion_eflux_bin_info_'+bintype+'_'+stryrs+'.sav'
    datafns = datadir+'/'+stryrs+'/polar_ion_eflux_bin_info_hydra_'+bintype+'_'+stryrs+'.sav'

    figfn = shomedir()+'/fig_data_coverage_hydra_'+bintype+'.pdf'


;---Load bin info, contains bin1s, bin2s, bin3s, bininfos.
    restore, datafns[0]
    binsz = size(bininfos,/dimensions)
    nbin1 = binsz[0]
    nbin2 = binsz[1]
    nbin3 = binsz[2]
    bincnts = dblarr(binsz)


    foreach tfn, datafns do begin
        if file_test(tfn) eq 0 then continue
        restore, tfn
        for i=0,nbin1-1 do for j=0,nbin2-1 do for k=0, nbin3-1 do begin
            if ~ptr_valid(bininfos[i,j,k].uts) then continue
            bincnts[i,j,k] += n_elements(*bininfos[i,j,k].uts)
        endfor
    endforeach



;---reduce bin dimension.
    case sliceidx of
        0: tbins = bin1s
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
    
    
;---Plot.
    ;figfn = 0
    ;stop

    sgopen, figfn, xsize=xsz, ysize=ysz, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    pos0 = [0.15,0.1,0.9,0.85]
    pos1 = pos0 & pos1[1] = pos0[3]+ychsz*1 & pos1[3] = pos1[1]+ychsz*0.5
    
    
    plot, xrang, yrang, position=pos0, $
        xstyle=5, xrange=xrang, xtitle=xtitl, $
        ystyle=5, yrange=yrang, ytitle=ytitl, $
        /nodata, /noerase, /isotropic
    
    
    device, decomposed=0
    loadct, ct0
    

    for i=0,tsz[0]-1 do for j=0,tsz[1]-1 do begin
        tc = tdat[i,j]
        if tc le zrang[0] then continue
    
        case sliceidx of
            0: begin
                txxrng = bininfos[i,j].range2
                tyyrng = bininfos[i,j].range3
                end
            1: begin
                txxrng = bininfos[i,j].range1
                tyyrng = bininfos[i,j].range3
                end
            2: begin
                txxrng = bininfos[i,j].range1
                tyyrng = bininfos[i,j].range2
                end
        endcase
        txs = [txxrng[0],txxrng[0],txxrng[1],txxrng[1],txxrng[0]]
        tys = [tyyrng[0],tyyrng[1],tyyrng[1],tyyrng[0],tyyrng[0]]
        polyfill, txs, tys, /data, color=tc
    endfor
    


    ; draw earth.
    device, decompose=1
    tcolor = sgcolor('black')
    tmp = findgen(51)/50*2*!dpi
    txs = cos(tmp)
    tys = sin(tmp)
    ;idx = where(txs ge 0)
    ;polyfill, txs[idx], tys[idx], color=sgcolor('white')
    idx = where(txs lt 0)
    ;polyfill, txs[idx], tys[idx], color=sgcolor('silver')
    polyfill, txs[idx], tys[idx], /line_fill, color=tcolor, orientation=45
    polyfill, txs[idx], tys[idx], /line_fill, color=tcolor, orientation=-45
    plots, txs, tys, /data, color=tcolor, thick=thick

    plot, xrang, yrang, position=pos0, $
        xstyle=1, xrange=xrang, xtitle=xtitl, $
        ystyle=1, yrange=yrang, ytitle=ytitl, $
        /nodata, /noerase, /isotropic


    ; draw colorbar.
    sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
        zrange=zrang, ztitle=ztitl, $
        ztickv=ztickv, ztickname=ztickn
    
    sgclose


end