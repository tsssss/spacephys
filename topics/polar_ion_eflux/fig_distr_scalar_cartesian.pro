

pro fig_distr_scalar_cartesian, dattype, bintype


;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180




;---Settings.
    if n_elements(dattype) eq 0 then message, 'no data type ...'
    if n_elements(bintype) eq 0 then bintype = 'gsm_xyz'
    if n_elements(statype) eq 0 then statype = 'average'

    ; slice on gsm_y.
    sliceidx = 1
    slicerng = [-1,1]*1.5

    ; plot settings.
    ct0 = 40
    top = 254
    thick = 2

    case dattype of
        'h_density': begin
            unit = 'Log10 H!U+!N density (cm!U-3!N)'
            zrang = [-1,2.5]
            end
        'o_density': begin
            unit = 'Log10 O!U+!N density (cm!U-3!N)'
            zrang = [-1,2.5]
            end
        'h_t_avg': begin
            unit = 'Log10 H!U+!N Temperature (eV)'
            zrang = [2.5,4]-1
            end
        'o_t_avg': begin
            unit = 'Log10 O!U+!N Temperature (eV)'
            zrang = [2.5,4]-1
            end
    endcase
    
    case bintype of
        'gsm_xyz': begin
            case sliceidx of
                1: begin
                    xtitl = 'X GSM (Re)'
                    ytitl = 'Z GSM (Re)'
                    ztitl = 'Y GSM in ('+sgnum2str(slicerng[0])+','+sgnum2str(slicerng[1])+') Re!C'+unit
                    xrang = [-7.5,6.5]
                    yrang = [-3,9]
                end
            endcase
            end
    endcase
    

    ztickv = smkarthm(zrang[0],zrang[1],0.5,'dx')
    nztick = n_elements(ztickv)
    ztickn = string(ztickv,format='(F4.1)')




    xsz = 6 ; inch
    ysz = abs((yrang[1]-yrang[0])/(xrang[1]-xrang[0])*xsz)

    ; filenames.
    stryrs = ['1996','1997','1998']
    nyr = n_elements(stryrs)
    rootdir = shomedir()+'/polar_ion_eflux'
    datadir = shomedir()+'/polar_ion_eflux/data'
    binfns = datadir+'/'+stryrs+'/polar_ion_eflux_bin_info_'+bintype+'_'+stryrs+'.sav'
    datfns = datadir+'/'+stryrs+'/polar_ion_eflux_'+dattype+'_'+stryrs+'.sav'
    utsfns = datadir+'/'+stryrs+'/polar_ion_eflux_uts_'+stryrs+'.sav'


    figfn = shomedir()+'/fig_distr_'+dattype+'_'+bintype+'.pdf'


;---Load bin info, contains bin1s, bin2s, bin3s, bininfos.
    restore, binfns[0]

    ; reduce bin dimension.
    case sliceidx of
        0: begin
            kidx = where(bin1s ge slicerng[0] and bin1s le slicerng[1], nslice)
            if nslice eq 0 then return
            bin1s = bin1s[kidx]
            bininfos = bininfos[kidx,*,*]
            end
        1: begin
            kidx = where(bin2s ge slicerng[0] and bin2s le slicerng[1], nslice)
            if nslice eq 0 then return
            bin2s = bin2s[kidx]
            bininfos = bininfos[*,kidx,*]
            end
        2: begin
            kidx = where(bin3s ge slicerng[0] and bin3s le slicerng[1], nslice)
            if nslice eq 0 then return
            bin3s = bin3s[kidx]
            bininfos = bininfos[*,*,kidx]
            end
    endcase
    binsz = size(bininfos,/dimensions)
    nbin1 = binsz[0]
    nbin2 = binsz[1]
    nbin3 = binsz[2]
    bincnts = dblarr(binsz)
    bindata = ptrarr(binsz)


;---Read data and bininfo.
    for l=0, nyr-1 do begin
        restore, binfns[l]  ; bininfos.
        case sliceidx of
            0: bininfos = bininfos[kidx,*,*]
            1: bininfos = bininfos[*,kidx,*]
            2: bininfos = bininfos[*,*,kidx]
        endcase


        restore, datfns[l]
        for i=0,nbin1-1 do for j=0,nbin2-1 do for k=0,nbin3-1 do begin
            if ~ptr_valid(bininfos[i,j,k].uts) then continue
            if ptr_valid(bindata[i,j,k]) then begin
                *bindata[i,j,k] = [*bindata[i,j,k],dat0[*bininfos[i,j,k].idx]]
            endif else begin
                bindata[i,j,k] = ptr_new(dat0[*bininfos[i,j,k].idx])
            endelse
            bincnts[i,j,k] += n_elements(*bininfos[i,j,k].idx)
        endfor
    endfor
    
    ; calc the quantity to be plotted.
    tdat = dblarr(binsz)
    for i=0,nbin1-1 do for j=0,nbin2-1 do for k=0,nbin3-1 do begin
        if ~ptr_valid(bindata[i,j,k]) then continue
        case statype of
            'average': tdat[i,j,k] = mean(*bindata[i,j,k],/nan)
        endcase
    endfor


    case sliceidx of
        0: bininfos = reform(bininfos[0,*,*])
        1: bininfos = reform(bininfos[*,0,*])
        1: bininfos = reform(bininfos[*,*,0])
    endcase
    
    tsz = size(bininfos,/dimensions)
    bincnts = total(bincnts, sliceidx+1)
    tdat = total(tdat, sliceidx+1)/nslice
    tdat = alog10(tdat)
    tdat = bytscl(tdat, min=zrang[0], max=zrang[1], top=top)
    mincnt = 10
    
    
;---Plot.
;    figfn = 0
;    stop

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
        if bincnts[i,j] le mincnt then continue
    
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
    idx = where(txs lt 0)
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


vars = ['t_avg','density']
vars = ['h_'+vars,'o_'+vars]
foreach tvar, vars do fig_distr_scalar_cartesian, tvar
end
