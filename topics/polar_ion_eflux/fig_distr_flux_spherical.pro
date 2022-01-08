

pro fig_distr_flux_spherical, dattype, bintypee, $
    bininfos=bininfos, tdat=tdat, bincnts=bincnts, noplot=noplot, slicerng=slicerng, _extra=extra


;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180
    strgamma = '!9'+string(71b)+'!X'
    strmu = '!9'+string(109b)+'!X'

;---Settings.
    if n_elements(dattype) eq 0 then message, 'no data type ...'
    if n_elements(bintype) eq 0 then bintype = 'mlt_mlat_dis'
    if n_elements(statype) eq 0 then statype = 'average'

    ; slice on mlt.
    sliceidx = 0
    dmlt = 2d
    if n_elements(slicerng) eq 0 then begin
        slicerng = [-1,1]*dmlt+12
    endif

    ; plot settings.
    ct0 = 70
    top = 255
    thick = 2
    
    case dattype of
        'h_energy_flux': begin
            unit = 'H!U+!N '+strgamma+'!I||!N @100km ('+strmu+'W/m!U-2!N)'
            zrang = [-1,1]*2
            end
        'o_energy_flux': begin
            unit = 'O!U+!N '+strgamma+'!I||!N @100km ('+strmu+'W/m!U-2!N)'
            zrang = [-1,1]*2
            end
        'h_number_flux': begin
            unit = 'H!U+!N '+'F!I||!N @100km (#/cm!U-2!N-s)'
            zrang = [-1,1]*6
        end
        'o_number_flux': begin
            unit = 'O!U+!N '+'F!I||!N @100km (#/cm!U-2!N-s)'
            zrang = [-1,1]*6
        end
    endcase
    
    case bintype of
        'mlt_mlat_dis': begin
            case sliceidx of
                0: begin
                    xtitl = 'X SM (Re)'
                    ytitl = 'Z SM (Re)'
                    ztitl = 'MLT in ('+sgnum2str(mean(slicerng))+','+sgnum2str(mean(slicerng)-12)+')+/-'+sgnum2str(dmlt)+' hr!C'+unit
                    xrang = [-1d,1]*7
                    yrang = [0d,1]*9
                end
            endcase
            end
    endcase
    

    nztick = 5
    ztickv = smkarthm(zrang[0],zrang[1],nztick,'n')
    ztickn = '10!U'+string(abs(ztickv),format='(I0)')
    idx = where(ztickv lt 0, cnt)
    if cnt ne 0 then ztickn[idx] = '-'+ztickn[idx]
    idx = where(ztickv eq 0, cnt)
    if cnt ne 0 then ztickn[idx] = '0'

    xsz = 6d ; inch
    ysz = abs((yrang[1]-yrang[0])/(xrang[1]-xrang[0])*xsz)

    ; filenames.
    stryrs = string(smkarthm(1996,1998,1,'dx'),format='(I4)')
    nyr = n_elements(stryrs)
    rootdir = shomedir()+'/polar_ion_eflux'
    datadir = shomedir()+'/polar_ion_eflux/data'
    binfns = datadir+'/'+stryrs+'/polar_ion_eflux_bin_info_'+bintype+'_'+stryrs+'.sav'
    datfns = datadir+'/'+stryrs+'/polar_ion_eflux_'+dattype+'_'+stryrs+'.sav'
    utsfns = datadir+'/'+stryrs+'/polar_ion_eflux_uts_'+stryrs+'.sav'
    mapfns = datadir+'/'+stryrs+'/polar_ion_eflux_mapcoef_'+stryrs+'.sav'

    figfn = shomedir()+'/fig_distr_'+dattype+'_'+bintype+'.pdf'




;---Load bin info, contains bin1s, bin2s, bin3s, bininfos.
    restore, binfns[0]
    
    ; reduce bin dimension.
    case sliceidx of
        0: begin
            bin1s[where(bin1s ge 18)] -= 24
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

        restore, mapfns[l]  ; dat0.
        mapc = dat0

        restore, datfns[l]  ; dat0.
        dat0 *= mapc        
        if strmid(dattype,10,/reverse_offset) eq 'energy_flux' then dat0 *= 1e3  ; change to uW/m^2.
        dat1 = alog10(abs(dat0))>1
        idx = where(dat0 lt 0, cnt)
        if cnt ne 0 then dat1[idx] = -dat1[idx]
        dat0 = dat1

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
            'average': tdat[i,j,k] = -mean(*bindata[i,j,k],/nan)
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
    tdat = bytscl(tdat, min=zrang[0], max=zrang[1], top=top)
    mincnt = 10
    
    if keyword_set(noplot) then return
    
    
;---Plot.
    ;figfn = 0
    ;stop

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
        if bincnts[i,j] le mincnt then continue
    
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
    
    fig_distr_flux_spherical, dattype, bintype, bininfos=bininfos, tdat=tdat, bincnts=bincnts, /noplot, slicerng=slicerng-12
    for i=0,tsz[0]-1 do for j=0,tsz[1]-1 do begin
        tc = tdat[i,j]
        if bincnts[i,j] le mincnt then continue

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
        ztickv=ztickv, ztickname=reverse(ztickn), zticks=nztick
    
    sgclose


end


vars = ['number_flux','energy_flux']
vars = ['h_'+vars,'o_'+vars]
foreach tvar, vars do fig_distr_flux_spherical, tvar
end
