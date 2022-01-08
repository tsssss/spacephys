
deg = 180d/!dpi
rad = !dpi/180


rootdir = shomedir()+'/polar_ion_eflux'
fns = rootdir+'/polar_ion_eflux_bin_info_cartesian_'+['1996','1997','1998']+'.sav'

restore, fns[0]
binsize = size(bininfos,/dimensions)
bincnts = dblarr(binsize)
bindens = dblarr(binsize)

foreach tfn, fns do begin
    restore, tfn
    for i=0,binsize[0]-1 do $
        for j=0,binsize[1]-1 do $
            for k=0, binsize[2]-1 do begin
                if ~ptr_valid(bininfos[i,j,k].uts) then continue
                bincnts[i,j,k] += n_elements(*bininfos[i,j,k].uts)
                bindens[i,j,k] += total(*bininfos[i,j,k].nhs)
            endfor

endforeach



    x_del = 0.4d      ; in re.
    x_rng = [-7.5d,6.5]
    y_del = 0.4d      ; in re.
    y_rng = [-7.5d,6.5]
    z_del = 0.4d      ; in re.
    z_rng = [-3,9]
    
    binxs = smkarthm(x_rng[0],x_rng[1],x_del,'dx')
    binys = smkarthm(y_rng[0],y_rng[1],y_del,'dx')
    binzs = smkarthm(z_rng[0],z_rng[1],z_del,'dx')


ct0 = 40
top = 255
thick = 2

xxrng = [-7.5,6.5]
yyrng = [-7.5,6.5]
zzrng = [-3,9]

zrng = [-2,2]
ztickv = smkarthm(zrng[0],zrng[1],1,'dx')
ztickn = '10!U'+string(ztickv,format='(I0)')
ztitl = 'Density (cm!U-3!N)!C '

minlat = 55

slicerngs = smkarthm(yyrng[0],yyrng[1],1,'dx')
slicerngs = [-1,1]*1.5
nslicerng = n_elements(slicerngs)
for k=0, nslicerng-2 do begin
    kidx = where(binys ge slicerngs[k] and binys le slicerngs[k+1], cnt)
    if cnt eq 0 then continue
    
    
    ofn = shomedir()+'/fig_data_coverage_y_'+string(slicerngs[k],format='(F4.1)')+'_'+string(slicerngs[k+1],format='(F4.1)')+'_re.pdf'
    ;ofn = 0
    ;stop

    sgopen, ofn, xsize=6, ysize=5, /inch
    
    device, decomposed=0
    loadct, ct0
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    pos0 = [0.15,0.1,0.9,0.85]
    pos1 = pos0 & pos1[1] = pos0[3]+ychsz*1 & pos1[3] = pos1[1]+ychsz*0.5
    
    
    plot, xxrng, zzrng, position=pos0, $
        xstyle=1, xrange=xxrng, xtitle='X GSM (Re)', $
        ystyle=1, yrange=zzrng, ytitle='Z GSM (Re)', $
        /nodata, /noerase, /isotropic
    
    
        ;kidx = where(binzzs ge 3 and binzzs le 5)
        tbincnts = total(bincnts[*,kidx,*],2)
        tbindens = total(bindens[*,kidx,*],2)
        
        tdat = tbindens/tbincnts
        ;tdat = tbinefxs/tbincnts*slicerngs[k]^3
        tdat = alog10(tdat)
        tdat = bytscl(tdat, min=zrng[0], max=zrng[1], top=top)
        
        for i=0,binsize[0]-1 do $
            for j=0,binsize[2]-1 do begin
                txxrng = bininfos[i,kidx[0],j].x_rng
                tyyrng = bininfos[i,kidx[0],j].z_rng
                txs = [txxrng[0],txxrng[0],txxrng[1],txxrng[1],txxrng[0]]
                tys = [tyyrng[0],tyyrng[1],tyyrng[1],tyyrng[0],tyyrng[0]]
                tc = tdat[i,j]
                polyfill, txs, tys, /data, color=tc
            endfor
            
        device, decompose=1
        tmp = findgen(51)/50*2*!dpi
        txs = cos(tmp)
        tys = sin(tmp)
        plots, txs, tys, /data, color=sgcolor('white'), thick=thick
                    
        sgcolorbar, findgen(top), /horizontal, position=pos1, /noerase, ct=ct0, $
            zrange=zrng, ztitle='Y GSM in '+sgnum2str(slicerngs[k])+' to '+sgnum2str(slicerngs[k+1])+' Re!C'+ztitl, $
            ztickv=ztickv, ztickname=ztickn
        
    sgclose

endfor

end
