;+
; Test to plot THEMIS E field on the vertical surface.
;-
;


;---Constants.
    secofday = 86400d


;---Settings.
    rootdir = shomedir()+'/efield_distr'
    ; time.
    date = '2008-02-02'
    tutr = time_double(date)+[0,secofday]
    stryr = time_string(tutr[0],tformat='YYYY')
    bintype = 'dis_mlt_dns'
    
    top = 254
    ezrng = [-1,1]*30
    ct0 = 70
    psym = 1
    symsz = 0.3
    scl = 0.1
    
    dis_rng = [8,12]
    mlt_rng = [-6,6]
    dns_rng = [-1,1]*4

    dis_del = 1 
    mlt_del = 1
    dns_del = 1
    
    bin1s = smkarthm(dis_rng[0],dis_rng[1],dis_del,'dx')
    bin2s = smkarthm(mlt_rng[0],mlt_rng[1],mlt_del,'dx')
    bin3s = smkarthm(dns_rng[0],dns_rng[1],dns_del,'dx')

    ; make the bins.
    nbin1 = n_elements(bin1s)-1
    nbin2 = n_elements(bin2s)-1
    nbin3 = n_elements(bin3s)-1
    binsz = [nbin1,nbin2,nbin3]
    
    bininfo0 = {$
        range1:[0d,0], $
        range2:[0d,0], $
        range3:[0d,0], $
        uts: ptr_new(), $
        idx: ptr_new(), $   ; the golabl index that maps all ut to each bin.
        dat: ptr_new(), $   ; e field.
        id:bintype}
    
    bininfos = replicate(bininfo0, nbin1, nbin2, nbin3)
    for i=0, nbin1-1 do bininfos[i,*,*].range1 = bin1s[i:i+1]
    for j=0, nbin2-1 do bininfos[*,j,*].range2 = bin2s[j:j+1]
    for k=0, nbin3-1 do bininfos[*,*,k].range3 = bin3s[k:k+1]
    
    
    ; spacecraft.
    probes = ['d','e']
    posvar0 = ['Epoch','SM_LCT_T','DNEUTS','RADIUS']
    posvar1 = ['epoch','mlt','dneuts','dis']
    efivar1 = ['utsec','egsm']
        
    foreach tprobe, probes do begin
        binfn = rootdir+'/data/thm'+tprobe+'_global_efield_'+bintype+'_'+stryr+'.sav'
        if file_test(binfn) eq 0 then begin
            dir = file_dirname(binfn)
            if file_test(dir,/directory) eq 0 then file_mkdir, dir
        endif; else continue
        
        pos = sread_thm_orbit(tutr, vars=posvar0, newname=posvar1, probes=tprobe)
        efi = sread_thm_efi_l2(tutr, newname=efivar1, probes=tprobe)
        ut1s = efi.utsec
        tuts = sfmepoch(pos.epoch,'unix')
        pos1 = sinterpol(pos.dis, tuts, ut1s)
        pos2 = sinterpol(pos.mlt, tuts, ut1s) & pos2[where(pos2 ge 12)] -= 24 & sdespike, ut1s, pos2
        pos3 = sinterpol(pos.dneuts, tuts, uts)
        dat0 = efi.egsm[*,2]
        
        for i=0, nbin1-1 do for j=0, nbin2-1 do for k=0, nbin3-1 do begin
            trng1 = minmax(bininfos[i,j,k].range1)
            trng2 = minmax(bininfos[i,j,k].range2)
            trng3 = minmax(bininfos[i,j,k].range3)

            idx = where(pos1 ge trng1[0] and pos1 lt trng1[1] $
                and pos2 ge trng2[0] and pos2 lt trng2[1] $
                and pos3 ge trng3[0] and pos3 lt trng3[1], cnt)

            if cnt eq 0 then continue
            bininfos[i,j,k].uts = ptr_new(ut1s[idx])
            bininfos[i,j,k].idx = ptr_new(idx)
            bininfos[i,j,k].dat = ptr_new(dat0[idx])
        endfor

        save, bin1s, bin2s, bin3s, bininfos, filename=binfn
        for i=0, nbin1-1 do for j=0, nbin2-1 do for k=0, nbin3-1 do begin
            ptr_free, bininfos[i,j,k].uts
            ptr_free, bininfos[i,j,k].idx
            ptr_free, bininfos[i,j,k].dat
        endfor
    endforeach


;---Generate plot. 
    datptr = ptrarr(binsz)
    foreach tprobe, probes do begin
        binfn = rootdir+'/data/thm'+tprobe+'_global_efield_'+bintype+'_'+stryr+'.sav'
        restore, filename=binfn
        for i=0, nbin1-1 do for j=0, nbin2-1 do for k=0, nbin3-1 do begin
            if ~ptr_valid(bininfos[i,j,k].dat) then continue
            if ptr_valid(datptr[i,j,k]) then *datptr[i,j,k] = [*datptr[i,j,k],*bininfos[i,j,k].dat] $
            else datptr[i,j,k] = ptr_new(*bininfos[i,j,k].dat)
        endfor
    endforeach
    
    dat0 = dblarr(binsz)
    for i=0, nbin1-1 do for j=0, nbin2-1 do for k=0, nbin3-1 do begin
        if ~ptr_valid(datptr[i,j,k]) then continue
        dat0[i,j,k] = mean(*datptr[i,j,k],/nan)
        ;dat0[i,j,k] = n_elements(*datptr[i,j,k])
    endfor
    
    dat1 = total(dat0, 1)   ; sum over distance.
    dat1 = bytscl(dat1, min=ezrng[0], max=ezrng[1], top=top, /nan)

    sgopen, ofn, xsize=6, ysize=4, /inch

    device, decomposed=0
    loadct, ct0

    pos0 = [.2,.2,.9,.9]
    
    sgtv, dat1, position=pos0, ct=ct0
    plot, mlt_rng, dns_rng, /nodata, /noerase, position=pos0, xstyle=1, ystyle=1
    
    
    foreach tprobe, probes do begin
        
        pos = sread_thm_orbit(tutr, vars=posvar0, newname=posvar1, probes=tprobe)
        efi = sread_thm_efi_l2(tutr, newname=efivar1, probes=tprobe)
        ut1s = efi.utsec
        tuts = sfmepoch(pos.epoch,'unix')
        pos1 = sinterpol(pos.dis, tuts, ut1s)
        pos2 = sinterpol(pos.mlt, tuts, ut1s) & pos2[where(pos2 ge 12)] -= 24 & sdespike, ut1s, pos2
        pos3 = sinterpol(pos.dneuts, tuts, uts)

        idx = where(pos2 ge mlt_rng[0] and pos2 le mlt_rng[1] $
            and pos3 ge dns_rng[0] and pos3 le dns_rng[1])
        pos2 = pos2[idx]
        pos3 = pos3[idx]
        oplot, pos2, pos3
    endforeach

    sgclose

end
