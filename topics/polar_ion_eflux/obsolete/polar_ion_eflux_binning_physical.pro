;+
; Bin data based on physical coord: distance, mlt, ilat, and dis.
;-


pro polar_ion_eflux_binning_physical, utr0

;---Constants.
    secofday = 86400d   ; sec.
    secofhour = 3600d   ; sec.

    re = 6378d          ; km.
    re1 = 1d/re         ; km^-1.
    h0 = 100d           ; km.
    r0 = h0*re1+1       ; Re.

    deg = 180d/!dpi
    rad = !dpi/180
    


;---Settings.

    ; bin settings.
    cns = -1    ; console.
    ; file names.
    rootdir = shomedir()+'/polar_ion_eflux'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    logfn = rootdir+'/polar_ion_eflux_binning_physical_'+time_string(utr0[0],tformat='YYYY')+'.log'
    if file_test(logfn) then file_delete, logfn & stouch, logfn
    datfn = rootdir+'/polar_ion_eflux_bin_info_physical_'+$
        time_string(utr0[0],tformat='YYYY')+'.sav'

    
    mlt_del = 0.5d      ; in hr.
    mlt_rng = [0d,24]
    lat_del = 1d        ; in deg.
    lat_rng = [55d,90]
    dis_del = 0.2d      ; in re.
    dis_rng = [1.6d,9]
    
    binmlts = smkarthm(mlt_rng[0],mlt_rng[1],mlt_del,'dx')
    binlats = smkarthm(lat_rng[0],lat_rng[1],lat_del,'dx')
    bindiss = smkarthm(dis_rng[0],dis_rng[1],dis_del,'dx')
    
    nbinmlt = n_elements(binmlts)-1
    nbinlat = n_elements(binlats)-1
    nbindis = n_elements(bindiss)-1
    
    
    bininfo0 = {$
        mlt_rng:[0d,0], $
        lat_rng:[0d,0], $
        dis_rng:[0d,0], $
        uts: ptr_new(), $
        idx: ptr_new(), $   ; the golabl index that maps all ut to each bin.
        id:0d}
    
    bininfos = replicate(bininfo0, nbinmlt, nbinlat, nbindis)
    for i=0, nbinlat-1 do bininfos[*,i,*].lat_rng = binlats[i:i+1]
    for i=0, nbinmlt-1 do bininfos[i,*,*].mlt_rng = binmlts[i:i+1]
    for i=0, nbindis-1 do bininfos[*,*,i].dis_rng = bindiss[i:i+1]

    ; mapping settings.
    par = 2
    t89 = 1
    t96 = 0
    t01 = 0
    ts04 = 0
    storm = 0



;---Loop through each day.
    nday = (utr0[1]-utr0[0])/secofday
    ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
    ut2s = ut1s+secofday

    var0s = 'Epoch_H'
    offset = 0ull           ; offset for the global index.
    
    for i=0, nday-1 do begin
        tutr = [ut1s[i],ut2s[i]]
        printf, cns, '****Processing '+time_string(tutr[0])+' ...'
        openw, lun, logfn, /get_lun, /append
        printf, lun, ''
        printf, lun, '****Processing '+time_string(tutr[0],tformat='YYYY_MMDD')+' ...'
        free_lun, lun
        pos = sread_polar_orbit(tutr, flag='def', /local_only)
        if size(pos,/type) ne 8 then pos = sread_polar_orbit(tutr, flag='pre', /local_only)
        if size(pos,/type) ne 8 then begin
            printf, cns, 'no position data ...'
            continue
        endif
        if size(pos,/type) ne 8 then continue
        
        tim = sread_polar_timas(tutr, vars=var0s, type='h0', /local_only)
        if size(tim,/type) ne 8 then continue           ; no data.
        if n_elements(tim.epoch_h) le 1 then continue   ; no real data.
        ets = tim.epoch_h
        uts = sfmepoch(ets,'unix')
        uts = polar_ion_eflux_fix_timas_time(uts, error=err)
        nrec = n_elements(uts)

        openw, lun, logfn, /get_lun, /append
        printf, lun, 'Original time rate: '+sgnum2str(sdatarate(ets)*1e-3,ndec=1)+' sec'
        printf, lun, 'Modified time rate: '+sgnum2str(sdatarate(uts),ndec=1)+' sec, error='+sgnum2str(err,ndec=1)+' sec'
        free_lun, lun
        
        ;---read position, and map it to the ionosphere.
            tuts = sfmepoch(pos.epoch, 'unix')
            rdis = snorm(pos.pos_gse)*re1  ; in-situ dis in Re.
            rdis = sinterpol(rdis, tuts, uts)
            rmlt = sinterpol(pos.mlt, tuts, uts)    ; in hr in [0,24].
            ilat = acos(sqrt(1d/pos.lshell))*deg    ; in deg.
            ilat = sinterpol(ilat, tuts, uts)
            


        ;---binning.
            for j=0, nbinmlt-1 do begin
                for k=0, nbinlat-1 do begin
                    for l=0, nbindis-1 do begin
                        tdisrng = bininfos[j,k,l].dis_rng
                        tmltrng = bininfos[j,k,l].mlt_rng
                        tlatrng = bininfos[j,k,l].lat_rng
                        idx = where(rdis ge tdisrng[0] and rdis lt tdisrng[1] $
                            and ilat ge tlatrng[0] and ilat lt tlatrng[1] $
                            and rmlt ge tmltrng[0] and rmlt lt tmltrng[1], cnt)
                        if cnt eq 0 then continue
                        if ptr_valid(bininfos[j,k,l].uts) then logflag = 1 else logflag = 0
                        if logflag eq 1 then begin
                            openw, lun, logfn, /get_lun, /append
                            cmd = 'Bin: ('+$
                                string(tdisrng[0],format='(F4.1)')+','+string(tdisrng[1],format='(F4.1)')+')('+$
                                string(tmltrng[0],format='(F4.1)')+','+string(tmltrng[1],format='(F4.1)')+')('+$
                                string(tlatrng[0],format='(F4.1)')+','+string(tlatrng[1],format='(F4.1)')+')('+$
                                string(n_elements(*bininfos[j,k,l].uts),format='(I8)')+','
                        endif
                        if ptr_valid(bininfos[j,k,l].uts) then begin
                            *bininfos[j,k,l].uts = [*bininfos[j,k,l].uts,uts[idx]]
                            *bininfos[j,k,l].idx = [*bininfos[j,k,l].idx,idx+offset]
                        endif else begin
                            bininfos[j,k,l].uts = ptr_new(uts[idx])
                            bininfos[j,k,l].idx = ptr_new(idx+offset)
                        endelse
                        if logflag eq 1 then begin
                            cmd = cmd+string(n_elements(*bininfos[j,k,l].uts),format='(I8)')+')'
                            printf, lun, cmd
                            free_lun, lun
                        endif
                    endfor
                endfor
            endfor

            offset += nrec
    endfor

    
;---save the results.
    save, binmlts, binlats, bindiss, bininfos, filename=datfn
    for i=0, nbinmlt-1 do for j=0, nbinlat-1 do for k=0, nbindis-1 do begin
        ptr_free, bininfos[i,j,k].uts
        ptr_free, bininfos[i,j,k].idx
    endfor
    
end

stop
polar_ion_eflux_binning_physical, time_double(['1996-03-17','1996-12-31'])
polar_ion_eflux_binning_physical, time_double(['1997-01-01','1997-12-31'])
polar_ion_eflux_binning_physical, time_double(['1998-01-01','1998-12-08'])
end
