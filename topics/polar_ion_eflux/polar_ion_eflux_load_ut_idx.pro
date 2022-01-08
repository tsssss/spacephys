;+
; Read all the ut in po_tim_xxx files. Find the indices that map
; the vut to each bin.
;-
pro polar_ion_eflux_load_ut_idx, utr0, reload=reload


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
    if n_elements(utr0) eq 0 then $
        utr0 = time_double(['1996-03-17','1998-12-08'])
    var0s = 'h_'+['ut_sec']
    
    rootdir = shomedir()+'/polar_ion_eflux'
    binfns = rootdir+'/polar_ion_eflux_bin_info_cartesian_'+['1996','1997','1998']+'.sav'
    utfn = rootdir+'/polar_ion_eflux_all_ut.sav'
    utidxfn = rootdir+'/polar_ion_eflux_ut_idx.sav'
    if keyword_set(reload) then begin
        if file_test(utidxfn) eq 1 then file_delete, utidxfn
    endif
    
    cns = -1

;---Read all the ut.
    if file_test(utfn) eq 0 then begin
        ut0s = []
        nday = (utr0[1]-utr0[0])/secofday
        ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
        ut2s = ut1s+secofday

        for i=0, nday-1 do begin
            tutr = [ut1s[i],ut2s[i]]
            printf, cns, '****Processing '+time_string(tutr[0])+' ...'

            tim = sread_polar_timas(tutr, /local_only, type='h0')
            if size(tim,/type) ne 8 then continue
            if n_elements(tim.epoch_h) le 1 then continue
            uts = sfmepoch(tim.epoch_h,'unix')
            tmp = polar_ion_eflux_fix_timas_time(uts)
            stop
            ut0s = [ut0s,uts]
        endfor
        
        save, ut0s, filename=utfn
    endif else restore, filename=utfn
    

;---Load the bin info.
    restore, binfns[0]
    binsize = size(bininfos,/dimensions)
    binuts = ptrarr(binsize)
    binutidx = ptrarr(binsize)
    
    ; combine the bin uts for all years into binuts.
    foreach tfn, binfns do begin
        restore, tfn
        for i=0, binsize[0]-1 do $
            for j=0, binsize[1]-1 do $
            for k=0, binsize[2]-1 do begin
            if ~ptr_valid(bininfos[i,j,k].uts) then continue
            if ptr_valid(binuts[i,j,k]) then begin
                *binuts[i,j,k] = [*binuts[i,j,k],*bininfos[i,j,k].uts]
            endif else begin
                binuts[i,j,k] = ptr_new(*bininfos[i,j,k].uts)
            endelse
        endfor
    endforeach
    
    
    
;---Loop though each bin, find the mapping coef.
    for i=0, binsize[0]-1 do $
        for j=0, binsize[1]-1 do $
        for k=0, binsize[2]-1 do begin
            if ~ptr_valid(binuts[i,j,k]) then continue
            tuts = *binuts[i,j,k]
            nrec = n_elements(tuts)
            tidx = ulon64arr(nrec)
            for l=0, nrec-1 do begin
                tmp = where(ut0s eq tuts[l], cnt)
                if cnt ne 1 then message, 'wrong time???'
                tidx[l] = tmp[0]
            endfor
            binutidx[i,j,k] = ptr_new(tidx)
        endfor
    save, binutidx, filename=utidxfn

end


polar_ion_eflux_load_ut_idx, /reload
end