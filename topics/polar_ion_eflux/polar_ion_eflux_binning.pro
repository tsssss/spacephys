;+
; Bin data based on gsm cartesian coord.
; 
; Read all the data to memory and find the index which maps all times to bin.
;-


pro polar_ion_eflux_binning, utr0, type=bintype

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
    if n_elements(bintype) eq 0 then message, 'must set bintype ...'
    case bintype of
        'gsm_xyz': begin
            x_del = 0.4d      ; in re.
            x_rng = [-7.5d,6.5]
            y_del = 0.4d      ; in re.
            y_rng = [-7.5d,6.5]
            z_del = 0.4d      ; in re.
            z_rng = [-3,9]

            bin1s = smkarthm(x_rng[0],x_rng[1],x_del,'dx')
            bin2s = smkarthm(y_rng[0],y_rng[1],y_del,'dx')
            bin3s = smkarthm(z_rng[0],z_rng[1],z_del,'dx')
            end
        'mlt_mlat_dis': begin
            mlt_del = 0.5d      ; in hr.
            mlt_rng = [0d,24]
            lat_del = 1.5d      ; in deg.
            lat_rng = [45d,90]
            dis_del = 0.2d      ; in re.
            dis_rng = [1.6d,9]

            bin1s = smkarthm(mlt_rng[0],mlt_rng[1],mlt_del,'dx')
            bin2s = smkarthm(lat_rng[0],lat_rng[1],lat_del,'dx')
            bin3s = smkarthm(dis_rng[0],dis_rng[1],dis_del,'dx')
            end
    endcase

    
    

    ; file names.
    stryr = time_string(utr0[0],tformat='YYYY')

    rootdir = shomedir()+'/polar_ion_eflux'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir

    datadir = shomedir()+'/polar_ion_eflux/data/'+stryr
    if file_test(datadir,/directory) eq 0 then file_mkdir, datadir

    logfn = rootdir+'/polar_ion_eflux_binning_'+bintype+'_'+stryr+'.log'
    if file_test(logfn) then file_delete, logfn & stouch, logfn

    binfn = datadir+'/polar_ion_eflux_bin_info_'+bintype+'_'+stryr+'.sav'
    utfn = datadir+'/polar_ion_eflux_uts_'+stryr+'.sav'
    posfn = datadir+'/polar_ion_eflux_pos_'+stryr+'_'+bintype+'.sav'


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
        id:bintype}
    
    bininfos = replicate(bininfo0, nbin1, nbin2, nbin3)
    for i=0, nbin1-1 do bininfos[i,*,*].range1 = bin1s[i:i+1]
    for j=0, nbin2-1 do bininfos[*,j,*].range2 = bin2s[j:j+1]
    for k=0, nbin3-1 do bininfos[*,*,k].range3 = bin3s[k:k+1]


    ; other settings.
    cns = -1    ; console.
    posvar0s = ['Epoch','GSM_POS','EDMLT_TIME','L_SHELL','MAG_LATITUDE']
    posvar1s = ['epoch','pos_gsm','mlt','lshell','mlat']
    var0s = 'utsec' ; the corrected ut, 'ut_sec' is the original timas ut.



;---Load all the times.
    nday = (utr0[1]-utr0[0])/secofday
    ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
    ut2s = ut1s+secofday

    load = 0
    if file_test(utfn) or file_test(posfn) then load = 1
    if keyword_set(reload) then load = 1
    
    if load then begin
        ut0s = []
        pos1 = []
        pos2 = []
        pos3 = []
        for i=0, nday-1 do begin
            tutr = [ut1s[i],ut2s[i]]
            printf, cns, '****Reading '+time_string(tutr[0])+' ...'
            pos = sread_polar_orbit(tutr, flag='def', /local_only, vars=posvar0s, newname=posvar1s)
            if size(pos,/type) ne 8 then $
                pos = sread_polar_orbit(tutr, flag='pre', /local_only, vars=posvar0s, newname=posvar1s)
            if size(pos,/type) ne 8 then begin
                printf, cns, 'no position data ...'
                continue
            endif
            if size(pos,/type) ne 8 then continue
            
            ;---read the timas times.
            tim = sread_polar_timas_sheng_moment(tutr, vars=var0s, /local_only)
            if size(tim,/type) ne 8 then continue           ; no data.
            if n_elements(tim.(0)) le 1 then continue       ; no real data.
            uts = tim.(0)
            ut0s = [ut0s,uts]
            
            ;---read position.
            case bintype of
                'gsm_xyz': begin
                    tuts = sfmepoch(pos.epoch, 'unix')
                    rgsm = pos.pos_gsm*re1  ; in-situ pos in gsm in Re.
                txs = [disrng[0]*cos(latrng[0]),disrng[1]*cos(latrng[0]),disrng[1]*cos(latrng[1]),disrng[1]*cos(latrng[1]),disrng[0]*cos(latrng[0])]
                    rgsm = float(sinterpol(rgsm, tuts, uts))
                    pos1 = [pos1,rgsm[*,0]]
                    pos2 = [pos2,rgsm[*,1]]
                    pos3 = [pos3,rgsm[*,2]]
                    end
                    
                'mlt_mlat_dis': begin
                    tuts = sfmepoch(pos.epoch, 'unix')
                    rdis = snorm(pos.pos_gsm)*re1  ; in-situ dis in Re.
                    pos1 = [pos1,float(sinterpol(pos.mlt, tuts, uts))]    ; in hr in [0,24].
                    pos2 = [pos2,float(abs(sinterpol(pos.mlat, tuts, uts)))]
                    pos3 = [pos3,float(sinterpol(rdis, tuts, uts))]
                    end
            endcase
        endfor
    
    
        ut00 = time_double(time_string(utr0[0],tformat='YYYY-01-01/00:00'))
        ut1s = ut0s-ut00
        
        if ~file_test(utfn) eq 1 then save, ut00, ut1s, filename=utfn
        if file_test(posfn) eq 1 then file_delete, posfn
        save, pos1, pos2, pos3, bintype, filename=posfn
    endif else begin
        restore, filename=utfn
        restore, filename=posfn
    endelse
    

;---Bin the data.
    openw, lun, logfn, /get_lun, /append
    printf, lun, 'Binning data in '+stryr+' according to '+bintype+' ...'
    printf, lun, ''
    printf, lun, '# of uts ('+string(n_elements(ut1s),format='(I16)')+')'
    printf, lun, '# of pos ('+string(n_elements(pos1),format='(I16)')+')'
    printf, lun, ''
    free_lun, lun

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

        openw, lun, logfn, /get_lun, /append
        cmd = 'Bin ('+$
            string(trng1[0],format='(F4.1)')+','+string(trng1[1],format='(F4.1)')+')('+$
            string(trng2[0],format='(F4.1)')+','+string(trng2[1],format='(F4.1)')+')('+$
            string(trng3[0],format='(F4.1)')+','+string(trng3[1],format='(F4.1)')+')('+$
            string(cnt,format='(I)')+')'
        printf, lun, cmd
        free_lun, lun
    endfor

    
;---save the results.
    save, bin1s, bin2s, bin3s, bininfos, filename=binfn
    for i=0, nbin1-1 do for j=0, nbin2-1 do for k=0, nbin3-1 do begin
        ptr_free, bininfos[i,j,k].uts
        ptr_free, bininfos[i,j,k].idx
    endfor
    
end

stop
bintype = 'mlt_mlat_dis'
polar_ion_eflux_binning, time_double(['1996-03-17','1996-12-31']), type=bintype
polar_ion_eflux_binning, time_double(['1997-01-01','1997-12-31']), type=bintype
polar_ion_eflux_binning, time_double(['1998-01-01','1998-12-08']), type=bintype
;stop


bintype = 'gsm_xyz'
polar_ion_eflux_binning, time_double(['1996-03-17','1996-12-31']), type=bintype
polar_ion_eflux_binning, time_double(['1997-01-01','1997-12-31']), type=bintype
polar_ion_eflux_binning, time_double(['1998-01-01','1998-12-08']), type=bintype


end
