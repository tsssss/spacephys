;+
; Bin THEMIS E field on the vertical surface.
;-


pro global_e_distr_bin_themis, utr0, type=bintype, reload=reload, probe=tprobe

;---Constants.
    secofday = 86400d   ; sec.


;---Settings.

    ; bin settings.
    if n_elements(bintype) eq 0 then bintype = 'dis_mlt_dns'
    case bintype of
        'dis_mlt_dns': begin
            dis_del = 1 
            dis_rng = [6,20]
            mlt_del = 1
            mlt_rng = [-6,6]
            dns_del = 1
            dns_rng = [-1,1]*4

            bin1s = smkarthm(dis_rng[0],dis_rng[1],dis_del,'dx')
            bin2s = smkarthm(mlt_rng[0],mlt_rng[1],mlt_del,'dx')
            bin3s = smkarthm(dns_rng[0],dns_rng[1],dns_del,'dx')
            end
    endcase

    ; file names.
    stryr = time_string(utr0[0],tformat='YYYY')

    rootdir = shomedir()+'/global_e_distr'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir

    datadir = rootdir+'/data/'+stryr
    if file_test(datadir,/directory) eq 0 then file_mkdir, datadir

    logfn = rootdir+'/global_e_distr_bin_th'+tprobe+'_'+bintype+'_'+stryr+'.log'
    if file_test(logfn) then file_delete, logfn & stouch, logfn

    binfn = datadir+'/global_e_distr_bin_th'+tprobe+'_'+bintype+'_'+stryr+'.sav'
    utfn = datadir+ '/global_e_distr_bin_th'+tprobe+'_'+stryr+'.sav'
    posfn = datadir+'/global_e_distr_bin_th'+tprobe+'_'+stryr+'_'+bintype+'.sav'
    

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
    posvar0 = ['Epoch','SM_LCT_T','DNEUTS','RADIUS']
    posvar1 = ['epoch','mlt','dneuts','dis']


;---Load all the times.
    nday = (utr0[1]-utr0[0])/secofday+1
    ut1s = smkarthm(utr0[0], secofday, nday, 'x0')
    ut2s = ut1s+secofday

    load = 0
    if file_test(utfn) eq 0 or file_test(posfn) eq 0 then load = 1
    if keyword_set(reload) then load = 1
    
    if load then begin
        ut0s = []
        pos1 = []
        pos2 = []
        pos3 = []
        for i=0, nday-1 do begin
            tutr = [ut1s[i],ut2s[i]]
            printf, cns, '****Reading '+time_string(tutr[0])+' ...'
            pos = sread_thm_orbit(tutr, vars=posvar0, newname=posvar1, probes=tprobe)
            if size(pos,/type) ne 8 then begin
                printf, cns, 'no position data ...'
                continue
            endif
            if size(pos,/type) ne 8 then continue

            ;---read the efield times.
            efi = sread_thm_efi_l2(tutr, probes=tprobe)
            if size(efi,/type) ne 8 then continue           ; no data.
            if n_elements(efi.(0)) le 1 then continue       ; no real data.
            uts = efi.(0)
            ut0s = [ut0s,uts]

            ;---read position.
            case bintype of
                'dis_mlt_dns': begin
                    tuts = sfmepoch(pos.epoch, 'unix')
                    dis = sinterpol(pos.dis, tuts, uts)
                    tmp = pos.mlt & tmp[where(pos.mlt ge 12)] -= 24
                    mlt = sinterpol(tmp, tuts, uts)
                    dneut = sinterpol(pos.dneuts, tuts, uts)
                    pos1 = [pos1,dis]
                    pos2 = [pos2,mlt]
                    pos3 = [pos3,dneut]
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
        openw, lun, logfn, /get_lun, /append
        cmd = 'Bin ('+$
            string(trng1[0],format='(F4.1)')+','+string(trng1[1],format='(F4.1)')+')('+$
            string(trng2[0],format='(F4.1)')+','+string(trng2[1],format='(F4.1)')+')('+$
            string(trng3[0],format='(F4.1)')+','+string(trng3[1],format='(F4.1)')+')('+$
            string(cnt,format='(I)')+')'
        printf, lun, cmd
        free_lun, lun
        
        if cnt eq 0 then continue
        bininfos[i,j,k].uts = ptr_new(ut1s[idx])
        bininfos[i,j,k].idx = ptr_new(idx)
    endfor
    
;---save the results.
    save, bin1s, bin2s, bin3s, bininfos, filename=binfn
    for i=0, nbin1-1 do for j=0, nbin2-1 do for k=0, nbin3-1 do begin
        ptr_free, bininfos[i,j,k].uts
        ptr_free, bininfos[i,j,k].idx
    endfor
    
end


stryrs = string(smkarthm(2013,2015,1,'dx'),format='(I4)')
stryrs = string(smkarthm(2008,2012,1,'dx'),format='(I4)')
stryrs = string(smkarthm(2012,2013,1,'dx'),format='(I4)')
probes = ['a','b','c','d','e']
foreach tprobe, probes do foreach tstryr, stryrs do begin
    tutr = time_double(tstryr+'/'+['01-01','12-31'])
    global_e_distr_bin_themis, tutr, probe=tprobe, /reload
endforeach

end