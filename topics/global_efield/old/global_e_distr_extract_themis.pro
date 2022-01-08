;+
;
;-


pro global_e_distr_extract_themis, utr0, var0, newname=var1, probe=tprobe

;---Constants.
    secofday = 86400d   ; sec.
    cns = -1


;---Settings.

    if n_elements(var0) eq 0 then message, 'no varname ...'
    if n_elements(var1) eq 0 then var1 = var0
    posvar0 = ['Epoch']
    posvar1 = ['epoch']
    
    ; file names.
    stryr = time_string(utr0[0],tformat='YYYY')
    
    rootdir = shomedir()+'/global_e_distr'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    
    datadir = rootdir+'/data/'+stryr
    if file_test(datadir,/directory) eq 0 then file_mkdir, datadir
        
    datfn = datadir+'/global_e_distr_th'+tprobe+'_'+stryr+'_'+var1+'.sav'
    

;---Loop through each day.
    nday = (utr0[1]-utr0[0])/secofday+1
    ut1s = smkarthm(utr0[0], secofday, nday, 'x0')
    ut2s = ut1s+secofday
    
    dat0 = []
    for i=0, nday-1 do begin
        tutr = [ut1s[i],ut2s[i]]
        printf, cns, '****Processing '+time_string(tutr[0])+' ...'
        pos = sread_thm_orbit(tutr, vars=posvar0, newname=posvar1, probes=tprobe)
        if size(pos,/type) ne 8 then begin
            printf, cns, 'no position data ...'
            continue
        endif
        if size(pos,/type) ne 8 then continue

        ;---read the efield times.
        efi = sread_thm_efi_l2(tutr, vars=var0, newname=var1, probes=tprobe)
        if size(efi,/type) ne 8 then continue           ; no data.
        if n_elements(efi.(0)) le 1 then continue       ; no real data.
        dat0 = [dat0,efi.(0)]
    endfor

    save, dat0, filename=datfn
end


stryrs = string(smkarthm(2013,2015,1,'dx'),format='(I4)')
probes = ['a','b','c','d','e']
foreach tprobe, probes do foreach tstryr, stryrs do begin
    var0 = 'th'+tprobe+'_efs_dot0_gsm'
    var1 = 'edot0_gsm'
    global_e_distr_extract_themis, tutr, var0, probe=tprobe, newname=var1
endforeach

end


