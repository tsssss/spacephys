
pro polar_ion_eflux_extract_hydra_var, utr0, var=var0

;---Constants.
    secofday = 86400d   ; sec.
    cns = -1
    

;---Settings.
    stryr = time_string(utr0[0],tformat='YYYY')
    rootdir = shomedir()+'/polar_ion_eflux/data/'+stryr
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    datfn = rootdir+'/polar_ion_eflux_hydra_'+var0+'_'+stryr+'.sav'


    var0s = var0

;---Loop through each day.
    nday = (utr0[1]-utr0[0])/secofday
    ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
    ut2s = ut1s+secofday
    dat0 = []
    for i=0, nday-1 do begin
        tutr = [ut1s[i],ut2s[i]]
        printf, cns, '****Processing '+time_string(tutr[0])+' ...'
        pos = sread_polar_orbit(tutr, flag='def', /local_only)
        if size(pos,/type) ne 8 then pos = sread_polar_orbit(tutr, flag='pre', /local_only)
        if size(pos,/type) ne 8 then begin
            printf, cns, 'no position data ...'
            continue
        endif
        if size(pos,/type) ne 8 then continue
        
        tim = sread_polar_hydra_moment(tutr, vars=var0s, /local_only)
        if size(tim,/type) ne 8 then begin
            tim = sread_polar_hydra_moment(tutr, vars='EPOCH', /local_only)
            if size(tim,/type) ne 8 then continue
            nrec = n_elements(tim.(0))
            if nrec le 1 then continue
            dat0 = [dat0,dblarr(nrec)+!values.d_nan]
        endif else begin
            if n_elements(tim.(0)) le 1 then continue       ; no real data.
            dat0 = [dat0,tim.(0)]
        endelse
    endfor

    save, dat0, filename=datfn

end




vars = ['DENSITY','BULK_VELOCITY','TPARL','TPERP']
vars = ['BULK_VELOCITY']
vars = [vars+'_ELE',vars+'_ION']

stryrs = string(smkarthm(1996,2008,1,'dx'),format='(I4)')

foreach tvar, vars do foreach tstryr,stryrs do begin
    tutr = time_double(tstryr+'/'+['01-01','12-31'])
    polar_ion_eflux_extract_hydra_var, var=tvar, tutr
endforeach

end
