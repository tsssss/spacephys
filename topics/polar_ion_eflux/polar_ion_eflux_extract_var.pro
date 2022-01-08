
pro polar_ion_eflux_extract_var, utr0, var=var0

;---Constants.
    secofday = 86400d   ; sec.
    cns = -1
    

;---Settings.
    stryr = time_string(utr0[0],tformat='YYYY')
    rootdir = shomedir()+'/polar_ion_eflux/data/'+stryr
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    datfn = rootdir+'/polar_ion_eflux_'+var0+'_'+stryr+'.sav'


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
        
        tim = sread_polar_timas_sheng_moment(tutr, vars=var0s, /local_only)
        if size(tim,/type) ne 8 then continue           ; no data.
        if n_elements(tim.(0)) le 1 then continue       ; no real data.
        
        dat0 = [dat0,tim.(0)]
    endfor

    save, dat0, filename=datfn

end


utr0 = time_double([['1996-03-17','1996-12-31'],$
    ['1997-01-01','1997-12-31'],$
    ['1998-01-01','1998-12-08']])
nutr = n_elements(utr0)/2
vars = ['density','bulk_velocity','t_avg','energy_flux']
vars = ['h_'+vars,'o_'+vars]
vars = ['mapcoef']
vars = ['h_','o_']+'number_flux'
foreach tvar, vars do for i=0, nutr-1 do $
    polar_ion_eflux_extract_hydra_var, var=tvar, utr0[*,i]

end
