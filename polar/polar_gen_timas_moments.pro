;+
; On daily basis, read Timas h0 data, calc moments, and save cdf.
;-

pro polar_gen_timas_moments, utr

;---Settings.
    strmu = 'u'
    m0 = 1.67d-27   ; proton mass in kg.
    species = ['h','o','he_1','he_2']
    nspecie = n_elements(species)
    epvars = 'epoch_'+species
    ionvars = 'flux_'+species
    
    m_es = m0/1.6d-19*1d6*[1,16,4,4]
    q0s = [1d,1,1,2]
    vsc = 0d
    
        
    rootdir = sdiskdir('Research')+'/sdata/polar/timas'
    datdir = rootdir+'/'+time_string(utr[0],tformat='YYYY')
    if file_test(datdir,/directory) eq 0 then file_mkdir, datdir
    cdffn = datdir+'/po_tim_moments_'+time_string(utr[0],tformat='YYYY_MMDD')+'.cdf'
    
    
    logfn = shomedir()+'/polar_gen_timas_moments.log'
    if file_test(logfn) eq 0 then stouch, logfn
    
    
;---Read Timas h0 data.
    date = time_string(utr[0],tformat='YYYY_MMDD')
    dat = sread_polar_timas(utr, type='h0')
    if size(dat,/type) ne 8 then begin
        openw, lun, logfn, /get_lun, /append
        printf, lun, ''
        printf, lun, 'No data on '+date+' ...'
        free_lun, lun
        return
    endif
    if n_elements(dat.epoch_h) le 1 then begin
        openw, lun, logfn, /get_lun, /append
        printf, lun, ''
        printf, lun, 'Little data on '+date+' ...'
        free_lun, lun
        return
    endif
    
    printf, -1, '**** Processing '+date+' ...'
    tags = strlowcase(tag_names(dat))
    
    
    ; save uts and gatt.
    if file_test(cdffn) eq 1 then file_delete, cdffn
    utname = 'ut_sec'
    ginfo = {$
        title: 'Polar TIMAS moments, calculated based on h0 data',$
        text: 'Calculated by Sheng Tian at the University of Minnesota'}
    if file_test(cdffn) eq 0 then scdfwrite, cdffn, gattribute=ginfo

    ; ut time.
    tval =  sfmepoch(dat.epoch_h, 'unix')
    ainfo = {$
        FIELDNAM:'Center UT for the moments',$
        UNITS:'s (UT)',$
        VAR_TYPE:'support_data'}
    scdfwrite, cdffn, utname, value=tval, cdftype='CDF_DOUBLE', attribute=ainfo
    
    
    ; ut time 2, correct irregularities in the original time.
    tval =  polar_ion_eflux_fix_timas_time(tval)
    ainfo = {$
        FIELDNAM:'Center UT for the moments, corrected by polar_ion_eflux_fix_timas_time',$
        UNITS:'s (UT)',$
        VAR_TYPE:'support_data'}
    scdfwrite, cdffn, 'utsec', value=tval, cdftype='CDF_DOUBLE', attribute=ainfo


;---Loop through each species.
    for i=0, nspecie-1 do begin

        uts = sfmepoch(dat.(where(tags eq epvars[i])), 'unix')
        dr0 = sdatarate(uts)
        nrec = n_elements(uts)

        fdat = dat.(where(tags eq ionvars[i]))
        sz0 = size(fdat,/dimensions)
        npol = sz0[2]
        nsec = 30d
        nen0 = sz0[1]


        theta = 90d -dat.angle      ; pitch angle, in deg.
        dtheta = 180d/npol
        dphi = 360d/nsec
        phi = findgen(nsec)*dphi

        en0s = dat.energy            ; energy bins, in eV.
        de_e = abs(shift(en0s,-1)-shift(en0s,1))*0.5/en0s
        de_e[0] = de_e[nen0-2]
        de_e[nen0-1] = de_e[nen0-2]
        den0 = en0s*de_e


        moms = replicate(smom3d_init(), nrec)
        for j=0, nrec-1 do begin
            tmp = reform(fdat[j,*,*])
            idx = where(finite(tmp,/nan),cnt)
            tmp[idx] = 0d
            tfdat = dblarr(nen0,npol,nsec)
            for k=0, nsec-1 do tfdat[*,*,k] = tmp

            ttheta = dblarr(nen0,npol,nsec)
            for k=0, npol-1 do ttheta[*,k,*] = theta[k]
            tdtheta = dblarr(nen0,npol,nsec)+dtheta

            tphi = dblarr(nen0,npol,nsec)
            for k=0, nsec-1 do tphi[*,*,k] = phi[k]
            tdphi = dblarr(nen0,npol,nsec)+dphi

            ten0s = dblarr(nen0,npol,nsec)
            for k=0, nen0-1 do ten0s[k,*,*] = en0s[k]
            den0s = dblarr(nen0,npol,nsec)
            for k=0, nen0-1 do den0s[k,*,*] = den0[k]

            tmp = nsec*npol
            tfdat = reform(tfdat, nen0,tmp)
            ten0s = reform(ten0s, nen0,tmp)
            den0s = reform(den0s, nen0,tmp)
            tphi = reform(tphi, nen0,tmp)
            tdphi = reform(tdphi, nen0,tmp)
            ttheta = reform(ttheta, nen0,tmp)
            tdtheta = reform(tdtheta, nen0,tmp)

            dat0 = smom3d_init(tfdat)
            dat0.time = uts[j]
            dat0.end_time = uts[j]+dr0
            dat0.charge = q0s[i]
            dat0.mass = m_es[i]
            dat0.species = species[i]
            dat0.magf = [0d,0,1]
            dat0.sc_pot = vsc
            dat0.units = 'eflux'
            dat0.energy = ten0s
            dat0.denergy = den0s
            dat0.phi = tphi
            dat0.dphi = tdphi
            dat0.theta = ttheta
            dat0.dtheta = tdtheta
            dat0.nenergy = nen0
            dat0.bins = 1d
            dat0.valid = 1d

            moms[j] = smom3d(dat0)
            ;tmom = moments_3d(dat0,/no_unit_conversion)
        endfor
        
    ;---save data.
        tspecie = strupcase(species[i])
        pre0 = strlowcase(tspecie)+'_'
        
        ; density.
        vname = pre0+'density'
        tval = moms.density
        ainfo = {$
            FIELDNAM:'Density '+tspecie,$
            UNITS:'cm!U-3!N',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
        
        ; average temperature.
        vname = pre0+'t_avg'
        tval = moms.tavg
        ainfo = {$
            FIELDNAM:'Average temperature for '+tspecie,$
            UNITS:'eV',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
        
        ; number flux.
        vname = pre0+'number_flux'
        tval = moms.nflux[2]
        ainfo = {$
            FIELDNAM:'Number flux for '+tspecie,$
            UNITS:'#/(cm!U2!N-s)',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
        
        ; energy flux.
        vname = pre0+'energy_flux'
        tval = moms.eflux[2]
        ainfo = {$
            FIELDNAM:'Energy flux for '+tspecie,$
            UNITS:'mW/m!U2!N',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
        
        ; enthalpy.
        vname = pre0+'enthalpy'
        tval = moms.enthpy
        ainfo = {$
            FIELDNAM:'Enthalpy for '+tspecie,$
            UNITS:'mW/m!U2!N',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
                
        ; bulk velocity.
        vname = pre0+'bulk_velocity'
        tval = moms.vbulk[2]
        ainfo = {$
            FIELDNAM:'Bulk velocity for '+tspecie,$
            UNITS:'km/s',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
        
        ; pressure tensor.
        vname = pre0+'p_tensor'
        tval = moms.ptens
        ainfo = {$
            FIELDNAM:'Pressure tensor for '+tspecie,$
            UNITS:'eV/cm!U-3!N',$
            VAR_TYPE:'data',$
            DEPEND_0:utname}
        scdfwrite, cdffn, vname, value=tval, attribute=ainfo
    endfor
end


utr0 = time_double(['1996-03-17','1998-12-08'])
;utr0 = time_double(['1996-09-04','1998-12-08'])
secofday = 86400d

nday = (utr0[1]-utr0[0])/secofday
ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
ut2s = ut1s+secofday

stop
for i=0, nday-1 do begin
    utr = [ut1s[i],ut2s[i]]
    polar_gen_timas_moments, utr
endfor

end