;+
; Read RBSP L2 data, calc moments, and save to cdf.
;
; time. Start and end times in UT sec.
; probe. A string 'a' or 'b'.
; ion_energy_range. The lower and upper limits of the ion energy bins.
; electron_energy_range. The lower and upper limits of the electron energy bins.
; filename=. A string of the CDF file, to which data are saved to.
; balance_pixel=. A boolean to balance opposite pixels, for better odd order fluxes like the energy flux.
; ncycle=. A number sets the number of cycles to be combined for a better count rate. It's 1 by default.
;-

pro rbsp_gen_hope_moments, time, probe=probe, errmsg=errmsg, $
    ion_energy_range=ion_erng, electron_energy_range=ele_erng, $
    oxygen_energy_range=o_erng, helium_energy_range=he_erng, $
    release=release, $
    filename=cdf_file, balance_pixel=balance_pixel, ncycle=ncycle


;---Check inputs.
    if n_elements(probe) ne 1 then begin
        errmsg = handle_error('Invalid probe ...')
        return
    endif

    if n_elements(ion_erng) ne 2 then ion_erng = [30,1e6]
    if n_elements(ele_erng) ne 2 then ele_erng = [200,1e6]
    if n_elements(o_erng) ne 2 then o_erng = ion_erng
    if n_elements(he_erng) ne 2 then he_erng = ion_erng
    if n_elements(ncycle) eq 0 then ncycle = 1

    test_cristian = 0

;---Constants.
    errmsg = ''
    re = 6378d & re1 = 1d/re
    p_mass = 1.67d-27   ; kg.
    e_mass = 0.91d-30   ; kg.
    xyz = ['x','y','z']
    uvw = ['u','v','w']
    fac = ['b','w','n']
    rgb = sgcolor(['red','green','blue'])
    species_info = dictionary()
    species_info['e']  = dictionary('q',-1, 'mass',e_mass, 'energy_range'  ,ele_erng)
    species_info['p']  = dictionary('q', 1, 'mass',p_mass, 'energy_range'  ,ion_erng)
    species_info['o']  = dictionary('q', 1, 'mass',p_mass*16,'energy_range',o_erng)
    species_info['he'] = dictionary('q', 1, 'mass',p_mass*4,'energy_range' ,he_erng)
    species = species_info.keys()
    ion_species = species[where(species ne 'e')]
    dt_phase_ratio = 180d/360     ; in terms of 360 deg.

    strmu = '!9'+string(109b)+'!X'
    pre0 = 'rbsp'+probe+'_'

;---Load Hope L2 data.
    foreach id, 'l2%'+['electron','ion'] do begin
        rbsp_read_hope, time, id=id, probe=probe, errmsg=errmsg, release=release
        if errmsg ne '' then return
    endforeach

;---Load quaternion data.
    rbsp_read_quaternion, time, probe=probe, errmsg=errmsg
    if errmsg ne '' then return


;---loop through each species.
    foreach the_species, species do begin
        the_type = (the_species eq 'e')? 'ele': 'ion'

    ;---Get the flux data and the energy bins.
        ; fdat, ut0s, ut1s.
        flux_var = 'f'+the_species+'du'
        get_data, flux_var, times, fdat
        fdat = fdat*1e-3    ; from 1/s-cm^2-sr-keV to 1/s-cm^2-sr-eV.
        nrec = n_elements(times)
        dims = size(reform(fdat[0,*,*,*]),/dimensions)
        nen0 = dims[0]
        nsec = dims[1]
        npol = dims[2]

        en_var = 'hope_energy_'+the_type
        en0s = get_var_data(en_var)     ; in eV.


    ;---Balancing the flux for oppsite pixels.
        if keyword_set(balance_pixel) then begin
            ; H and He will use that of H.
            if the_species eq 'p' or the_species eq 'e' then begin
                ; ratio of efficiency of the opposite pxiel.
                reff = dblarr(nrec,nen0,nsec,(npol+1)/2)

                idx1 = findgen(nsec)
                idx2 = shift(idx1,nsec/2)
                coef1 = dblarr(nen0,(npol+1)/2)
                coef2 = dblarr(nen0,npol)       ; weight for each pixel.

                for jj=0, nrec-1 do begin
                    ; flux, energy and denergy.
                    tfdat = reform(fdat[jj,*,*,*])
                    ten0s = reform(en0s[jj,*])
                    for k=0, nen0-1 do tfdat[k,*,*] *= ten0s[k]

                    for k=0, npol/2 do begin
                        v1 = tfdat[*,idx1,k]
                        v2 = tfdat[*,idx2,npol-1-k]
                        reff[jj,*,*,k] = (v1-v2)/(v1+v2)
                    endfor
                endfor

                for jj=0, nen0-1 do begin
                    for k=0, npol/2 do $
                        coef1[jj,k] = mean(reff[*,jj,*,k],/nan)
                    idx = where(finite(coef1,/nan),cnt)
                    if cnt ne 0 then coef1[idx] = 0

                    coef2[jj,4] = 1d/(1-coef1[jj,0])
                    coef2[jj,0] = 1d/(1+coef1[jj,0])
                    coef2[jj,3] = 1d/(1-coef1[jj,1])
                    coef2[jj,1] = 1d/(1+coef1[jj,1])
                    coef2[jj,2] = 1d
                endfor
                store_data, pre0+the_species+'_coef', coef1, coef2
            endif else get_data, pre0+'p_coef', coef1, coef2    ; use h.

            ; Apply pixel balancing.
            for jj=0, nen0-1 do for k=0, npol-1 do fdat[*,jj,*,k] *= coef2[jj,k]
        endif


        ; Get the start/end times of an integration cycle.
        dt_var = 'epoch_'+the_type+'_delta'
        dtimes = get_var_data(dt_var)*2d-3
        ut0s = times-dtimes*(dt_phase_ratio)
        ut1s = times+dtimes*(1-dt_phase_ratio)

        ; Combine integration cycles.
        if ncycle gt 1 then begin
            nrec = floor(nrec/ncycle)
            for jj=0, nrec-1 do begin
                j1 = jj*ncycle
                j2 = j1+ncycle-1
                fdat[jj,*,*,*] = total(fdat[j1:j2,*,*,*],1)/ncycle
                en0s[jj,*] = total(en0s[jj,*],1)/ncycle
                ut0s[jj] = ut0s[j1]
                ut1s[jj] = ut1s[j2]
            endfor
            fdat = fdat[0:nrec-1,*,*,*]
            en0s = en0s[0:nrec-1,*]
            ut0s = ut0s[0:nrec-1]
            ut1s = ut1s[0:nrec-1]
        endif

    ;---loop through each integration cycle to calculate the moments.
    ; need fdat, en0s, ut0s, ut1s, nrec.
        moms = replicate(smom3d_init(), nrec)
        q0 = species_info[the_species].q
        m0 = species_info[the_species].mass
        m_e = m0/1.6d-19*1d6
        vsc = 0d
        for jj=0, nrec-1 do begin
            ; flux, energy and denergy.
            tfdat = reform(fdat[jj,*,*,*])
            ten0s = reform(en0s[jj,*])
            for k=0, nen0-1 do tfdat[k,*,*] *= ten0s[k]

            de_e = abs(shift(ten0s,1)-shift(ten0s,-1))*0.5/ten0s
            de_e[[0,nen0-1]] = de_e[nen0-2]
            den0s = de_e*ten0s

            tmp = dblarr(nen0,nsec,npol)
            for k=0, nen0-1 do tmp[k,*,*] = ten0s[k]
            ten0s = tmp
            tmp = dblarr(nen0,nsec,npol)
            for k=0, nen0-1 do tmp[k,*,*] = den0s[k]
            den0s = tmp

            dphi = 360d/nsec
            dephi = dphi/nen0
            phi = dblarr(nen0,nsec,npol)
            for k=0, nsec-1 do phi[0,k,*] = k*dphi
            for k=1, nen0-1 do phi[k,*,*] = phi[0,*,*]+k*dephi

            dtheta = 180d/npol
            theta = dblarr(nen0,nsec,npol)
            for k=0, npol-1 do theta[*,*,k] = k*dtheta
            theta = theta-90    ; [-90,90], latitude.

            dphi = dphi+dblarr(nen0,nsec,npol)
            dtheta = dtheta+dblarr(nen0,nsec,npol)

            ; convert flux, energy, denergy, phi, dphi, theta, dtheta to 2d.
            tmp = nsec*npol
            tfdat = reform(tfdat, nen0,tmp)
            ten0s = reform(ten0s, nen0,tmp)
            den0s = reform(den0s, nen0,tmp)
            phi = reform(phi, nen0,tmp)
            dphi = reform(dphi, nen0,tmp)
            theta = reform(theta, nen0,tmp)
            dtheta = reform(dtheta, nen0,tmp)


            dat0 = smom3d_init(tfdat)
            dat0.time = ut0s[jj]
            dat0.end_time = ut1s[jj]
            dat0.charge = q0
            dat0.mass = m_e
            dat0.species = the_species
            dat0.magf = [0d,0,0]
            dat0.sc_pot = vsc
            dat0.units = 'eflux'
            dat0.energy = ten0s
            dat0.denergy = den0s
            dat0.phi = phi
            dat0.dphi = dphi
            dat0.theta = theta
            dat0.dtheta = dtheta
            dat0.nenergy = nen0
            dat0.bins = 1d
            dat0.valid = 1d

            erng = species_info[the_species].energy_range
            moms[jj] = smom3d(dat0, erange=erng)
            ;tmom = moments_3d(dat0, /no_unit_conversion)
            ;if moms[jj].tavg lt 0 then stop
        endfor
        store_data, the_species+'_moms', 0, moms
    endforeach

;---Save data to cdf file.
    if n_elements(cdf_file) eq 0 then begin
        errmsg = handle_error('No output CDF file ...')
        return
    endif
    if file_test(cdf_file) ne 0 then file_delete, cdf_file
    ginfo = {$
        title: 'RBSP HOPE moments, calculated based on l2 data',$
        text: 'Calculated by Sheng Tian at the University of Minnesota, email:tianx138@umn.edu'}
    scdfwrite, cdf_file, gattribute=ginfo

    ; ut time.
    ainfo = {$
        FIELDNAM:'UT sec',$
        UNITS:'s (UT)',$
        VAR_TYPE:'support_data'}
    foreach the_species, ['e','p'] do begin
        the_type = (the_species eq 'e')? 'ele': 'ion'
        utname = 'ut_'+the_type
        moms = get_var_data(the_species+'_moms')
        uts = moms.ut
        scdfwrite, cdf_file, utname, value=uts, cdftype='CDF_DOUBLE', attribute=ainfo
        get_data, pre0+'q_uvw2gse', tuts, quvw2gse
        tval = qslerp(quvw2gse, tuts, uts)
        store_data, pre0+'q_uvw2gse_'+the_type, uts, tval
    endforeach

    foreach the_species, species do begin
        the_type = (the_species eq 'e')? 'ele': 'ion'
        utname = 'ut_'+the_type
        moms = get_var_data(the_species+'_moms')

        ; Prepare the rotation matrix from UVW to GSE.
        get_data, pre0+'q_uvw2gse_'+the_type, uts, quvw2gse
        muvw2gse = transpose(qtom(quvw2gse))
        nrec = n_elements(uts)

        ; Density.
        vname = the_species+'_density'
        tval = moms.density
        ainfo = {$
            fieldnam: 'Density',$
            units: 'cm!U-3!N',$
            var_type: 'data',$
            depend_0: utname}
        scdfwrite, cdf_file, vname, value=tval, attribute=ainfo

        ; Average temperature.
        vname = the_species+'_t_avg'
        tval = moms.tavg
        ainfo = {$
            fieldnam: 'Average temperature',$
            units: 'eV',$
            var_type: 'data',$
            depend_0: utname}
        scdfwrite, cdf_file, vname, value=tval, attribute=ainfo

        ; Bulk velocity.
        vname = the_species+'_vbulk'
        tval = transpose(moms.vbulk)
        ainfo = {$
            fieldnam: 'Bulk velocity',$
            units: 'km/s',$
            var_type: 'data',$
            coordinate: 'GSE',$
            depend_0: utname}
        if keyword_set(test_cristian) then tval[*,0:1] *= -1
        for jj=0, nrec-1 do tval[jj,*] = tval[jj,*] # muvw2gse[*,*,jj]
        scdfwrite, cdf_file, vname, value=transpose(tval), attribute=ainfo, dimensions=[3], dimvary=[1]

        ; Number flux.
        vname = the_species+'_number_flux'
        tval = transpose(moms.nflux)
        ainfo = {$
            fieldnam: 'Number flux',$
            units: '#/(cm!U2!N-s)',$
            var_type: 'data',$
            coordinate: 'GSE',$
            depend_0: utname}
        if keyword_set(test_cristian) then tval[*,0:1] *= -1
        for jj=0, nrec-1 do tval[jj,*] = tval[jj,*] # muvw2gse[*,*,jj]
        scdfwrite, cdf_file, vname, value=transpose(tval), attribute=ainfo, dimensions=[3], dimvary=[1]

        ; Energy flux.
        vname = the_species+'_energy_flux'
        tval = transpose(moms.eflux)
        ainfo = {$
            fieldnam: 'Energy flux',$
            units: 'mW/m!U2!N',$
            var_type: 'data',$
            coordinate: 'GSE',$
            depend_0: utname}
        if keyword_set(test_cristian) then tval[*,0:1] *= -1
        for jj=0, nrec-1 do tval[jj,*] = tval[jj,*] # muvw2gse[*,*,jj]
        scdfwrite, cdf_file, vname, value=transpose(tval), attribute=ainfo, dimensions=[3], dimvary=[1]

        ; Enthalpy.
        vname = the_species+'_enthalpy'
        tval = transpose(moms.enthpy)
        ainfo = {$
            fieldnam:'Enthalpy',$
            units:'mW/m!U2!N',$
            var_type:'data',$
            coordinate: 'GSE',$
            depend_0:utname}
        if keyword_set(test_cristian) then tval[*,0:1] *= -1
        for jj=0, nrec-1 do tval[jj,*] = tval[jj,*] # muvw2gse[*,*,jj]
        scdfwrite, cdf_file, vname, value=transpose(tval), attribute=ainfo, dimensions=[3], dimvary=[1]
    endforeach

end


pro rbsp_read_hope_moments, time, probe=probe, errmsg=errmsg, $
    species=the_species, energy_range=energy_range, $
    renew_file=renew_file, version=version, local_root=local_root, _extra=extra

    compile_opt idl2
    on_error, 0
    errmsg = ''

    xyz = ['x','y','z']
    rgb = sgcolor(['red','green','blue'])
    string_gamma = '!9'+string(71b)+'!X'


;---Check inputs.
    if n_elements(time) eq 0 then begin
        errmsg = handle_error('No input time ...')
        return
    endif
    if size(time,/type) eq 7 then time = time_double(time)
    cadence = 86400d    ; sec.
    file_times = break_down_times(time, cadence)

    if n_elements(probe) ne 1 then begin
        errmsg = handle_error('Invalid probe ...')
        return
    endif
    if n_elements(the_species) eq 0 then begin
        errmsg = handle_error('No species ...')
        return
    endif
    the_type = (the_species eq 'e')? 'ele': 'ion'
    species = ['e','p','o','he']
    index = where(species eq the_species, count)
    if count eq 0 then begin
        errmsg = handle_error('Invalid species: '+species+' ...')
        return
    endif
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'sdata','rbsp'])

    if n_elements(energy_range) eq 2 then begin
        renew_file = 1
        case the_species of
            'e': electron_energy_range = energy_range
            'p': ion_energy_range = energy_range
            'o': oxygen_energy_range = energy_range
            'he': helium_energy_range = energy_range
            else: ;
        endcase
    endif


;---Prepare for the file name.
    base_pattern = 'rbsp'+probe+'_hope_moments_%Y_%m%d_'+version+'.cdf'
    local_paths = [local_root,'rbsp'+probe,'hope_moment','%Y']
    files = list()
    foreach file_time, file_times do begin
        file = apply_time_to_pattern(join_path([local_paths,base_pattern]), file_time)
        if keyword_set(renew_file) then if file_test(file) eq 1 then file_delete, file
        if file_test(file) eq 0 then rbsp_gen_hope_moments, file_time+[0,cadence], probe=probe, filename=file, errmsg=errmsg, $
            electron_energy_range=electron_energy_range, ion_energy_range=ion_energy_range, $
            oxygen_energy_range=oxygen_energy_range, helium_energy_range=helium_energy_range, _extra=extra
        if file_test(file) eq 1 then files.add, file
    endforeach
    files = files.toarray()


    in_vars = the_species+'_'+['density','t_avg','vbulk','number_flux','energy_flux','enthalpy']
    out_vars = 'rbsp'+probe+'_'+in_vars
    time_var_name = 'ut_'+the_type
    time_var_type = 'unix'
    request = dictionary($
        'var_list', list($
            dictionary($
                'in_vars', in_vars, $
                'out_vars', out_vars, $
                'time_var_name', time_var_name, $
                'time_var_type', time_var_type )))
    read_files, time, files=files, request=request


;---Configurate the variables.
    var = 'rbsp'+probe+'_'+the_species+'_density'
    add_setting, var, /smart, {$
        display_type: 'scalar', $
        unit: 'cm!U-3!N', $
        short_name: 'N!D'+the_species+'!N', $
        ylog: 1}

    var = 'rbsp'+probe+'_'+the_species+'_t_avg'
    add_setting, var, /smart, {$
        display_type: 'scalar', $
        unit: 'eV', $
        short_name: 'T!D'+the_species+'!N', $
        ylog: 1}

    var = 'rbsp'+probe+'_'+the_species+'_vbulk'
    add_setting, var, /smart, {$
        display_type: 'vector', $
        unit: 'km/s', $
        short_name: 'V!S!U'+the_species+'!N!R', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}

    var = 'rbsp'+probe+'_'+the_species+'_number_flux'
    add_setting, var, /smart, {$
        display_type: 'vector', $
        unit: '#/cm!U2!N-s', $
        short_name: 'F!S!U'+the_species+'!N!R', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}

    var = 'rbsp'+probe+'_'+the_species+'_energy_flux'
    add_setting, var, /smart, {$
        display_type: 'vector', $
        unit: 'mW/m!U2!N', $
        short_name: string_gamma+'!S!U'+the_species+'!N!R', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}

    var = 'rbsp'+probe+'_'+the_species+'_enthalpy'
    add_setting, var, /smart, {$
        display_type: 'vector', $
        unit: 'mW/m!U2!N', $
        short_name: 'H!S!U'+the_species+'!N!R', $
        coord: 'GSM', $
        coord_labels: xyz, $
        colors: rgb}

end

;time_range0 = time_double(['2012-10-25','2016-12-31'])
;time_range0 = time_double(['2016-05-17','2016-12-31'])
;;time_range0 = time_double(['2014-08-27','2014-08-28'])
;secofday = 86400d
;
;nday = (time_range0[1]-time_range0[0])/secofday
;ut1s = smkarthm(time_range0[0], time_range0[1], nday+1, 'n')
;ut2s = ut1s+secofday
;
;
;for i=0, nday-1 do begin
;    time = [ut1s[i],ut2s[i]]
;    rbsp_gen_hope_moments, time
;endfor


time = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
probe = 'a'
pre0 = 'rbsp'+probe+'_'
species = 'o'
rbsp_read_hope_moments, time, probe=probe, species=species;, /renew_file
;rbsp_gen_hope_moments, time, probe='a', file=shomedir()+'/test.cdf'
rbsp_read_orbit, time, probe=probe
rbsp_read_bfield, time, probe=probe
bvar = pre0+'b_gsm'
rvar = pre0+'r_gsm'
define_fac, bvar, rvar
vars = pre0+species+'_'+['vbulk','number_flux','energy_flux','enthalpy']
foreach var, vars do to_fac, var

end
