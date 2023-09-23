;+
; Calculate moments for given species and time range.
;
; input_time_range. The time range.
; probe=. 'a' or 'b'.
; errmsg=. The output error massage. By default is ''.
; species=. By default is ['e','p','o','he'].
; coord=. By default is 'gse', can be 'gsm', 'sm', etc.
; ion_energy_range=. By default is >30 eV.
; electron_energy_range=. By default is >200 eV.
; oxygen_energy_range=. By default is the same as ion_energy_range.
; helium_energy_range=. By default is the same as ion_energy_range.
; release=. By default is 'rel04'.
; balance_pixel=. Set to fine tune (balance) the flux in opposite direction.
; ncycle=. By default is 1. The # of spin used to calculate moments.
; vsc_var=. A string storing the SC potential.
;-
pro rbsp_calc_hope_moments, input_time_range, probe=probe, errmsg=errmsg, $
    species=input_species, coord=coord, $
    ion_energy_range=ion_energy_range, electron_energy_range=electron_energy_range, $
    oxygen_energy_range=oxygen_energy_range, helium_energy_range=helium_energy_range, $
    release=release, $
    balance_pixel=balance_pixel, ncycle=ncycle, $
    vsc_var=vsc_var

    if n_elements(ion_energy_range) ne 2 then ion_energy_range = [30,1e6]
    if n_elements(electron_energy_range) ne 2 then electron_energy_range = [200,1e6]
    if n_elements(oxygen_energy_range) ne 2 then oxygen_energy_range = ion_energy_range
    if n_elements(helium_energy_range) ne 2 then helium_energy_range = ion_energy_range
    if n_elements(ncycle) eq 0 then ncycle = 1
    if n_elements(coord) eq 0 then coord = 'gse'
    coord_str = strupcase(coord)
    time_range = time_double(input_time_range)

;---Read data.
    files = rbsp_load_hope(time_range, probe=probe, id='l2%sector', release=release, errmsg=errmsg)
    if errmsg ne '' then return

    prefix = 'rbsp'+probe+'_'
    var_list = list()

    var_list.add, dictionary($
        'in_vars', [ 'Epoch_Ele_DELTA','HOPE_ENERGY_Ele','FEDU'], $
        'time_var_name', 'Epoch_Ele', $
        'time_var_type', 'Epoch' )
    var_list.add, dictionary($
        'in_vars', ['Epoch_Ion_DELTA','HOPE_ENERGY_Ion','FPDU','FODU','FHEDU'], $
        'time_var_name', 'Epoch_Ion', $
        'time_var_type', 'Epoch' )
    var_list.add, dictionary($
        'in_vars', ['Sector_Collapse_Cntr','Energy_Collapsed','Epoch','Epoch_Ele','Epoch_Ion'], $
        'generic_time', 1 )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return


;---Load quaternion data.
    rbsp_read_quaternion, time_range, probe=probe, errmsg=errmsg
    if errmsg ne '' then return
    q_var = prefix+'q_uvw2gse'



;---Constants.
    errmsg = ''
    re = constant('re')
    re1 = 1d/re
    p_mass = 1.67d-27   ; kg.
    e_mass = 0.91d-30   ; kg.
    xyz = ['x','y','z']
    uvw = ['u','v','w']
    fac = ['b','w','n']
    rgb = constant('rgb')
    species_info = dictionary()
    species_info['e']  = dictionary('q',-1, 'mass',e_mass, 'energy_range'  ,electron_energy_range )
    species_info['p']  = dictionary('q', 1, 'mass',p_mass, 'energy_range'  ,ion_energy_range )
    species_info['o']  = dictionary('q', 1, 'mass',p_mass*16,'energy_range',oxygen_energy_range )
    species_info['he'] = dictionary('q', 1, 'mass',p_mass*4,'energy_range' ,helium_energy_range )
    all_species = species_info.keys()
    dt_phase_ratio = 180d/360     ; in terms of 360 deg.

    strmu = '!9'+string(109b)+'!X'
    prefix = 'rbsp'+probe+'_'

;---loop through each species.
    if n_elements(input_species) eq 0 then input_species = all_species
    foreach the_species, input_species do begin
        charge_type = (the_species eq 'e')? 'Ele': 'Ion'

    ;---Get the flux data and the energy bins.
        ; fdat, ut0s, ut1s.
        flux_var = strupcase('f'+the_species+'du')
        get_data, flux_var, times, fdat
        fdat = fdat*1e-3    ; from 1/s-cm^2-sr-keV to 1/s-cm^2-sr-eV.
        index = where(fdat le 0, count)
        if count ne 0 then fdat[index] = 0
        nrec = n_elements(times)
        dims = size(reform(fdat[0,*,*,*]),/dimensions)
        nen0 = dims[0]
        nsec = dims[1]
        npol = dims[2]

        en_var = 'HOPE_ENERGY_'+charge_type
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
                store_data, prefix+the_species+'_coef', coef1, coef2
            endif else get_data, prefix+'p_coef', coef1, coef2    ; use h.

            ; Apply pixel balancing.
            for jj=0, nen0-1 do for k=0, npol-1 do fdat[*,jj,*,k] *= coef2[jj,k]
        endif


        ; Get the start/end times of an integration cycle.
        dt_var = 'Epoch_'+charge_type+'_DELTA'
        dtimes = get_var_data(dt_var)*2d-3
        ; The start and end of the integration times.
        ; By default, for 1 spin, the time is the center of the spin, i.e., at phi=180 deg.
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
        if n_elements(vsc_var) eq 0 then begin
            vscs = fltarr(nrec)+0d
        endif else begin
            vscs = get_var_data(vsc_var, at=ut0s+dtimes(dt_phase_ratio))
        endelse
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
            dat0.sc_pot = vscs[jj]
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
        store_data, prefix+the_species+'_moms', 0, moms

        ; Time and Quaternion.
        times = moms.ut
;        get_data, q_var, uts, q_uvw2gse
;        q_uvw2gse = qslerp(q_uvw2gse, uts, times)

        ; Density.
        density_var = prefix+the_species+'_density'
        store_data, density_var, times, moms.density
        add_setting, density_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'N!D'+the_species+'!N', $
            'unit', 'cm!U-3!N' )

        ; Average temperature.
        temp_var = prefix+the_species+'_t_avg'
        store_data, temp_var, times, moms.tavg
        add_setting, temp_var, smart=1, dictionary($
            'display_type', 'scalar', $
            'short_name', 'T!D'+the_species+'!N', $
            'unit', 'eV' )

        ; Bulk velocity.
        vel_var = prefix+the_species+'_vbulk_'+coord
        v_uvw = transpose(moms.vbulk)
        cotran_str = 'uvw2'+coord
        v_coord = cotran(v_uvw, times, cotran_str, probe=probe)
        store_data, vel_var, times, v_coord
        add_setting, vel_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'V!S!U'+the_species+'!N!R', $
            'unit', 'km/s', $
            'coord', coord_str )

        ; Number flux.
        nflux_var = prefix+the_species+'_nflux_'+coord
        nflux_uvw = transpose(moms.nflux)
        nflux_coord = cotran(nflux_uvw, times, cotran_str, probe=probe)
        store_data, nflux_var, times, nflux_coord
        add_setting, nflux_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'F!S!U'+the_species+'!N!R', $
            'unit', '#/cm!U2!N-s', $
            'coord', coord_str )

        ; Energy flux.
        eflux_var = prefix+the_species+'_eflux_'+coord
        eflux_uvw = transpose(moms.eflux)
        eflux_coord = cotran(eflux_uvw, times, cotran_str, probe=probe)
        store_data, eflux_var, times, eflux_coord
        add_setting, eflux_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', tex2str('Gamma')+'!S!U'+the_species+'!N!R', $
            'unit', 'mW/m!U2!N', $
            'coord', coord_str )

        ; Enthalpy.
        enthalpy_var = prefix+the_species+'_enthalpy_'+coord
        enthalpy_uvw = transpose(moms.enthpy)
        enthalpy_coord = cotran(enthalpy_uvw, times, cotran_str, probe=probe)
        store_data, enthalpy_var, times, enthalpy_coord
        add_setting, enthalpy_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'H!S!U'+the_species+'!N!R', $
            'unit', 'mW/m!U2!N', $
            'coord', coord_str )
    endforeach


end

time_range = ['2013-06-07','2013-06-08']
time_range = ['2015-02-18','2015-02-19']
time_range = ['2015-02-17/22:35','2015-02-18/07:25']
time_range = time_double(time_range)
probe = 'a'
prefix = 'rbsp'+probe+'_'

rbsp_calc_hope_moments, time_range, probe=probe, species=['e','p','o'], coord='gsm'
rbsp_read_orbit, time_range, probe=probe, coord='gsm'
rbsp_read_bfield, time_range, probe=probe, coord='gsm'
define_fac, prefix+'b_gsm', prefix+'r_gsm', time_var=prefix+'r_gsm'
rbsp_read_spice_var, time_range, probe=probe


vars = prefix+['p','o']+'_vbulk'
foreach var, vars do begin
    var_gsm = var+'_gsm'
    get_data, var_gsm, times, vec
    width = 600d/sdatarate(times)
    for ii=0,2 do begin
        vec[*,ii] = smooth(vec[*,ii], width, edge_zero=1, nan=1)
    endfor
    store_data, var_gsm, times, vec
    var_fac = var+'_fac'
    to_fac, var_gsm, to=var_fac
end

get_data, prefix+'p_density', times, p_dens
get_data, prefix+'o_density', times, o_dens
get_data, prefix+'p_vbulk_fac', times, vp_fac
get_data, prefix+'o_vbulk_fac', times, vo_fac, limits=lim

p_mass = 1d
o_mass = 16d
o_ratio = o_dens/(p_dens+o_dens)
p_ratio = 1-o_ratio
avg_mass = p_ratio*p_mass+o_ratio*o_mass
v_fac = vp_fac
for ii=0,2 do begin
    v_fac[*,ii] = (vp_fac[*,ii]*p_ratio*p_mass + vo_fac[*,ii]*o_ratio*o_mass)/avg_mass
endfor
store_data, prefix+'vbulk_fac', times, v_fac, limits=lim
options, prefix+'vbulk_fac', 'labels', 'MHD V'+['b','w','o']

tplot_options, 'constant', 0
get_data, prefix+'r_gsm', times, r_gsm
store_data, prefix+'dis', times, snorm(r_gsm)
get_data, prefix+'mlt', times, mlt
index = where(mlt ge 12, count)
if count ne 0 then mlt[index] -= 24
store_data, prefix+'mlt', times, mlt


tplot, prefix+[['o','p']+'_vbulk_fac','vbulk_fac','vf_fac','mlt','mlat','dis'], $
    trange=time_range
end
