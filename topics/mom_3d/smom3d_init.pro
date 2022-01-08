
function smom3d_init, data, size=size

    if n_elements(data) ne 0 then begin
    sz = size(data,/dimensions)
    if keyword_set(size) then sz = data ; the input is actually the size.
    tdat = dblarr(sz)
    
    dat3d = {$
        time: 0d, $         ; ut in sec.
        end_time: 0d, $     ; end time in sec.
        charge: 0d, $       ; charge, C.
        mass: 0d, $         ; mass per eV*1e6.
        species: '', $      ; species.
        magf: dblarr(3), $  ; B field vector, nT.
        sc_pot: 0d, $       ; spacecraft potential, V.
        units: '', $        ; the type of input data.
        data: tdat, $       ; the input data.
        energy: tdat, $     ; energy bins in eV.
        denergy: tdat, $    ; width of the energy bins in eV.
        phi: tdat, $        ; azimuth angle in deg.
        dphi: tdat, $       ; width of the azimuth angle in deg.
        theta: tdat, $      ; elevation angle in deg, should in [-90,90].
        dtheta: tdat, $     ; width of the elevation angle in deg.
        nenergy: 0d, $      ; # of energy bins.
        bins: 1d, $         ; flags for which bins are used? obsoleted?
        valid: 0d}
    if ~keyword_set(size) then dat3d.data = data
        
    
    return, dat3d
    endif

    mom3d = {$
        ut: 0d, $           ; ut in sec.
        vsc: 0d, $          ; sc potential, V.
        mass: 0d, $         ; mass per eV*1e6.
        charge: 0d, $       ; charge, C.
        bvec:dblarr(3), $   ; B field vector, nT.
        density: 0d, $      ; density, 1/cm^3.
        nflux:dblarr(3), $  ; number flux, 1/cm^2-s.
        jsc: 0d, $          ; total current into sc, 1/cm^2-s.
        mftens:dblarr(6), $ ; mass flow tensor, eV/cm^3.
        ptens:dblarr(6), $  ; pressure tensor, eV/cm^3.
        eflux:dblarr(3), $  ; total energy flux, mW/m^2.
        keflux:dblarr(3), $ ; kinetic energy flux, mW/m^2.
        hflux:dblarr(3),$   ; heat flux, mW/m^2.
        enthpy:dblarr(3),$  ; enthalpy, mW/m^2.
        tavg: 0d, $         ; average temperature, eV.
        vthermal: 0d, $     ; thermal velocity, km/s.
        vbulk: dblarr(3), $ ; bulk velocity, km/s.
        valid: 0d}          ; valid when density is finite.

    return, mom3d

end