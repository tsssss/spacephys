
pro stplot_enspec2density, vname, mass = mass, newname = newname, keV = keV, electron = electron
    
    if n_elements(mass) eq 0 then mass = 1  ; proton.
    if keyword_set(electron) then mass = 1d/1836.15
    if n_elements(newname) eq 0 then newname = vname+'_density'
    
    get_data, vname, t0, dat, ens       ; dat in #/cm^2-s-sr-eV.
    dat*= 4*!dpi                        ; dat in #/cm^2-s-eV.
    
    if keyword_set(keV) then ens*= 1e3  ; make sure energy in eV.
    
    nen = n_elements(ens)               ; # of energy bands.
    nrec = n_elements(t0)               ; # of records.
    den = ens[1:*]-ens[0:nen-2] & den = [0,den] ; energy steps.
    
    cc = sqrt(2*1.6/1.67)*1e6           ; convert sqrt(eV/mp) to cm/s.
    v1s = 1d/(sqrt(ens/mass)*cc)        ; inverse of velocity, in s/cm.
    for j = 0, nrec-1 do begin
;        dat[j,*]*= den                  ; dat in #/cm^2-s.
        dat[j,*]*= v1s                  ; dat in #/cm^3.
    endfor
    store_data, newname, t0, total(dat,2)
end

