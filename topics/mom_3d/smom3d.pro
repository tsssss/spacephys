;+
; Passed checks on moments_3d, for density to efluxes.
; 
; This program calculates moments in a despun coord.
; Ideally, it needs B vector and sc potential.
;-
function smom3d, data, vsc=vsc, dmom=dmom, erange=erange

;---Constants.
    mp = 1.67d-27   ; kg.
    qe = 1.6d-19    ; Coulume.
    idx6 = [0,4,8,1,2,5]    ; map from matrix[3x3] to vec[6].
    idx3x3 = [[0,3,4],[3,1,5],[4,5,2]]  ; from vec[6] to matrix[3x3].
    deg = 180d/!dpi
    rad = !dpi/180d
    
    mom3d = smom3d_init()


;---Collect basic info from data.
    dat3d = data
    if dat3d.valid eq 0 then return, mom3d
    
    ut0 = dat3d.time    ; ut in sec.
    m_e = dat3d.mass    ; mass(kg) over (qe*1e6), so that sqrt(E(eV)/m_e) in km/s.
    q0 = dat3d.charge   ; q/qe.
    bvec = dat3d.magf   ; 3D B field, in nT.

    ; sc potential in V.
    if n_elements(vsc) eq 0 then vsc = 0d
    if stagexist('sc_pot',dat3d) then vsc = dat3d.sc_pot
    if ~finite(vsc) then vsc = 6d

    ; differential energy flux, in eV/cm^2-eV-sr-s.
    fdat = dat3d.data
    idx = where(finite(fdat,/nan), cnt)
    if cnt ne 0 then fdat[idx] = 0

    ; energy bins and related.
    e0s = dat3d.energy  ; E, energy bins in eV.
    idx = where(e0s le 0, cnt)
    if cnt ne 0 then e0s = e0s>0.1  ; copied from L168 in moments_3d.
    ne0 = dat3d.nenergy

    if stagexist('denergy',dat3d) then begin
        des = dat3d.denergy
        de_e = des/e0s
    endif else begin
        de_e = abs(shift(e0s,1)-shift(e0s,-1))*0.5/e0s
        de_e[0,*] = de_e[ne0-2,*]
        de_e[ne0-1,*] = de_e[ne0-2,*]
        des = de_e*e0s
    endelse
    e1s = e0s+q0*vsc    ; E after sc potential correction.
    weight = 0d> (e0s+q0*vsc)/des+0.5 <1d

    if n_elements(erange) eq 2 then begin   ; energy range in eV.
        idx = where(e0s lt erange[0] or e0s gt erange[1], cnt)
        if cnt ne 0 then fdat[idx,*] = 0d
    endif 



    ; integration volumes.
    if n_elements(domgs) eq 0 then $
        domgs = smom3d_domega(dat3d.theta,dat3d.phi,dat3d.dtheta,dat3d.dphi)
    do0 = domgs[0,*,*,*]
    dox = domgs[1,*,*,*]
    doy = domgs[2,*,*,*]
    doz = domgs[3,*,*,*]



;---Calculate the moments.
    ; The essential quantities need be integrated are:
    ; dens (density), nflux (number flux), vften (velocity flux tensor), eflux (energy flux).
    ; jsc_e (total current over q)
    ; 
    ; Other quantites are derived quantities:
    ; mftens (mass flux tensor, from vftens)
    ; v_bk (bulk velocity, from nflux and dens)
    ; ptens (pressure tensor, from vftens, dens, and nflux)
    ; tavg (averaged temperature, from ptens and dens)
    ; vthermal (thermal velocity, from tavg)
    ; t3 (temperature eigen values, from ptens and dens)
    ; t3_mag (temperature eigen values in FAC, from ptens, dens, v_bk, and bvec)
    

    ; density in 1/cm^3.
    dens = total(fdat*de_e*weight*e1s^0.5/e0s*do0,nan=1)*(sqrt(m_e*0.5)*1e-5)
    if dens lt 0 then dens = !values.f_nan
    if ~finite(dens) then return, mom3d

    ; total current over q into the s/c, in 1/cm^2-s.
    jsc_e = total(fdat*de_e*weight*do0)

    ; number flux, in 1/cm^2-s.
    nflux = [$
        total(fdat*de_e*weight*e1s/e0s*do0*dox), $
        total(fdat*de_e*weight*e1s/e0s*do0*doy), $
        total(fdat*de_e*weight*e1s/e0s*do0*doz)]

    ; velocity flux, in 1/cm-s^2.
    vften = [$
        total(fdat*de_e*weight*e1s^1.5/e0s*do0*dox^2), $
        total(fdat*de_e*weight*e1s^1.5/e0s*do0*doy^2), $
        total(fdat*de_e*weight*e1s^1.5/e0s*do0*doz^2), $
        total(fdat*de_e*weight*e1s^1.5/e0s*do0*dox*doy), $
        total(fdat*de_e*weight*e1s^1.5/e0s*do0*dox*doz), $
        total(fdat*de_e*weight*e1s^1.5/e0s*do0*doy*doz)]/(sqrt(m_e*0.5)*1e-5)

    ; mass flux tensor, in eV/cm^3.
    mftens = vften*m_e*1e-10

    ; energy flux, in mW/m^2. Originally in eV/cm^2-s.
    eflux = [$
        total(fdat*de_e*weight*e1s^2/e0s*do0*dox), $
        total(fdat*de_e*weight*e1s^2/e0s*do0*doy), $
        total(fdat*de_e*weight*e1s^2/e0s*do0*doz)]*1.6e-12

    ; bulk velocity, in km/s.
    v_bk = nflux/dens*1e-5

    ; pressure tensor, in eV/cm^3.
    pt3x3 = mftens[idx3x3]-(v_bk # nflux)*m_e*1e-5
    ptens = pt3x3[idx6]
    pavg = mean(pt3x3[[0,4,8]])

    ; average temperature, in eV.
    tavg = pavg/dens
    
    ; thermal velocity, in km/s.
    vthermal = sqrt(2d*tavg/m_e)

    ; ke flux, in mW/m^2.
    keflux = 0.5*dens*m_e*total(v_bk^2)*v_bk*1.6e-7 ; 1.6e-7 = 1e5*1.6e-12.

    ; enthalpy, in mW/m^2.
    enthalpy = 2.5*pavg*v_bk*1.6e-7

    ; heat flux, in mW/m^2.
    hflux = eflux-keflux-enthalpy


;---Collect results.
    mom3d.valid = 1
    mom3d.ut = ut0
    mom3d.vsc = vsc
    mom3d.mass = m_e
    mom3d.charge = q0
    mom3d.bvec = bvec
    mom3d.density = dens
    mom3d.nflux = nflux
    mom3d.jsc = jsc_e
    mom3d.mftens = mftens
    mom3d.ptens = ptens
    mom3d.eflux = eflux
    mom3d.keflux = keflux
    mom3d.enthpy = enthalpy
    mom3d.hflux = hflux
    mom3d.tavg = tavg
    mom3d.vthermal = vthermal
    mom3d.vbulk = v_bk



;---Calculate uncertainty.
    if ~arg_present(dmom) then return, mom3d
    dmom = mom3d

    ; need to convert 1 count to eflux.
    ddat = conv_units(fdat, 'eflux')    ; TODO, this is just a placeholder.
end
