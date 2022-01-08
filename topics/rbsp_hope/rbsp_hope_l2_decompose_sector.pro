;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro rbsp_hope_l2_decompose_sector, hope, px, type = vname0

    vars = strlowcase(tag_names(hope))

    t0 = sfmepoch(hope.(where(vars eq 'epoch')),'unix')
    smode = hope.(where(vars eq 'sector_collapse_cntr'))

    ion_t0 = sfmepoch(hope.(where(vars eq 'epoch_ion')),'unix')
    ion_dt = hope.(where(vars eq 'epoch_ion_delta'))*1d-3
    v0 = hope.(where(vars eq vname0))
    e0 = hope.(where(vars eq 'hope_energy_ion'))

    
    nrec = n_elements(ion_t0)
    nen = n_elements(e0)/nrec


    t1 = []
    v1 = []
    e1 = []
    
    ; loop through each ion spin.
    for i = 0ull, nrec-1 do begin

        tdt = ion_dt[i]             ; half spin in sec.
        tt0 = ion_t0[i]             ; start of the spin?
        tv0 = reform(v0[i,*,*,px])  ; [72,16].
        te0 = e0[i,*]               ; [1,72].

        tsmode = (smode[where(t0 eq tt0),px])[0]
        nsec = 16/tsmode
        secidx = fix(smkarthm(0,16/nsec,nsec,'x0'))

        tt1 = smkarthm(tt0,tdt*2/nsec,nsec,'x0')
        tv1 = transpose(tv0[*,secidx])
        te1 = te0 ## (intarr(nsec)+1)

        nan = te0 & nan[*] = !values.d_nan
        t1 = [t1,tt0-tdt*2/nsec,tt1,tt0+tdt*2]
        v1 = [v1,nan,tv1,nan]
        e1 = [e1,te0,te1,te0]

    endfor

    tvar = vname0+'_px'+string(px,format='(I0)')
    store_data, tvar, t1, v1, e1
    options, tvar, 'ytitle', 'Pixel#'+string(px,format='(I0)')+'!CEnergy!C(eV)'
    options, tvar, 'ztitle', 's!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N'
    options, tvar, 'spec', 1
    options, tvar, 'no_interp', 1
    options, tvar, 'ylog', 1
    options, tvar, 'zlog', 1
    ylim, tvar, 1, 1e5, 1
    zlim, tvar, 1e2, 1e12, 1
    
end
