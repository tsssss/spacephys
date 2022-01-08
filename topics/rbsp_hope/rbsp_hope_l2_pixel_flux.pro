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

; calc diff flux between opposite pixels.

pro rbsp_hope_l2_pixel_flux, hope, pxs, type = vname0, index = enidx

    vars = strlowcase(tag_names(hope))

    t0 = sfmepoch(hope.(where(vars eq 'epoch')),'unix')
    smode = hope.(where(vars eq 'sector_collapse_cntr'))

    ion_t0 = sfmepoch(hope.(where(vars eq 'epoch_ion')),'unix')
    ion_dt = hope.(where(vars eq 'epoch_ion_delta'))*1d-3
    v0 = hope.(where(vars eq vname0))
    e0 = hope.(where(vars eq 'hope_energy_ion'))

    if n_elements(enidx) eq 0 then enidx = indgen(72)
    v0 = v0[*,enidx,*,*]
    e0 = e0[*,enidx]

    nrec = n_elements(ion_t0)
    nen = n_elements(e0)/nrec

    px1 = pxs[0]
    px2 = pxs[1]

    for j = 0, nen-1 do begin

        v1 = reform(v0[*,j,*,px1])      ; [n,16].
        v2 = reform(v0[*,j,*,px2])      ; [n,16].

        t1 = []
        f1 = []
        f2 = []
        e1 = []

        for i = 0ull, nrec-1 do begin

            tdt = ion_dt[i]
            tdt = ion_dt[i]             ; half spin in sec.
            tt0 = ion_t0[i]             ; start of the spin?
            te0 = e0[i,j]               ; energy.


            tsmode = (smode[where(t0 eq tt0),px1])[0]
            nsec = 16/tsmode
            secidx = fix(smkarthm(0,16/nsec,nsec,'x0'))

            tt1 = smkarthm(tt0,tdt*2/nsec,nsec,'x0')
            tv1 = reform(v1[i,*])       ; [16].
            tv2 = reform(v2[i,*])       ; [16].

            tv2 = shift(tv2,8)          ; shitf the 2nd pixel by 180 deg.
            te1 = te0 ## (intarr(nsec)+1)
            
            t1 = [t1,tt1]
            f1 = [f1,(tv1-tv2)[secidx]] ; diff flux.
            f2 = [f2,(tv1+tv2)[secidx]] ; total flux, use mean?
            e1 = [e1,te1[secidx]]

        endfor

        tvar = vname0+'_en'+string(enidx[j],format='(I02)')+ $
            '_df'+string(px1,format='(I0)')+string(px2,format='(I0)')
        store_data, tvar, t1, f1/f2, te1
        
        options, tvar, 'labels', snum2str(te1[0])+'eV'
        options, tvar, 'ytitle', 'Pixel '+string(px1,format='(I0)')+'-'+string(px2,format='(I0)')+'!CCommon Mode!CRejection Ratio'

    endfor

end
