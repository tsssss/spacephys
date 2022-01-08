;+
; Type: procedure.
; Purpose: Rotate Polar E/B fields from SPC to FAC coord.
; Parameters: bmod_spc, in, dblarr[n,3], required. Model B SPC.
;   db_spc, in, dblarr[n,3], required. dB SPC.
;   de_spc, in, dblarr[n,3], required. dE SPC.
;   db_fac, out, dblarr[n,3], required. dB FAC.
;   de_fac, out, dblarr[n,3], required. dE FAC.
;   lon, out, dblarr[n], optional. Angle rotate to dipole plane.
;   lat, out, dblarr[n], optional. Angle rotate to B.
; Keywords: none.
; Notes: {SPC. xy:away sun, z:north, 56:spin axis}{FAC. v:perp B, along v sense,
;   p:bxv, b:parallel B}. For noon-midnight meridian, xy~v, z~b, 56~p.
; Dependence: none.
; Author: Sheng Tian.
; History: 2013-12-01, Sheng Tian, create.
;-
pro polar_sdt_spc2fac, bmod_spc, db_spc, de_spc, $
    db_fac, de_fac, b1 = b1, b2 = b2
    
    bmod = bmod_spc
    ; b1: eliminate bmod_56.
    b1 = atan(bmod[*,2],bmod[*,0])
    srotate, bmod, b1, 1
    ; b2: eliminate bmod_xy.
    b2 = atan(bmod[*,0],bmod[*,1])
    srotate, bmod, b2, 2

    if n_elements(db_spc) eq 0 then return
    
    tmp = db_spc
    srotate, tmp, b1, 1
    srotate, tmp, b2, 2
    db_fac = [[tmp[*,0]],[-tmp[*,2]],[tmp[*,1]]]
    
    tmp = de_spc
    srotate, tmp, b1, 1
    srotate, tmp, b2, 2
    de_fac = [[tmp[*,0]],[-tmp[*,2]],[tmp[*,1]]]
    
;    ; in N-hem, z is -b, need rotation 180 around 56.
;    a1 = (ilat[0] gt 0)*2*!dpi
;    if a1 ne 0 then srotate, tmp, a1, 2
;    ; rotate 180 around z to make xy in the v sense.
;    ; only works for cusp, direct use v is better.
;    a2 = (ilat[0] lt ilat[-1])*2*!dpi
;    if a2 ne 0 then srotate, tmp, a2, 1

end