;+
; Function: sescapevel.
; Purpose: Calculate the escape velocity at given altitude.
; Parameters:
;   dis, in, double, req. Distance in Re.
;   species, in, int, opt. Molecular mass, default is electron, otherwise ion.
; Keywords:
;   energy, in, boolean, opt. Set to return escape energy in eV.
; Return: double. Escape velocity in km/s or energy in eV.
; Notes: none
; Dependence: none.
; History:
;   2016-08-14, Sheng Tian, create.
;-

function sescapevel, dis, species, energy = energy

  compile_opt idl2
  on_error, 2
  
    g0 = 6.67e-11   ; G.
    re = 6378e3     ; radius of earth in m.
    me = 5.972e24   ; mass of earth in kg.

    r = dis*re

    vesc = sqrt(2*g0*me/r)  ; escape velocity in m/s.
    if ~keyword_set(energy) then return, vesc*1e-3

    if n_elements(species) eq 0 then species = 1d/1836
    mp = 1.67e-27   ; proton mass in kg.
    qe = 1.6e-19    ; electron charge in C.
    eesc = 0.5*species*mp*vesc^2/qe
    return, eesc

end
