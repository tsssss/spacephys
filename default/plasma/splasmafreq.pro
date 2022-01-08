;+
; Function:
;    splasmafreq.
;
; Purpose:
;   I calculate the plasma frequency of
;   a plasma, given density and species.
;    
; Parameters:
;   density, in, type = double, required.
;     Density of a plasma, default unit is #/cm^-3.
;
;   species, in, type = int, optional.
;     Molecular mass, default is electron, otherwise ion.
;       
; Keywords:
;   charge, in, type = int, optional.
;     Set the charge state Z.
;
; Return:
;   return, out, type = double.
;     Plasma frequency of the plasma, unit is rad/sec.
;
; Example:
;   omega_poxygen = splasmafreq(10, 16)
;
; Notes:
;   * The density is in #/cm^-3 
;     and the species is in molecular mass.
;
; Dependence:
;   none.
;
; Author:
;   Sheng Tian.
;
; History:
;   2012-06-08, Sheng Tian, create.
;
;-

function splasmafreq, density, species

  compile_opt idl2

  on_error, 2
  
  ; omega_p = \sqrt{\frac{n q^2}{\varepsilon_0 m}}
  coeff = 5.64e4
    
  if n_elements(charge) eq 0 then $
    charge = 1
  
  if n_elements(species) eq 0 then $
    species = 1 $
  else coeff = 2.10e2
    
  return, coeff * charge * sqrt(density/species)

end 
