;+
; Function:
;    sdebyelength.
;
; Purpose:
;   This function calculates the Debye length of a plasma, given
;   given density, magnetic field strength.
;   
;   Note that the density is in #/cm^-3 
;   and the temperature is in eV.
;    
; Parameters:
;   density, in, type = double, required.
;     Density of a plasma, default unit is #/cm^-3.
;   
;   temperature, in, type = double, required.
;     Temperature of a plasma, default unit is eV.
;       
; Keywords:
;   kelvin, in, type = boolean, optional.
;     Set the temperature in K.
;   
;   charge, in, type = int, optional.
;     Set the charge state Z.
;
; Return:
;   return, out, type = double.
;     Debye length of the plasma, default unit is cm.
;
; Example:
;   debyelength = sdebyelength(10, 10, /Kelvin)
;
; Notes:
;   none.
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

function sdebyelength, density, temperature, $
  kelvin = kelvin, charge = charge

  compile_opt idl2

  on_error, 2
  
  ; lambda_D = \sqrt{\frac{\varepsilon_0 k_B T}{n e^2}}
  coeff = 743d
  
  if n_elements(kelvin) ne 0 then $
    temperature = temperature / 11604.505d
    
  if n_elements(charge) eq 0 then $
    charge = 1
    
  return, coeff * sqrt(temperature/(density * charge^2))

end 
