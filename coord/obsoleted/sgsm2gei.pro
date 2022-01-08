;+
; Type: function.
; Name: sgsm2gei.
; Purpose: Convert vector in GSM coord to GEI.
; Parameters:
;   vec0, in, dblarr[3]/dblarr[n,3], req. Vector(s) in GSM.
;   et, in, double/dblarr[n], req. The UT epoch(s).
; Keywords: none.
; Return: dblarr[3]/dblarr[3,n]. Vector(s) in GEI. n is the num of records.
; Dependence: slib/coord.
; Notes: Keep the unit of input vector in Cartesian coordinate.
; History:
;   2019-01-03, Sheng Tian, create.
;-
function sgsm2gei, vec0, et
    compile_opt idl2 & on_error, 2
    
    return, sgse2gei(sgsm2gse(vec0, et), et)

end
