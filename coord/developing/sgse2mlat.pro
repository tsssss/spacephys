;+
; Type: function.
; Purpose: Get mlat from vector in GSE.
; Parameters:
;   vec0, in, dblarr[3]/dblarr[n,3], req. Vector(s) in GSE.
;   et, in, double/dblarr[n], req. The UT epoch(s).
; Keywords:
;   degree, in, boolean, opt. Set to return in degree.
; Return: double/dblarr[n]. mlat in radian or degree.
; Dependence: slib/coord.
; Notes: none.
; History:
;   2014-06-17, Sheng Tian, create.
;-
function sgse2mlat, vec0, et, degree = degree
    compile_opt idl2 & on_error, 2
    
    tmp = sgsm2sm(sgse2gsm(vec0, et), et)
    tmp = cv_coord(from_rect = transpose(tmp), /to_sphere, degrees = degree)
    return, reform(tmp[1,*])

end