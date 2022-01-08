;+
; Type:
;   function.
; 
; Name:
;   sgse2gei.
;
; Purpose:
;   Convert vector in GSE coord to GEI.
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3], required.
;       Vector(s) in GSE.
;     
;   et, in, type = double/dblarr[n], required.
;       The UT epoch(s).
;   
; Keywords:
;   none.
;   
; Return:
;   return, out, type = dblarr[3]/dblarr[3,n].
;       Vector(s) in GEI. n is the num of records.
;   
; Example:
;   vec1 = sgse2gei(vec0, epoch).
;   
; Dependence:
;   ssundir.
;   
; Notes:
;   Keep the unit of input vector in Cartesian coordinate.
;  
; Author:
;   Sheng Tian.
; 
; History:
;   2012-07-13, Sheng Tian, create.
;-

function sgse2gei, vec0, et

    compile_opt idl2
    on_error, 2
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get transpose(t2).
    ssundir, et, e, l
    sine = sin(e)
    cose = cos(e)
    sinl = sin(l)
    cosl = cos(l)
    
    ; vectorized, so should be faster than matrix ##.
    tmp =  sinl*vx0 + cosl*vy0
    vx1 =  cosl*vx0 - sinl*vy0
    vy1 =  cose*tmp - sine*vz0
    vz1 =  sine*tmp + cose*vz0
    
    vec1[0:n1-1] = temporary(vx1)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1

end

; <l,z> = [[ cosl, sinl, 0D], $ 
;          [-sinl, cosl, 0D], $
;          [   0D,   0D, 1D]]
; <e,x> = [[1D,    0D,   0D], $
;          [0D,  cose, sine], $
;          [0D, -sine, cose]]
; t2 = <l,z> * <e,x>
;    = [[ cosl, sinl *cose, sinl *sine], $
;       [-sinl, cosl *cose, cosl *sine], $
;       [   0D,      -sine,       cose]]
