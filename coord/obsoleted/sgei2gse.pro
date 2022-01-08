;+
; Type: function.
; 
; Purpose: Convert vector in GEI coord to GSE.
;
; Parameters:
;   vec0, in, dblarr[3]/dblarr[n,3], required. Vector(s) in GEI.
;   et, in, double/dblarr[n], required. The UT epoch(s).
;   
; Keywords: none.
;   
; Return: dblarr[3]/dblarr[3,n]. Vector(s) in GSE. n is the num of records.
;   
; Notes: Keep the unit of input vector in Cartesian coordinate.
;
; Dependence: ssundir.
;  
; History:
;   2012-07-13, Sheng Tian, create.
;-

function sgei2gse, vec0, et

    compile_opt idl2
    on_error, 2
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get t2.
    ssundir, et, e, l
    sine = sin(e)
    cose = cos(e)
    sinl = sin(l)
    cosl = cos(l)
    
    ; vectorized, so should be faster than matrix ##.
    tmp =  cose*vy0 + sine*vz0
    vx1 =  cosl*vx0 + sinl*tmp
    vy1 = -sinl*vx0 + cosl*tmp
    vz1 = -sine*vy0 + cose*vz0
    
    vec1[0:n1-1] = temporary(vx1)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1

end

vec1 = sgei2gse(vec0, epoch).
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
