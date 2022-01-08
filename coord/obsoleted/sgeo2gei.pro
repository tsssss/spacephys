;+
; Type:
;   function.
; 
; Name:
;   sgeo2gei.
;
; Purpose:
;   Convert vector in GEO coord to GEI.
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3], required.
;       Vector(s) in GEO.
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
;       vec1 = smag2geo(vec0, epoch).
;   
; Dependence:
;       sgmst.
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

function sgeo2gei, vec0, et

    compile_opt idl2
    on_error, 2
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get transpose(t1).
    sgmst, et, gmst, /radian
    sing = sin(gmst)
    cosg = cos(gmst)
    
    ; vectorized, so should be faster than matrix ##.
    vx1 =  cosg*vx0 - sing*vy0
    vy1 =  sing*vx0 + cosg*vy0
    vz1 =  vz0
    
    vec1[0:n1-1] = temporary(vx1)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1

end

;  t1 = [[ cosg, sing, 0], $
;        [-sing, cosg, 0], $
;        [    0,    0, 1]]
