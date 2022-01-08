;+
; Type:
;   function.
; 
; Name:
;   sgse2mgse.
;
; Purpose:
;   Convert vector in GSE coord to MGSE.
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3], required.
;       Vector(s) in GSE.
;     
;   et0, in, type = double/dblarr[n], required.
;       The UT epoch(s).
;   
; Keywords:
;   probe = probe, in, type = 'a' or 'b', optional.
;       Select rbspa or rbspb. Default is 'a'.
;       Deal with ONLY one probe, ['b','a'] is equivalent to 'a'.
;       
;  wsc = wsc, in, type = dblarr[n,3], optional.
;      spice runs slow, spacecraft spin direction changes slow too.
;      Thus we can get wsc at several epochs and interpolate.
;   
; Return:
;   return, out, type = dblarr[3]/dblarr[3,n].
;       Vector(s) in MGSE. n is the num of records.
;   
; Example:
;   vec1 = sgse2mgse(vec0, et).
;   
; Dependence:
;   none.
;   
; Notes:
;   Keep the unit of input vector in Cartesian coordinate.
;  
; Author:
;   Sheng Tian.
; 
; History:
;   2013-09-13, Sheng Tian, create.
;-

function sgse2mgse, vec0, et0, probe = probe0, wsc = wsc

    compile_opt idl2
    on_error, 2
    
    if n_elements(probe0) eq 0 then probe = 'a' else probe = probe0
    probe = strupcase(probe[0])
    
    vec1 = double(vec0)
    n1 = n_elements(vec1)/3 & n2 = n1+n1 & n3 = n2+n1
    vx0 = vec1[0:n1-1]
    vy0 = vec1[n1:n2-1]
    vz0 = vec1[n2:n3-1]
    
    ; get x_mgse, i.e., w_sc in gse.
    if n_elements(wsc) eq 0 then begin
        cdf_epoch, et0, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        tmp = string(yr, mo, dy, hr, mi, sc, msc, $
            format = '(I04,"-",I02,"-",I02,"T",I02,":",I02,":",I02,".",I06)')
        cspice_str2et, tmp, et      ; bad name, this et is spice epoch.
        cspice_pxform, 'RBSP'+probe+'_SCIENCE', 'GSE', et, pxform
        wx = pxform[2,0] & wy = pxform[2,1] & wz = pxform[2,2]
    endif else begin
        wx = wsc[*,0] & wy = wsc[*,1] & wz = wsc[*,2]
    endelse
    
    ; do rotation.
    p = atan(double(wy),wx)     ; this way p (phi) in [0,2pi].
    cosp = cos(p)
    sint = wx/cosp
    sinp = wy/sint
    cost = double(wz)
    
    vx1 =  cosp*vx0 + sinp*vy0
    vy1 = -sinp*vx0 + cosp*vy0
    vz1 =  vz0
    vx2 =  sint*vx1 + cost*vz1
    vz2 = -cost*vx1 + sint*vz1
    
    vec1[0:n1-1] = temporary(vx2)
    vec1[n1:n2-1] = temporary(vy1)
    vec1[n2:n3-1] = temporary(vz2)
    return, vec1
    
end