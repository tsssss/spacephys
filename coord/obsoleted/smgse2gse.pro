;+
; Type:
;   function.
; 
; Name:
;   smgse2gse.
;
; Purpose:
;   Convert vector in mGSE coord to GSE.
;
; Parameters:
;   vec0, in, type = dblarr[3]/dblarr[n,3], required.
;       Vector(s) in mGSE.
;     
;   et0, in, type = double/dblarr[n], optional.
;       The UT epoch(s).
;       When wsc is set, et0 is optional.
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
;       Vector(s) in GSE. n is the num of records.
;   
; Example:
;   vec1 = smgse2gse(vec0, et).
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

function smgse2gse, vec0, et0, probe = probe0, wsc = wsc

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
    
    vx1 = sint*vx0 - cost*vz0
    vz1 = cost*vx0 + sint*vz0
    vy1 = vy0
    vx2 = cosp*vx1 - sinp*vy1
    vy2 = sinp*vx1 + cosp*vy1
    
    vec1[0:n1-1] = temporary(vx2)
    vec1[n1:n2-1] = temporary(vy2)
    vec1[n2:n3-1] = temporary(vz1)
    return, vec1
    
end

rbsp_efw_init
!rbsp_efw.user_agent = ''

probe = 'a'

date = '2013-01-14'
duration = 1
timespan, date, duration

type = 'calibrated'

rbsp_load_spice_kernels
rbsp_load_spice_state,probe='a',coord='gse',/no_spice_load  
rbsp_load_spice_state,probe='b',coord='gse',/no_spice_load

time2=time_double(date)
time2 = [time2,time2+86399.]
time3=time_string(time2, prec=6)
strput,time3,'T',10
cspice_str2et,time3,et2

cspice_pxform,'RBSP'+strupcase(probe)+'_SCIENCE','GSE',et2[0],pxform1
cspice_pxform,'RBSP'+strupcase(probe)+'_SCIENCE','GSE',et2[1],pxform2
        
wsc = [0d,0,1]
wsc_GSE1 = pxform1 ## wsc
wsc_GSE2 = pxform2 ## wsc
        
rbsp_gse2mgse,'rbspa_state_pos_gse',reform(wsc_GSE1),$
    newname='rbspa_state_pos_mgse'

tplot, 'rbspa_state_pos_'+['gse','mgse']

get_data, 'rbspa_state_pos_gse', data = dat0
get_data, 'rbspa_state_pos_mgse', data = dat1

i = 0
ut = dat0.x[i]
et = stoepoch(ut,'unix')
v0 = dat0.y[i,*]
v1 = sgse2mgse(v0, et)
v2 = dat1.y[i,*]
v3 = smgse2gse(v1, et)
print, reform(v0), sqrt(total(v0^2))
print, reform(v1), sqrt(total(v1^2))
print, reform(v2), sqrt(total(v2^2))
print, reform(v3), sqrt(total(v3^2))
end