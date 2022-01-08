;+
; Type:
;   procedure.
;
; Name:
;   sdipoledir.
;
; Purpose:
;   Return dipole direction (unit vector) in GEO, work within 1900 to 2015.
;
; Parameters:
;   et, in, type = double/dblarr[n], required.
;       Time in epoch UTC.
;   
;   v1, v2, v3, out, type = double/dblarr[n], required.
;       If t0 is double or 1-element array, they are in double,
;       if t0 is dblarr[n], they are in dblarr[n].
;       When degree or radian is set, v1 is for latitude, v2 longitude.
;       Otherwise, [v1,v2,v3] is [x,y,z].
;
; Keywords:
;   degree = degree, in, type = boolean, optional.
;       Set to return dipole direction in latitude, longitude in degree.
;
;   radian = radian, in, type = boolean, optional.
;       Set to return dipole direction in latitude, longitude in radian.
;       
;   interp = interp, in, type = boolean, optional.
;       Set to use method (a), see notes. To use method (a), t0 must be
;       between [1900,2015].
;
; Return:
;   none.
;
; Example:
;   sdipoledir, t0, dir, /degree
;
; Notes:
;   See Hapgood 1992, there are two ways to calculate the dipole GEO
;   latitude and longitude: (a) interpolate from IGRF model coefficients
;   (g01, g11, h11); (b) estimate from, e.g., the two equation below
;   equation (9) in Hapgood 1992.
;   
;   Method (b) is used by default. Set interp to use method (a). First, 
;   according to test_dipole_dir.pro, (a)'s error is < 0.01 deg, or 2e-4 rad,
;   (b)'s error is < ~0.3 deg, or 5e-3 rad. The accuray of (b) is generally
;   enough, plus (b) is much simpler and quicker to run.
;   
;   Method (a) use igrf11coeffs.txt at 
;   http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, which lists semi-
;   normalized spherical harmonics coefficient. Only n=1 harmonics (g01, g11, 
;   and h11) are used, they are the beginning coefficients, thus no Schidt
;   normalization is needed.
;   
;   NOTE: To vertorize the code, I modified the last column of table 
;   (i.e., 2010-2015). g01s is originally [-31543D, ..., -29496.5D, 11.4D], 
;   11.4D is changed to -29439.5 = -29496.5+11.4*5. Update similarly when
;   new IGRF coefficients are released.
;
; Dependence:
;   none.
;
; Author:
;   Sheng Tian.
;
; Hisotry:
;   2013-07-15, Sheng Tian, document.
;-

pro sdipoledir, et, v1, v2, v3, interp = interp, $
    degree = degree, radian = radian
    
    rtod = 180D/!dpi    
    mjd = (1D/86400000D)*et-678941D
    
    if ~keyword_set(interp) then begin
        t0 = mjd - 46066D
        v1 = 1.3753194505715316D + 0.0000020466107099D*t0    ; in radian.
        v2 = 5.0457468675156072D - 0.0000006751951087D*t0    ; in radian.
        if keyword_set(degree) or keyword_set(radian) then begin
            if keyword_set(degree) then begin
                v1 *= rtod
                v2 *= rtod
            endif
        endif else begin
            t = 0.5*!dpi - v1
            p = v2
            v1 = sin(t)*cos(p)
            v2 = sin(t)*sin(p)
            v3 = cos(t)
        endelse
        return
    endif
    
    ; IGRF coefficients, from 'igrf11coeffs.txt' 
    g01s = [-31543D, -31464D, -31354D, -31212D, -31060D, $
            -30926D, -30805D, -30715D, -30654D, -30594D, $
            -30554D, -30500D, -30421D, -30334D, -30220D, $
            -30100D, -29992D, -29873D, -29775D, -29692D, $
            -29619.4D, -29554.63D, -29496.5D, -29439.5D]
    
    g11s = [-2298D, -2298D, -2297D, -2306D, -2317D, $
            -2318D, -2316D, -2306D, -2292D, -2285D, $
            -2250D, -2215D, -2169D, -2119D, -2068D, $
            -2013D, -1956D, -1905D, -1848D, -1784D, $
            -1728.2D, -1669.05D, -1585.9D, -1502.4D]
    
    h11s = [5922D, 5909D, 5898D, 5875D, 5845D, $
            5817D, 5808D, 5812D, 5821D, 5810D, $
            5815D, 5820D, 5791D, 5776D, 5737D, $
            5675D, 5604D, 5500D, 5406D, 5306D, $
            5186.1D, 5077.99D, 4945.1D, 4801.1D]
            
    ; modified Julian date for the 1st of year 1900, 1905, ..., 2010, 2015.
    yrs = [15020L, 16846L, 18672L, 20498L, 22324L, 24151L, $
           25977L, 27803L, 29629L, 31456L, 33282L, 35108L, $
           36934L, 38761L, 40587L, 42413L, 44239L, 46066L, $
           47892L, 49718L, 51544L, 53371L, 55197L, 57023L]         
    
    ; only deal with 1900-2015.
    imin = min(mjd, max = imax)
    imax = (where(yrs ge imax))[0]
    imin = (where(yrs lt imin))[-1]
    if imin eq -1 or imin ge imax then message, 'invalid year ...'
    
    nrec = n_elements(mjd)
    g01 = dblarr(nrec)
    g11 = dblarr(nrec)
    h11 = dblarr(nrec)
    ; linear interpolate between 5 years.
    for i = imin, imax-1 do begin
        idx = where(mjd ge yrs[i] and mjd le yrs[i+1])
        if idx[0] eq -1 then continue
        tmp = (1D/(yrs[i+1]-yrs[i]))*(mjd[idx]-yrs[i])
        g01[idx] = g01s[i] + (g01s[i+1]-g01s[i])*tmp  ; -z.
        g11[idx] = g11s[i] + (g11s[i+1]-g11s[i])*tmp  ; -x.
        h11[idx] = h11s[i] + (h11s[i+1]-h11s[i])*tmp  ; -y.
    endfor
    
    if nrec eq 1 then begin
        g01 = g01[0] & g11 = g11[0] & h11 = h11[0]
    endif
    
    if keyword_set(radian) or keyword_set(degree) then begin
        v2 = atan(h11, g11)+!dpi
        v1 = atan(sin(v2)*g01/h11)
        if keyword_set(degree) then begin
            v1 *= rtod & v2 *= rtod
        endif
    endif else begin
        tmp = -1d/sqrt(g01^2+g11^2+h11^2)
        v1 = g11*tmp & v2 = h11*tmp & v3 = g01*tmp
    endelse
  
end