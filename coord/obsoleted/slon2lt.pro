;+
; From longitude to local time. Eq (93) in https://arxiv.org/pdf/1611.10321.pdf.
;
; lon, input longitude, by default in hr, otherwise set degree or raidan.
; 
; default is for geo for glon to glt, set mag for mlon to mlt.
; lon in hour by default.
; 
; Return local time in hour in [-12,12].
;-


function slon2lt, lon, et, mag = mag, $
    degree = degree, radian = radian
    
    deg = 180d/!dpi
    d2h = 1d/15
    
    ssunlon, et, slon, mag=mag      ; slon in radian.
    slon = slon*deg*d2h             ; convert to hour.
    
    if keyword_set(degree) then begin
        lon0 = lon*d2h
    endif else if keyword_set(radian) then begin
        lon0 = lon*deg*d2h
    endif else begin
        lon0 = lon
    endelse
    
    lct = (lon0-slon+12) mod 24     ; in [0,24].
    idx = where(lct gt 12, cnt)
    if cnt ne 0 then lct[idx] -= 24 ; in [-12,12].
    
    return, lct
    
end