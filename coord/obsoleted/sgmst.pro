;+
; Type: procedure.
;
; Purpose: Calc the GMST (Greenwich Mean Sidereal Time) at given time in hour.
;
; Parameters:
;   t0, in, double/dblarr[n], required. Time in epoch UTC.
;   gmst, out, double/dblarr[n], required. GMST, default in hr, [0,24).
;
; Keywords:
;   degree, in, boolean, optional. Set to return angle in degree.
;   radian, in, boolean, optional. Set to return angle in radian.
;  
; Dependence: none.
;
; Notes: See details (notation and variable names) in Hapgood 1992.
;   doi: 10.1016/0032-0633(92)90012-D.
;   Coefficients are from Almanac for computers 1990.
;
; History:
;   Sheng Tian, 2012-07-03, create.
;   Sheng Tian, 2013-07-16, modify.
;-

pro sgmst, t0, gmst, radian = radian, degree = degree

  mjd = (1d/86400000d)*t0-678941d
  mjd0 = floor(mjd)         ; modified julian date for date only.
  ut = (mjd-mjd0)*24        ; fraction of day, convert to hour.

  a = 6.69737456d
  b = 0.0657098243942505d   ; = 2400.051336D/36525d
  c = 1.002737909d
  d = 51544.5d              ; J2000.0 in modified julian day.
  
  gmst = (((a+b*(mjd0-d)+c*ut) mod 24)+24) mod 24   ; in hour.
  if keyword_set(degree) then gmst *= 15            ; in degree.
  if keyword_set(radian) then gmst *= (!dpi/12)     ; in radian.

end