function slprtime, t0
  
  t1 = t0
  while t1 lt 0 do t1 += 24
  hr = fix(t1)
  t1 = (t1-hr)*60d
  mi = fix(t1)
  t1 = (t1-mi)*60d
  sc = fix(t1)
  t1 = (t1-sc)*10000
  msc = t1
  return, string(hr, mi, sc, msc, $
    format = '(I02,":",I02,":",I02,".",I04)')

end

; ==========
; time info.
; ==========
; source of reference, gmst.
; http://www.zeitladen.de/time.html
yr = 2013 & mo = 07 & dy = 01 & hr = 18 & mi = 00 & sc = 00
; gmst = 12:30:47.7.

; astronomical_almanac.pdf, for 2004.
yr = 2004 & mo = 04 & dy = 20 & hr = 00 & mi = 00 & sc = 00
; gmst = 13:53:39.9923.

; http://www2.arnes.si/~gljsentvid10/sidereal.htm
yr = 1994 & mo = 06 & dy = 16 & hr = 18 & mi = 00 & sc = 00
; gmst = 11:39:05.0675.

yr = 2004 & mo = 09 & dy = 14 & hr = 15 & mi = 59 & sc = 43.9729

; ======================
; prepare various times.
; ======================
cmj = 2400000.5D
; J2000.0: j2000 and mj2000.
j2000 = julday(01,01,2000,12)   ; in Julian date.
mj2000 = j2000-cmj              ; in modified Julian date.

; epoch: et.
cdf_epoch, et, yr, mo, dy, hr, mi, sc, /compute_epoch

; Julian date: jd0.
jd0 = julday(mo,dy,yr,0)

; Julian date time: jd.
jd = julday(mo,dy,yr,hr,mi,sc)

; modified Julian data time: mjd.
mjd = julday(mo,dy,yr,hr,mi,sc)-2400000.5

; modified Julian date: mjd0.
mjd0 = julday(mo,dy,yr,0)-2400000.5

; UT: ut.
ut = hr+(mi+sc/60D)/60D   ; in hour.

; ===========
; begin test.
; ===========
print, 'Input time: '+sfmepoch(et)
print, ''

; =============================
; Greenwich Mean Sidereal Time.
; =============================
  print, '**gmst...'
; this function.
  sgmst, et, gmst
  print, 'sgetgmst: '+slprtime(gmst)
  
; https://sites.google.com/site/physics135/Home/announcements/
; universaltimetolocalsiderealtimeconversion
  a = 6.656306D
  b = 0.0657098242D
  c = 1.0027379093D
  d = 2445700.5D
  gmst = (a+b*(jd0-d)+c*ut) mod 24   ; in hour.
  print, 'phys 135: '+slprtime(gmst)
  
; http://aa.usno.navy.mil/faq/docs/GAST.php.
  a = 18.697374558D
  b = 24.06570982441908D
  d = 2451545.0
  gmst = (a+b*(jd-d)) mod 24    ; in hour.
  print, 'usnonavy: '+slprtime(gmst)
  
; http://www2.arnes.si/~gljsentvid10/sidereal.htm.
  a = 280.46061837D
  b = 360.98564736629D
  gmst = (a+b*(jd-j2000)) mod 360      ; in degree.
  print, 'arnes.si: '+slprtime(gmst/15)
  
; Hapgood 1992.
  a = 100.461D
  b = 36000.770D/36525.0D
  c = 15.04107D
  d = 51544.5D
  gmst = (a+b*(jd0-j2000)+c*ut) mod 360   ; in degree.
  print, 'mhapgood: '+slprtime(gmst/15)
  
; almanac 2004.
  a = 24110.54841D/3600D
  b = 8640184.812866D/36525D/3600D
  gmst = (a+b*(jd-j2000)+ut) mod 24      ; in hour.
  print, 'almanc04: '+slprtime(gmst)

; almanac 1990.
  a = 6.69737456D
  b = 2400.051336D/36525D
  c = 1.002737909D
  gmst = (a+b*(jd0-j2000)+c*ut) mod 24      ; in hour.
  print, 'almanc90: '+slprtime(gmst) 
  
; =====================
; Sun direction angles.
; =====================
print, '**ssundir...'

; this function original.
  d = mjd-mj2000
  e = 23.439D - 0.00000036D*d
  g = 357.529D + 0.98560028D*d
  q = 280.459D + 0.98564736D*d
  gr = g*!dpi/180
  l = q + 1.915*sin(gr) + 0.020*sin(2*gr)
  print, 'getsundir: ', e, g, q, l

; Hapgoog.
  t0 = (mjd0-mj2000)/36525D
  e = 23.439D - 0.013*t0
  g = 357.528D + 35999.050D*t0 + 0.04107*ut
  gr = g*!dpi/180
  q = 280.460D + 36000.772D*t0 + 0.04107*ut
  l = q + (1.915D - 0.0048D*t0)*sin(gr) + 0.020D*sin(2*gr)
  print, 'mahapgood: ', e, g, q, l

; this function.
  ssundir, et, e, l, g, q
  print, 'radsundir: ', e*180/!dpi, g*180/!dpi, q*180/!dpi, l*180/!dpi
  ssundir, et, e, l, g, q, /degree
  print, 'degsundir: ', e, g, q, l


end