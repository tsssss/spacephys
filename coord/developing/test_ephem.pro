;+
;PROCEDURE: mag_to_geo
;
;PURPOSE:   Converts lattitude and longitude between MAG 
;       and GEO coordinates.
;       Uses a simple transformation matrix from Kivelson
;       and Russell, "Intro to Space Physics" which is not
;       very accurate in the polar regions.
;
;PARAMETERS:    lat The array of lattitudes.(In radians unless
;           degrees keyword set.)
;       lon The array of longitudes.(in radians unless
;           degrees keyword set.)
;
;KEYWORDS:  degrees Set this if both input and output are to be
;           in degrees.
;       mag Set this to do the inverse transformation,
;           GEO to MAG coordinates.
;Created by:    J.Rauchleiba    1/7/97
;-
pro mag_to_geo, lat, lon, degrees=deg, mag=mag

; Convert to radians if necessary

if keyword_set(deg) then begin
    lat = lat*!dtor
    lon = lon*!dtor
endif

; The transformation matrix

maggeo = [  [ .32110, .94498,  .06252], $
        [-.92756, .32713, -.18060], $
        [-.19112,      0,  .98157]  ]

if keyword_set(mag) then maggeo = transpose(maggeo)

lat = !pi/2. - lat  ; theta is measured from pole in spherics.

; Create array of column vectors (lats and lons in cartesian coords.)

Vmag = [    [sin(lat)*cos(lon)], $
        [sin(lat)*sin(lon)], $
        [cos(lat)     ] ]

; Transform each column vector

Vgeo = maggeo ## Vmag   ; array of 3-element column vectors

lat = acos( Vgeo(*,2) )
lon = atan( Vgeo(*,1)/Vgeo(*,0) )

flip = where(Vgeo(*,0) LT 0)    ; 0 < lon < pi, Vmag(1) > 0
if flip(0) NE -1 then $
    lon(flip) = lon(flip) + !pi ; pi < lon < 2pi, Vmag(1) < 0

lat = !pi/2. - lat  ; theta is measured from equator in GEO.

if keyword_set(deg) then begin
    lat = lat*!radeg
    lon = lon*!radeg
endif

return
end

yr = 2004 & mo = 09 & dy = 14 & hr = 15 & mi = 59 & sc = 43.9729
yr = 2000 & mo = 12 & dy = 03 & hr = 13 & mi = 48 & sc = 58.275
yr = 2000 & mo = 12 & dy = 03 & hr = 15 & mi = 10 & sc = 43.277

; J2000.0: j2000 and mj2000.
j2000 = julday(01,01,2000,12)   ; in Julian date.
mj2000 = 51544.5D               ; in modified Julian date.

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


ephem, jd0, ut, gha, dec, eqtime    ; in degree.
slon = (12. - ut - (eqtime/15.))*15.
print, 'ephempro: ', dec, slon, eqtime, ut
apexfile = '~/Dropbox/idl/idl82/Default/aurora/image/support/mlatlon.1997a.xdr'
apexfile = file_search(apexfile)
geotoapex, dec, slon, apexfile, mlat, mlon
print, 'imagefuv: ', slon, mlon

ssundir, et, e, l, g, q
dec = asin(sin(e)*sin(l))
ra = atan(cos(e)*tan(l))
eqtime = (q-ra-!dpi)*180/!dpi mod 360
slon = -ut*15 - eqtime + 180
print, 'thisfunc: ', dec*180/!dpi, slon, eqtime, ut

slon = ssunlon(et, /degree, /geo)
print, 'ssunlon.: ', slon, ssunlon(et, /degree)

mlat = dec*180/!dpi & mlon = slon
mag_to_geo, mlat, mlon, /mag, /degree
print, 'magtogeo: ', slon, mlon

tmp = scotrant5(et)

end

