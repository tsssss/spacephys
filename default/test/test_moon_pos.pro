function vec3d, phi, theta, r
    if n_elements(r) eq 0 then r = 1d
    cosel = cos(theta)
    vec0 = r*cos(phi)*cosel
    vec1 = r*sin(phi)*cosel
    vec2 = r*sin(theta)
    return, transpose([vec0,vec1,vec2])
end

function vec3d_calcpolarangle, vec
    vec0 = sunitvec(vec)
    phi = atan(vec0[1],vec0[0])
    theta = atan(vec0[2],sqrt(vec0[0]^2+vec0[1]^2))
    return, {phi:phi, theta:theta}
end


pro equ2hor, dec, tau, lat, h, az
    e_equ = vec3d(tau,dec)
    e_hor = e_equ
    srotate, e_hor,-(!dpi*0.5-lat), 1
    
    tmp = vec3d_calcpolarangle(e_hor)
    az = tmp.phi
    h  = tmp.theta
end


pro test_moon_pos, et0, glat, glon

    rtod = 180d/!dpi
    dtor = !dpi/180d
    
    glat = 45*dtor     ; rad.
    glon =-93*dtor     ; rad.

    cdf_epoch, et0, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    ; eqn (4). for after 1900-03.
    jd = 367*yr-7*(yr+(mo+9)/12)/4+275*mo/9+dy+1721014
    print, jd
    mjd = sfmepoch(et0,'mjd')
;    
;    lm = 0.606434d + 0.03660110129d*t
;    gm = 0.374897d + 0.03629164709d*t
;    fm = 0.259091d + 0.03674819520d*t
;    d  = 0.827362d + 0.03386319198d*t
;    om = 0.347343d + 0.00014709391d*t
    
    ; adopted from sec 3.2, astronomy on the personal computer.
    j2000 = 2451545d        ; 2000 JAN 1.5.
    mj2000 = 51544.5         ; MJD_J2000.
    pi2 = 2*!dpi            ; 2*pi.
    arcs = 3600d*180d/!dpi  ; arcs.
    rad = !dpi/180d
    eps = 23.43929111d*rad  ; in rad.

    t = (mjd-mj2000)/36525  ; time in julian centuries since J2000.
    l0 = 0.606433d + 1336.855225d*t ; mean lon.
    l  = 0.374897d + 1325.552410d*t ; moon's mean anomaly.
    ls = 0.993133d +   99.997361d*t ; sun's mean anomaly.
    d  = 0.827361d + 1236.853086d*t ; diff lon moon-sun.
    f  = 0.259086d + 1342.227825d*t ; dist from ascending node.

    l0 = l0 mod 1d
    l  = pi2*(l mod 1d)
    d  = pi2*(d mod 1d)
    f  = pi2*(f mod 1d)

    ; perturbations in lon and lat.
    dl = +22640d*sin(l) - 4586d*sin(l-2*d) + 2370d*sin(2*d) + 769d*sin(2*l) - $
        668d*sin(ls) - 412d*sin(2*f) - 212d*sin(2*l-2*d) - 206d*sin(l+ls-2*d) + $
        192d*sin(l+2*d) - 165d*sin(ls-2*d) - 125d*sin(d) - 110d*sin(l+ls) + $
        148d*sin(l-ls) - 55d*sin(2*f-2*d)
    s = f + (dl+412d*sin(2*f)+541d*sin(ls))/arcs
    h = f-2*d;
    n = -526d*sin(h) + 44d*sin(l+h) - 31d*sin(-l+h) - 23d*sin(ls+h) + $
        11d*sin(-ls+h) - 25d*sin(-2*l+f) + 21d*sin(-l+f);

    ; ecliptic longitude and latitude.
    l_moon = pi2*((l0+dl/1296d3) mod 1) ; in rad.
    b_moon = (18520d*sin(s)+n)/arcs     ; in rad.

    ; equatorial coordinates.
    cosb = cos(b_moon) & sinb = sin(b_moon)
    cosl = cos(l_moon) & sinl = sin(l_moon)
    vec = transpose([cosb*cosl,cosb*sinl,sinb])
    
    ; rotate around x by -eps.
    srotate, vec, eps, 0
    
    ; ra and dec in rad.
    ra = atan(vec[1],vec[0])
    dec = asin(vec[2])
    tmp = (ra*rtod/15+24) mod 24
    print, 'ra: (HrAng) ', floor(tmp), (tmp-floor(tmp))*60
    tmp = (dec*rtod)
    print, 'dec: (deg)  ', floor(tmp), (tmp-floor(tmp))*60
    
    ; gmst.
    sgmst, et0, gmst, /radian
   
    ; hour angle tau, Eqn (3.8), in rad.
    tau = gmst + glon - ra
    
    ; altitude h, Eqn (3.5).
    clat = !dpi/2-glat  ; co-lat.
    sinh = cos(glat)*cos(tau)*cos(dec)+sin(glat)*sin(dec)
    equ2hor, dec, tau, glat, h, az
    
    print, 'MJD_J2000   ', mj2000
    print, 'MJD_t       ', jd
    print, 'RA (rad)    ', ra
    print, 'Dec (rad)   ', dec
    print, 'GMST (rad)  ', gmst
    print, 'glon (deg)  ', glon*rtod
    print, 'glat (deg)  ', glat*rtod
    print, 'az (deg)    ', (az*rtod+360) mod 360
    print, 'az (deg)    ', ((az*rtod+360) mod 360)+180
    print, 'alt (deg)   ', h*rtod
    print, 'alt2(deg)   ', asin(sinh)*rtod
    
    
;    ; test altitude.
;    rad = !dpi/180
;    deg = 180/!dpi
;    lon = 11.6d*rad          ; in rad.
;    lat = 48.1d*rad             ; in rad.
;    ra = 15d*(20+41d/60+12.8d/3600)*rad ; 20h 41' 12.8" -> rad.
;    dec = (45+15d/60+25.8d/3600)*rad    ; 45d 15' 25.8" -> rad.
;    et0 = stoepoch('1993-08-01/21')
;    mjd = sfmepoch(et0,'mjd')
;    sgmst, et0, tau0, /radian   ; in rad.
;    tau = tau0 + lon - ra
;    equ2hor, dec, tau, lat, h, az
;    print, 'MJD_t       ', mjd
;    print, 'GMST (rad)  ', tau0
;    print, 'glat (rad)  ', lat
;    print, 'glon (rad)  ', lon
;    print, 'al (deg)    ', h*deg
;    print, 'az (deg)    ', (az*deg+360) mod 360
end

et0 = stoepoch('1998-03-24')
et0 = stoepoch('2000-01-03')
et0 = stoepoch('2015-07-01/05:21')
;ut0 = sfmdate('2015-07-03/14:23 CDT','%Y-%m-%d/%H:%M %Z')
ut0 = sfmdate('2015-07-03/02:22 CDT', '%Y-%m-%d/%H:%M:%S %Z')
;ut0 = sfmdate('2015-02-01/03:00:00 GMT', '%Y-%m-%d/%H:%M:%S %Z')
et0 = stoepoch(ut0,'unix')
print, sfmepoch(et0)
test_moon_pos, et0

; test equiv. to MJD, JD. Passed.
; test equiv. to R_[xyz], Vec3D. Passed.
;pi = !dpi
;rad = pi/180d
;deg = 180d/pi
;eps = 23.5*rad
;e = transpose([1d,2,3])
;e = sunitvec(e) & a = e
;srotate, a, eps, 0
;print, 'l   ', atan(e[1],e[0])*deg
;print, 'b   ', asin(e[2])*deg
;print, 'RA  ', atan(a[1],a[0])*deg/15
;print, 'Dec ', asin(a[2])*deg
; test equiv. to Moon's RA and Dec. Passed.
; test equiv. to GMST. Passed.
end