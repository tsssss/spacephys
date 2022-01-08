pro polar_uvi_corr, et, orbit, system, lat, lon, img, not_lbhl = not_lbhl

    ; **** line of sight correction.
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    time = [yr*1000+stodoy(yr,mo,dy), 1000d*(sc+60d*(mi+60d*hr))+msc]
    rsat = sqrt(total(double(orbit)^2))
    ic_gci_to_geo,time,rotm
    o_geo = transpose(transpose(rotm) # orbit)
    o_geo = double(o_geo)/sqrt(total(double(orbit)^2))
    sslat = asin(o_geo(2))/!dtor
    sslon = atan(o_geo(1),o_geo(0))/!dtor
    
    dtor = !dpi/180.d
    re = 6371.d
    norm_rsat = rsat/re
    
    case system of
        1: solid_ang = 0.03380d*dtor*0.03987d*dtor ; from swift
        2: solid_ang = 0.03511d*dtor*0.04178d*dtor ; 10/1/97
    endcase
    nadir = (rsat - re)*(rsat - re)*solid_ang
    
    the = replicate(0.,200,228)
    phi = replicate(0.,200,228)
    
    ;comp_area = 90. * dtor
    area = replicate(0.,200,228)
    
    ;comp_los = 1.5097585 ; =86.5 degrees, max of function
    los = replicate(1.,200,228)
    inv_los = replicate(0.,200,228)
    
    indx = where(lat gt -90., count)
    
    if (count gt 0) then begin
        lat1 = sslat*dtor
        lat2 = lat*dtor
        lon1 = sslon*dtor
        lon2 = lon*dtor
        
        the = double(acos(sin(lat1)*sin(lat2) + $
            cos(lat1)*cos(lat2)*cos(lon1-lon2)))
        phi = double( atan( sin(the) , norm_rsat - cos(the) ) )
        
        ang = the + phi
        
        area(indx) = double(nadir/cos(ang(indx)))
        
        ; simple formula
        ; indx_los = where(ang lt comp_los)
        ; los(indx) = double( 1./cos(ang(indx)) )
        
        ; glynn's empirical fit
        ; los(indx) = double(1.d/cos(ang(indx)))
        los(indx) = double(exp(0.061*(1.-1./cos(ang(indx))))/cos(ang(indx)))
        
        ; find inverse of los.
        maxang = 78.*dtor
        maxlos = 5.  ;(maximum angle from nadir is 78 deg)
        gndx = where( (lat gt -90.) and (ang lt maxang), ngndx )
        if (ngndx gt 0) then inv_los[gndx] = 1./los[gndx]
        
        ; if not lbhl, don't do los correction, but remove large angle pixels
        if (keyword_set(not_lbhl)) then $
          if (ngndx gt 0) then inv_los[gndx] = 1.
    endif

    img /= inv_los

    ; **** dayglow correction.

    ; need gha, dec.
    greg_to_jdaynum,yr,mo,dy,jd ; julian day number.
    ut = time[1]/24d/3600000d   ; ut hour.
    ephem,jd,ut,gha,dec
    print, gha
    
;    ssunlon, et, slon, /deg & gha = -slon

    ;calculate hour angle in degrees
    ha = lon + gha

    ;check range of ha 
    ndx=where(ha lt 0,count)
    if(count gt 0) then ha(ndx)=ha(ndx)+360.
    ndx=where(ha gt 360.,count)
    if(count gt 0) then ha(ndx)=ha(ndx)-360.

    ;calculate solar zenith angle in radian.
    t1  = sin(dec*!dtor)*sin(lat*!dtor)
    t2  = cos(dec*!dtor)*cos(lat*!dtor)*cos(ha*!dtor)
    sza = acos(t1 + t2)
    sza = sza * !radeg  ; convert to deg.

    phi = 0.83d
    amp = 15d
    szang = findgen(151)
    model = amp*cos(phi*szang*!dtor)^2
    dum = min(model,subs)
    model[subs:150] = 0.
    spl = spl_init(szang,model)
    dayglo = spl_interp(szang,model,spl,sza)
    dayglo = reverse(dayglo, 2)
    img = img/2.5
    img = (img- dayglo)
    img[where(img lt 0)] = 0

end
