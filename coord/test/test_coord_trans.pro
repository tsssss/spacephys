; test coordinate transformations.
; (a) scotran.
; (b) tdas.
; (c) sscweb.

pro test_coord_trans, vec0, et

    rtod = 180D/!dpi
    dtor = !dpi/180D

    ; smag2geo and smag2geo.
    vgeo = vec0
    geo = [asin(vgeo[2]),atan(vgeo[1],vgeo[0])]
    glt = slon2lt(et, geo[1])
    vmag = sgeo2mag(vgeo, et)
    mag = sgeo2mag(geo, et, /radian)
    mlt = slon2lt(et, mag[1], /mag)
    
    ; sgeo2gei and sgei2geo.
    vgei = sgeo2gei(vgeo, et)
    gei = geo
    gei[0] = asin(vgei[2])
    gei[1] = atan(vgei[1],vgei[0])
    
    ; sgei2gse and sgse2gei.
    vgse = sgei2gse(vgei, et)
    gse = geo
    gse[0] = asin(vgse[2])
    gse[1] = atan(vgse[1],vgse[0])
    
    ; sgse2gsm and sgsm2gse.
    vgsm = sgse2gsm(vgse, et)
    gsm = geo
    gsm[0] = asin(vgsm[2])
    gsm[1] = atan(vgsm[1],vgsm[0])
    
    ; sgsm2sm and ssm2gsm.
    vsm = sgsm2sm(vgsm, et)
    sm = geo
    sm[0] = asin(vsm[2])
    sm[1] = atan(vsm[1],vsm[0])

    print, 'GEI    ', [vgei]
    print, 'GEO    ', [vgeo, glt]
    print, 'MAG    ', [vgeo, mlt]
    print, 'GSE    ', [vgse, 12+gse[1]*12/!dpi]
    print, 'GSM    ', [vgsm]
    print, 'SM     ', [vsm]

    ; =====
    ; tdas.
    ; =====
    cotrans_lib
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    tu = {year:yr, doy:stodoy(yr,mo,dy), hour:hr, min:mi, sec: sc, fsec: msc*0.001}
    tu = time_struct(sfmepoch(et))
    vgeo = transpose(vec0)
    subgeo2gei, tu, vgeo, vgei
    subgei2gse, tu, vgei, vgse
    subgse2gsm, tu, vgse, vgsm
    subgsm2sm, tu, vgsm, vsm
    print, 'GEI    ', transpose(vgei)
    print, 'GEO    ', transpose(vgeo)
    ; print, 'MAG     ', transpose(vmag)
    print, 'GSE    ', transpose(vgse)
    print, 'GSM    ', transpose(vgsm)
    print, 'SM     ', transpose(vsm)

end


et = stoepoch('1998-01-01 01:06:30')
et = stoepoch('2000-12-03 15:10:43')
et = stoepoch('2000-01-01 15:10:43')

mag = [80,50]*!dpi/180       ; in radian.
vec = [cos(mag[0])*cos(mag[1]),cos(mag[0])*sin(mag[1]),sin(mag[0])]

print, ''
print, '**time: ', sfmepoch(et)

test_coord_trans, vec, et

end