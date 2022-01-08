;+
; Fix LANL orbit.
;-
;

pro fix_lanl_orbit, rvar_gsm, bad_time=bad_time, to=rvar_fix

    get_data, rvar_gsm, times, rgsm, limits=lim
    index = lazy_where(times, bad_time, count=count)
    if count eq 0 then return

    nan = !values.f_nan
    rgsm[index,*] = nan
    rgeo = cotran(rgsm, times, 'gsm2geo')
    dis = mean(snorm(rgeo),/nan)
    glon = atan(rgeo[*,1],rgeo[*,0])
    glat = asin(rgeo[*,2]/dis)

    period = 86400*0.5  ; sec.


    index = where(finite(glon))
    omega = 2*!dpi/period
    wt0 = omega*(times-times[0])
    wt1 = omega*(times[index]-times[0])
    x1 = sin(wt1)
    x2 = cos(wt1)
    x = [transpose(x1),transpose(x2)]

    y = glon[index]
    result = regress(x, y, const=const)
    glon1 = result[0]*sin(wt0)+result[1]*cos(wt0)+const


    index = where(finite(glon))
    omega = 1*!dpi/period
    wt0 = omega*(times-times[0])
    wt1 = omega*(times[index]-times[0])
    x1 = sin(wt1)
    x2 = cos(wt1)
    x = [transpose(x1),transpose(x2)]

    y = glat[index]
    result = regress(x, y, const=const)
    glat1 = result[0]*sin(wt0)+result[1]*cos(wt0)+const

;    colors = sgcolor(['black','red'])
;    labels = ['orig','fixed']
;    store_data, 'test_glon', times, [[glon],[glon1]], $
;        limits={colors:colors, labels:labels, ytitle:'GLon (rad)'}
;    store_data, 'test_glat', times, [[glat],[glat1]], $
;        limits={colors:colors, labels:labels, ytitle:'GLat (rad)'}
;    tplot, 'test_'+['glon','glat']

    rgeo[*,0] = dis*cos(glat1)*cos(glon1)
    rgeo[*,1] = dis*cos(glat1)*sin(glon1)
    rgeo[*,2] = dis*sin(glat1)

    rgsm = cotran(rgeo, times, 'geo2gsm')
    if n_elements(rvar_fix) eq 0 then rvar_fix = rvar_gsm
    store_data, rvar_fix, times, rgsm, limits=lim

end

time = time_double(['2014-08-28','2014-08-29'])
bad_time = time_double(['2014-08-28/05:30','2014-08-28/16:30'])
probe = '1991-080'
lanl_read_orbit, time, probe=probe
var = probe+'_r_gsm'
lanl_fix_orbit, var, bad_time=bad_time, to=var+'_fix'
end
