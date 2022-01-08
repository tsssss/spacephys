;+
; Calculate the time ranges in between perigees. Return time ranges in [n,2].
;
; time_range. A time range in UT.
; probe=. A string 'a' or 'b'.
;-

function pflux_grant_calculate_orbit_time_range, time_range, probe=probe, use_apogee=use_apogee

    orbit_period = 9*constant('secofhour')
    perigee_dis = 2
    retval = !null

    prefix = 'rbsp'+probe+'_'
    r_var = prefix+'r_gsm'
    if check_if_update(r_var, time_range) then rbsp_read_orbit, time_range, probe=probe
    dis = snorm(get_var_data(r_var,times=times))
    index = where(dis le perigee_dis, count)
    if count eq 0 then return, retval

    orbit_time_step = 60.
    perigee_time_ranges = time_to_range(times[index], time_step=orbit_time_step)
    nperigee_time_range = n_elements(perigee_time_ranges)/2
    perigee_times = dblarr(nperigee_time_range)
    for ii=0, nperigee_time_range-1 do perigee_times[ii] = mean(perigee_time_ranges[ii,*])
    norbit = n_elements(perigee_times)-1
    orbit_time_ranges = dblarr(norbit,2)
    orbit_time_ranges[*,0] = perigee_times[0:norbit-1]
    orbit_time_ranges[*,1] = perigee_times[1:norbit]

    if keyword_set(use_apogee) then begin
        apogee_times = (perigee_times[0:norbit-1]+perigee_times[1:norbit])*0.5
        apogee_times = [apogee_times[0]-orbit_period, $
            apogee_times, $
            apogee_times[norbit-1]+orbit_period]
        index = lazy_where(apogee_times, '[]', time_range)
        apogee_times = apogee_times[index]
        norbit = n_elements(apogee_times)-1
        orbit_time_ranges = dblarr(norbit,2)
        orbit_time_ranges[*,0] = apogee_times[0:norbit-1]
        orbit_time_ranges[*,1] = apogee_times[1:norbit]
    endif

    return, orbit_time_ranges

end

    time_range = time_double(['2012-10-01','2012-10-02'])
    orbit_time_ranges = pflux_grant_calculate_orbit_time_range(time_range, probe='a')
    print, time_string(orbit_time_ranges)
    orbit_time_ranges = pflux_grant_calculate_orbit_time_range(time_range, probe='a', /use_apogee)
    print, time_string(orbit_time_ranges)

end
