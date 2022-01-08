;+
; Generate the common_times for all variables.
;-

function cusp_ml_polar_test_common_time

    time_range = time_double(['1996-03-20','2008-06-14:23:59'])
    time_step = 60.
    common_times = make_bins(time_range, time_step)
    return, common_times

end
