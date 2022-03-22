;+
; Read RBSP HOPE moments.
;-

pro rbsp_read_hope_moments, input_time_range, probe=probe, errmsg=errmsg

    files = rbsp_load_hope_moments(input_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return

    cdf2tplot, files

end
