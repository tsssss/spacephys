

    event_list = list()
    event_list.add, dictionary($
        'time_range', time_double(['2015-02-18/01:30','2015-02-18/04:00']), $
        'probes', ['a','b'] )
    
    
    foreach event, event_list do begin
        time_range = event.time_range
        probes = event.probes
        foreach probe, probes do begin
            pflux_grant_read_pflux, time_range, probe=probe
            pflux_grant_read_efield, time_range, probe=probe
            pflux_grant_read_bfield, time_range, probe=probe
            rbsp_read_en_spec, time_range, probe=probe
            rbsp_read_pa_spec, time_range, probe=probe
        endforeach
        stop
    endforeach
    
    
end