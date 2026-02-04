function micro_injection_get_event_time_range_file, time_range, mission_probe=mission_probe

    mission = mission_probe[0]
    probe = mission_probe[1]

    project_id = 'micro_injection'
    project = project_load(project_id)
    data_dir = join_path([project['data_dir'],'statistical_study'])

    ; This is time ranges when MMS were in [5,13] Re.
    orbit_trs = micro_injection_load_survey_time_range()
    norbit = n_elements(orbit_trs[*,0])
    probes = ['1','2','3','4']
    probes = ['1','2','3']

    foreach probe, probes do begin
        for ii=0,norbit-1 do begin
            orbit_tr = reform(orbit_trs[ii,*])
            var = mms_read_kev_electron(orbit_tr, probe=probe)
        endfor
    endforeach

stop
end



function micro_injection_get_event_time_range, time_range, mission_probe=mission_probe

    mission = mission_probe[0]
    probe = mission_probe[1]
    files = micro_injection_get_event_time_range_file(time_range, mission_probe=mission_probe)


end


tr = ['2015-09-01','2015-11-01']
event_trs = micro_injection_get_event_time_range(tr, mission_probe=['mms','4'])
end