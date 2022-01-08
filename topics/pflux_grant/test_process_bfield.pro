;+
; Test all steps to process B field.
;-

;---Input.
    probe = 'a'
    day = time_double('2012-10-11')     ; spikes.
    day = time_double('2013-06-07')     ; storm.
    day = time_double('2012-10-01')     ; storm.
    
;    probe = 'b'
;    day = time_double('2012-12-04')    ; gap.
;    day = time_double('2012-10-02')    ; spikes.
;    day = time_double('2014-08-28')    ; spin tone.
    
    
    local_root = join_path([homedir(),'pflux_grant_test_process_bfield'])

    secofday = constant('secofday')
    prefix = 'rbsp'+probe+'_'
    
    time_range = day+[0,secofday]
    pflux_grant_read_b_gsm, time_range, probe=probe, local_root=local_root
    
    foreach time, day+[-1,1]*secofday do begin
        time_range = time+[0,secofday]
        pflux_grant_read_b_gsm, time_range, probe=probe, local_root=local_root
    endforeach
    
    time_range = day+[0,secofday]
    pflux_grant_read_bfield, time_range, probe=probe, local_root=local_root
    
    pflux_grant_survey_on_bfield, time_range, probe=probe, test=1, local_root=local_root
    

end
