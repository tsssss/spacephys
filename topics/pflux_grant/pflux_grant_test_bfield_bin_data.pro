;+
; Test B field quality by binning data in spatial grids.
;-

pro pflux_grant_test_bfield_bin_data


    ; empty bin.
    mlt_range = [-10]+[-1,1]*0.5
    mlat_range = -20+[0,5]
    dis_range = 4.0+[0,0.5]

    ; weird bin.
    mlt_range = [8]+[-1,1]*0.5
    mlat_range = -20+[0,5]
    dis_range = 5.5+[0,0.5]
    test_bin = list()
    test_bin.add, mlt_range
    test_bin.add, mlat_range
    test_bin.add, dis_range

    var_type = 'b0_sm'
    probes = ['a']
    foreach probe, probes do pflux_survey_bin_data, var_type, probe=probe;, test_bin=test_bin
    
end