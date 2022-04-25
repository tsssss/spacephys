;+
; Detect bad bfield data in the combined_file.
;
; Check for:
;   1. Missing data point, should be 1 min cadence.
;   2. Check for discontinuities, through checking large velocity.
;-

pro global_efield_detect_and_fix_bad_bfield_data, mission_probe, project=project

    global_efield_load_data, 'b_gsm', probe=mission_probe, project=project
    global_efield_load_data, 'bmod_t89_gsm', probe=mission_probe, project=project
    global_efield_load_data, 'r_gsm', probe=mission_probe, project=project

    prefix = project[mission_probe].prefix
    var = prefix+'b_gsm'

    lprmsg, 'Detecting bad data in '+var+' ...'
    ; Calculate the difference and ratio.
    get_data, var, times, bgsm
    get_data, prefix+'bmod_t89_gsm', times, bmodgsm
    rbmag = snorm(bgsm)/snorm(bmodgsm)
    store_data, prefix+'rbmag', times, rbmag
    dbgsm = bgsm-bmodgsm
    store_data, prefix+'db_gsm', times, dbgsm
    add_setting, prefix+'db_gsm', /smart, {$
        display_type: 'vector', $
        short_name: 'dB', $
        unit: 'nT', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}

    get_data, prefix+'r_gsm', times, rgsm
    dis = snorm(rgsm)
    store_data, prefix+'dis', times, dis

    r_psphere = 4.
    index = where(dis le 4 and abs(rbmag-1) gt 0.5, count)
    stop



end

global_efield_detect_and_fix_bad_bfield_data, 'rbspa', project=project
end