;+
; Write roi_list to file.
; Adopted from azim_df_search_candidate_write_file.
;-

pro azim_dp_roi_write, roi, filename=file

    tab = constant('4space')
    tformat = 'YYYY-MM-DD/hh:mm'

    msg = ''
    ; Time range.
    msg += strjoin(time_string(roi.time_range,tformat=tformat),' ')+tab
    ; Duration, in hr.
    msg += string(total(roi.time_range*[-1,1])/3600.,format='(F5.1)')+tab
    ; # of sections.
    msg += string(roi.nsection,format='(I2)')+tab
    ; The probes.
    msg += strjoin(roi.all_probes,',')+tab
    lprmsg, msg, file

    foreach time_range, roi.time_range_list, roi_id do begin
        probes = roi.probe_list[roi_id]
        msg = ''
        msg += tab+string(roi_id+1,format='(I2)')
        msg += tab+strjoin(time_string(time_range,tformat=tformat),' ')
        msg += tab+strjoin(probes,',')
        lprmsg, msg, file
    endforeach

end
