;+
; Write candidate to file.
; Adopted from azim_df_search_candidate_write_file.
; 
; candidate.
; filename=.
;-

pro azim_dp_candidate_write, candidate, filename=file

    tab = constant('4space')
    tformat = 'YYYY-MM-DD/hh:mm'

    msg = ''
    ; Time range.
    msg += strjoin(time_string(candidate.time_range,tformat=tformat),' ')+tab
    ; Duration, in hr.
    msg += string(total(candidate.time_range*[-1,1])/3600.,format='(F5.1)')+tab
    ; MLT range.
    msg += strjoin(string(candidate.mlt_range,format='(F5.1)'),' ')+tab
    ; The probes.
    msg += strjoin(candidate.probes,',')+tab
    lprmsg, msg, file

end
