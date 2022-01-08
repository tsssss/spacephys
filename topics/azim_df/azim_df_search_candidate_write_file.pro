
pro azim_df_search_candidate_write_file, candidate, filename=file

    tab = constant('4space')

    msg = string(candidate.id,format='(I4)')+tab+$
        strjoin(time_string(candidate.time_range,tformat='YYYY-MM-DD/hh:mm'),' ')+tab+$
        candidate.region+tab+$
        candidate.search_name+tab+$
        string(candidate.duration,format='(I4)')+tab+$
        string(candidate.nsection,format='(I2)')+tab+$
        strjoin(candidate.all_probes,',')
    lprmsg, msg, file
    foreach time_range, candidate.time_range_list, jj do begin
        available_probes = candidate.probe_list[jj]
        msg = tab+string(jj+1,format='(I0)')+tab+$
            strjoin(time_string(time_range,tformat='YYYY-MM-DD/hh:mm'),' ')+tab+$
            strjoin(available_probes,',')
            lprmsg, msg, file
    endforeach

end