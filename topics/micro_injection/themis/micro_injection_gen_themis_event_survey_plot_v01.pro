;+
; Survey plot for Themis events. Save in different folder.
;-

    event_list = list()
;    event_list.add, dictionary($
;        'tr', time_double('2009-01-07/'+['07:00','13:00']), $
;        'probes', ['a','d','e','c'] )
;    event_list.add, dictionary($
;        'tr', time_double('2009-06-19/'+['09:00','12:00']), $
;        'probes', ['a','d','e'] )
    event_list.add, dictionary($
        'tr', time_double('2009-06-19/'+['00:00','06:00']), $
        'probes', ['b','c','d','e'] )
        
    foreach event, event_list do begin
        tr = event['tr']
        probes = event['probes']
        plot_dir = join_path([googledir(),'works','2024_microinjection','plot','themis_event'])
        
        foreach probe, probes do begin
            print, 'Processing '+probe+' for '+time_string(time_range[0])+' ...'
            
            base = 'themis_microinjection_survey_plot_'+strjoin(time_string(tr,tformat='YYYY_MMDD_hh'),'_')+'_th'+probe+'_v01.pdf'
            plot_file = join_path([plot_dir,time_string(tr[0],tformat='YYYY_MMDD_hh'),base])
            ;if file_test(plot_file) eq 1 then continue
            del_data, '*'
            plot_file = micro_injection_gen_themis_survey_plot_v01(tr, probe=probe, test=0, filename=plot_file)
            print, plot_file
        endforeach
    endforeach
end
