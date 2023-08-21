;+
; These times are from the short list in alfven_arc_survey_manual, https://docs.google.com/spreadsheets/d/1z5BdXuCVXaiTYLfpOMbgM8O5gtG9hJCqdDDPyEe9pKY/edit#gid=1293495516.
;-
test_times = time_double([$
    '2013-04-14/07:00', $
    '2013-04-14/08:00', $
    '2013-04-26/04:50', $
    '2013-04-26/10:10', $
    '2013-05-01/06:20', $
    '2013-05-01/07:35', $
    '2013-05-01/08:50', $
    '2013-05-02/02:40', $
    '2013-05-06/07:50', $
    '2013-05-12/04:10', $
    '2013-06-07/02:50', $
    '2013-06-07/05:00', $
    '2013-06-07/06:00', $
    '2013-06-10/05:10', $
    '2013-06-10/06:00', $
    '2013-06-29/06:50', $
    '2013-06-29/07:10', $
    '2013-07-15/05:50', $
    '2013-08-16/03:10', $
    '2013-08-31/02:35', $
    '2014-12-07/12:25', $
    '2014-12-12/12:35', $
    '2014-12-25/13:55', $
    '2014-12-25/23:30', $
    '2015-01-04/23:15', $
    '2015-01-05/00:40', $
    '2015-01-05/01:10', $
    '2015-01-08/01:00', $
    '2015-01-22/13:10', $
    '2015-01-22/13:40', $
    '2015-01-22/10:30', $
    '2015-01-22/11:05', $
    '2015-01-22/13:05', $
    '2015-01-26/13:30', $
    '2015-01-30/09:40', $
    '2015-02-02/10:05', $
    '2015-02-05/10:00', $
    '2015-02-17/10:00', $
    '2015-02-18/10:40', $
    '2015-02-18/11:15', $
    '2015-03-02/11:00', $
    '2015-03-06/08:00', $
    '2015-03-08/11:00', $
    '2015-03-12/08:55', $
    '2015-03-12/09:35', $
    '2015-03-13/12:35', $
    '2015-03-17/07:00', $
    '2015-03-17/07:50', $
    '2015-03-17/08:15', $
    '2015-03-19/11:40', $
    '2015-05-13/03:20', $
    '2015-05-13/04:00', $
    '2015-05-13/05:00' ])
probes = 'f'+['16','17','18','19']
pad_time = 3600*2
local_root = join_path([default_local_root(),'survey_plot','alfven_arc'])

foreach test_time, test_times do begin
    the_time_range = test_time+[-1,1]*pad_time
    date_str = time_string(test_time, tformat='YYYY_MMDD')
    plot_dir = join_path([local_root,date_str])
    foreach probe, probes do begin
        print, 'Processing '+strupcase(probe)+' ...'
        files = dmsp_gen_polar_region_survey_plot(the_time_range, probe=probe, errmsg=errmsg)
        if errmsg eq '' then begin
            foreach file, files do begin
                base = fgetbase(file)
                source_file = file
                target_file = join_path([plot_dir,base])
                if file_test(target_file) eq 1 then continue
                if file_test(plot_dir) eq 0 then file_mkdir, plot_dir
                file_copy, source_file, target_file
            endforeach
        endif
    endforeach
endforeach


stop
