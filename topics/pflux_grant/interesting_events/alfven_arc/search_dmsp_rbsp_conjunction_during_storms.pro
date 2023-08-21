;+
; Generate survey plots during storms
;-

function read_storm_time_ranges

    file = join_path([srootdir(),'data','storm_list_v1.txt'])
    
    if file_test(file) eq 0 then message, 'Inconsistency ...'
    lines = read_all_lines(file)
    
    ntime = n_elements(lines)
    storm_time_ranges = dblarr(ntime,2)
    for ii=0,ntime-1 do begin
        info = strsplit(lines[ii],' ',extract=1)
        storm_time_ranges[ii,*] = time_double(info)
    endfor
    return, storm_time_ranges
end

function search_dmsp_rbsp_conjunction_during_storms, plot_dir=plot_dir, test=test

    if n_elements(plot_dir) eq 0 then plot_dir = join_path([srootdir(),'dmsp_rbsp_during_storms'])
    storm_time_ranges = read_storm_time_ranges()
    nstorm = n_elements(storm_time_ranges[*,0])

    dmsp_probes = 'f'+['16','17','18','19']
    rbsp_probes = ['a','b']
    valid_time_range = time_double(['2012-11-01','2015-10-01'])
    for ii=0,nstorm-1 do begin
        time_range = reform(storm_time_ranges[ii,*])
        if time_range[0] lt valid_time_range[0] then continue
        if time_range[1] gt valid_time_range[1] then continue
        storm_id = time_string(time_range[0],tformat='YYYY_MMDD')
        the_plot_dir = join_path([plot_dir,'rbsp_survey_plots'])
        all_files = list()
        foreach rbsp_probe, rbsp_probes do begin
            files = rbsp_gen_polar_region_survey_plot(time_range, probe=rbsp_probe, plot_dir=the_plot_dir, test=test)
            if n_elements(files) eq 0 then continue
            all_files.add, files, extract=1
        endforeach
        all_files = all_files.toarray()
    endfor

    valid_time_range = time_double(['2015-06-29','2015-10-01'])
    valid_time_range = time_double(['2015-08-24','2015-10-01'])
    for ii=0,nstorm-1 do begin
        time_range = reform(storm_time_ranges[ii,*])
        if time_range[0] lt valid_time_range[0] then continue
        if time_range[1] gt valid_time_range[1] then continue
        storm_id = time_string(time_range[0],tformat='YYYY_MMDD')
        the_plot_dir = join_path([plot_dir,storm_id])
        all_files = list()
        foreach dmsp_probe, dmsp_probes do begin
            files = dmsp_gen_polar_region_survey_plot(time_range, probe=dmsp_probe, plot_dir=the_plot_dir, test=test)
            if n_elements(files) eq 0 then continue
            all_files.add, files, extract=1
        endforeach
        foreach rbsp_probe, rbsp_probes do begin
            files = rbsp_gen_polar_region_survey_plot(time_range, probe=rbsp_probe, plot_dir=the_plot_dir, test=test)
            if n_elements(files) eq 0 then continue
            all_files.add, files, extract=1
        endforeach
        all_files = all_files.toarray()
        files = file_search(join_path([plot_dir,storm_id]),'*.pdf')
        foreach file, files do begin
            index = where(all_files eq file, count)
            if count ne 0 then continue
            file_delete, file
        endforeach
    endfor

    return, 'Done'

end


tmp = search_dmsp_rbsp_conjunction_during_storms()
end