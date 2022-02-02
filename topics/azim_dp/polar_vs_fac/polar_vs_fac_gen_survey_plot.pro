;+
; Generate survey plots (ewograms for FACs and MLT image) for all events.
;-

pro polar_vs_fac_gen_survey_plot

test = 1
    mlat_range = [60,70]
    mlt_range = [-1,1]*9

    fig_xsize = 5
    fig_ysize = 5

    plot_dir = join_path([homedir(),'polar_vs_fac','survey_plot'])

    event_time_ranges = polar_vs_fac_read_fac_time_ranges()
    nevent = n_elements(event_time_ranges)*0.5
    for event_id=0,nevent-1 do begin
        the_time_range = reform(event_time_ranges[event_id,*])
        themis_read_upward_current_ewo, the_time_range, mlat_range=mlat_range, mlt_range=mlt_range
        themis_read_downward_current_ewo, the_time_range, mlat_range=mlat_range, mlt_range=mlt_range
        polar_read_ewo, the_time_range, mlat_range=mlat_range, mlt_range=mlt_range

        vars = ['thg_j_up_ewo','thg_j_down_ewo','po_ewo']

        base = 'polar_vs_fac_gen_survey_plot_'+strjoin(time_string(the_time_range,tformat='YYYY_MMDD_hhmm'),'_')+'_v01.pdf'
        plot_file = join_path([plot_dir,base])
        if keyword_set(test) then plot_file = 0
        sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize
        tplot, vars, trange=the_time_range
        if keyword_set(test) then stop
        sgclose
    endfor

end