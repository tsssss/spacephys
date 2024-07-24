    time_range = ['2014-03-29/09:00','2014-03-29/10:00']
    ;sites = ['inuv','fsim','fsmi']
    sites = ['inuv','fykn','fsmi','gill']
    merge_method = 'merge_elev'
    min_elev = 2.5
; fsim skymap is wrong?

    ;time_range = time_double(['2013-03-17/08:30','2013-03-17/09:30'])
    ;sites = ['chbg']

    mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, merge_method=merge_method, min_elev=min_elev)

    root_dir = srootdir()
    base_name = 'thg_asf_mlt_image_movie_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hhmm'),'_to_')+'_v01.mp4'
    movie_file = join_path([root_dir,base_name])

    mlt_range = [-1d,1]*6
    mlat_range = [55d,90]
    zrange = [0,2e3]

    themis_gen_mlt_image_movie, mlt_image_var, filename=movie_file, $
        mlt_range=mlt_range, mlat_range=mlat_range, fig_xsize=6, zrange=zrange


    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, sites=sites, merge_method=merge_method, min_elev=min_elev)


end