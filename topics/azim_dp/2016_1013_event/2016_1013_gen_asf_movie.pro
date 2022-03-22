
time_range = time_double(['2016-10-13/12:00','2016-10-13/13:00'])
sites = ['kian','gako','whit','mcgr']
movie_file = join_path([googledir(),'works','azim_dp','2016_1013_event','2016_1013_asf_movie_v01.mp4'])

themis_gen_asf_mltimg_circ_movie, time_range, filename=movie_file, sites=sites, $
    zrange=zrange, reset_images=1
end