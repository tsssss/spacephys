
time_range = time_double(['2013-05-01/07:30','2013-05-01/09:30'])
sites = ['fsmi','gako','atha','tpas','whit','fsim']
min_elev = [2.5d,2.5,2.5,2.5,10,7.5]
sites = ['fsmi','gako','atha','tpas','fsim']
min_elev = [1,1,1,1,10]
sites = ['fsmi','gako','atha','tpas','fsim']
min_elev = [0.5,0.5,0.5,0.5,10]
merge_method = 'merge_elev'
mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, min_elev=min_elev, merge_method=merge_method)

end