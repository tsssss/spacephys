
time_range = time_double(['2015-03-17/08:30','2015-03-17/09:10'])
sites = ['whit','atha'];,'whit','fsmi','atha'];,'rank','snkq','kuuj']
min_elev = [1,2.5]
merge_method = 'max_elev'
mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, $
    min_elev=min_elev, merge_method=merge_method, calibration_method='simple')

get_data, mlt_image_var, times, mlt_images
sgopen, 0, size=[1,1]*6
stop
foreach time, times, time_id do sgtv, bytscl(reform(mlt_images[time_id,*,*]),top=254, min=0, max=5000), position=[0,0,1,1], ct=49

end