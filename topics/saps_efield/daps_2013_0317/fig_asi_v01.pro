
function fig_rbsp_overview_v01, test=test

    version = 'v01'
    id = '2013_0317'

    root_dir = join_path([googledir(),'works','2024_daps',id])
    plot_dir = join_path([root_dir,'plot'])

    time_range = time_double(['2013-03-17/08:30','2013-03-17/09:30'])
    sites = ['kapu','snkq','gill','pina','fsmi','fsim']
    nsite = n_elements(sites)
    min_elevs = fltarr(nsite)+5
    index = where(sites eq 'rank', count)
    if count ne 0 then min_elevs[index] = 10
    index = where(sites eq 'fsmi', count)
    if count ne 0 then min_elevs[index] = 10
    merge_method = 'merge_elev'
    calibration_method = 'simple'
    
    mlt_image_var = themis_asf_read_mlt_image(time_range, sites=sites, $
        min_elev=min_elevs, merge_method=merge_method, calibration_method=calibration_method)
    stop

    mlt_images = get_var_data(mlt_image_var, times=times, settings=settings)
    sgopen, 0, size=[6,6]
    tpos = [0d,0,1,1]
    zrange = [-1,1]*5e3
    ct = 70
    foreach time, times, time_id do sgtv, bytscl(reform(mlt_images[time_id,*,*]), min=zrange[0], max=zrange[1], top=254), position=tpos, ct=ct

    stop


end

print, fig_rbsp_overview_v01(test=1)
end