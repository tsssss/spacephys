;+
; Preprocess the features and prediction.
;-

pro cusp_mlt_polar_test_preprocess_data


;---Check data file.
    root_dir = join_path([googledir(),'works','cusp_ml','data'])
    file = join_path([root_dir,'polar_test_data.cdf'])


;---Load data.
    routines = 'cusp_ml_polar_test_preprocess_'+['density','orbit','dst']
    foreach routine, routines do call_procedure, routine


;---Remove bad data.
    get_data, 'polar_density', common_times, dens
    good_index = where(finite(dens))
    
    new_names = 'polar_r'+['x','y','z']+'_sm'
    stplot_split, 'polar_r_sm', newnames=new_names
    
    vars = ['dst','polar_'+['r'+['x','y','z']+'_sm','density']]
    foreach var, vars do begin
        get_data, var, common_times, data
        common_times = common_times[good_index]
        data = data[good_index]
        store_data, var, common_times, data
    endforeach


;---Save data.
    time_var = 'time'
    stplot2cdf, vars, filename=file, time_var=time_var, istp=1

end



cusp_mlt_polar_test_preprocess_data
end