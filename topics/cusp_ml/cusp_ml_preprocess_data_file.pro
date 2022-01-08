;+
; Return the data file saves the preprocessed data.
;-
function cusp_ml_preprocess_data_file

    root_dir = join_path([googledir(),'works','works','cusp_ml','data'])
    data_file = join_path([root_dir,'cusp_ml_preprocessed_data.cdf'])
    return, data_file

end
