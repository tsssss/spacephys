;+
; Test to bin Polar orbit.
;-

bin_info = dictionary($
    'resolution', 60, $     ; in sec.
    'based_on', 'r_gsm', $
    'generate_method', 'binning_prepare_rgsm', $    ; call this procedure to generate the orbital quantity for binning.
    'components', ['x','y','z']+'_gsm')

bin_info['x_gsm'] = dictionary($
    'boundary', make_bins([-10,10], 1))
