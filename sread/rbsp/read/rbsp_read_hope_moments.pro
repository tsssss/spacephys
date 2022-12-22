;+
; Read RBSP HOPE moments.
;-

function rbsp_read_hope_moments, input_time_range, probe=probe, errmsg=errmsg, $
    species=species
    

    errmsg = ''
    retval = ''
    prefix = 'rbsp'+probe+'_'
    
    files = rbsp_load_hope_moments(input_time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval

    cdf2tplot, files
    return, cdf_vars(files[0])

end


time_range = ['2013-05-01','2013-05-02']
probe = 'b'
species = 'p'
vars = rbsp_read_hope_moments(time_range, probe=probe, species=species)
end