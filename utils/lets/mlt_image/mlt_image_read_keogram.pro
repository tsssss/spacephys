;+
; Read keogram from MLT image.
;-

function mlt_image_read_keogram, mlt_image_var, $
    mlt_range=mlt_range, mlat_range=mlat_range, mlat_step=bin_step, $
    mlat_bins=bins, errmsg=errmsg, var_info=var_info

    errmsg = ''
    retval = !null
    
    if n_elements(var_info) eq 0 then var_info = mlt_image_var+'_keo'

    if tnames(mlt_image_var) eq '' then begin
        errmsg = 'MLT image variable not defined.'
        return, retval
    endif

    mlt_images = get_var_data(mlt_image_var, times=times, settings=settings)
    ntime = n_elements(times)
    pixel_mlat = settings.pixel_mlat
    pixel_mlt = settings.pixel_mlt
    unit = settings.unit
    
    
    if n_elements(bins) eq 0 then begin
        bin_range = minmax(pixel_mlat)
        if n_elements(mlat_range) eq 2 then bin_range = mlat_range
        if n_elements(bin_step) eq 0 then bin_step = 0.5
        bins = make_bins(bin_range, bin_step)
    endif
    nbin = n_elements(bins)
    
    if n_elements(mlt_range) eq 0 then mlt_range = [-1,1]*12d
    index = where(mlt_range ge 12, count)
    if count ne 0 then mlt_range[index] -= 24   ; convert to -12,12.
    wanted_index = where_pro(pixel_mlt, '[]', mlt_range, count=count)
    if count eq 0 then begin
        errmsg = 'No pixel in the given MLT range ...'
        return, retval
    endif
    
    keos = fltarr(ntime,nbin-1)
    ;counts = fltarr(ntime,nbin-1)
    wanted_pixel_mlat = pixel_mlat[wanted_index]
    for tid=0,ntime-1 do begin
        the_image = reform(mlt_images[tid,*,*])
        the_image = the_image[wanted_index]
        for bid=0,nbin-2 do begin
            the_range = bins[bid:bid+1]
            index = where_pro(wanted_pixel_mlat, '[]', the_range, count=count)
            ;counts[tid,bid] = count
            if count eq 0 then continue
            keos[tid,bid] = mean(the_image[index],nan=1)
        endfor
    endfor
    bin_centers = (bins[0:nbin-2]+bins[1:nbin-1])*0.5
    store_data, var_info, times, keos, bin_centers
    add_setting, var_info, smart=1, dictionary($
        'display_type', 'spec', $
        'yrange', bin_range, $
        'unit', unit, $
        'ytitle', 'MLat (deg)' )

    return, var_info
    
end