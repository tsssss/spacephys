;+
; Calculate E/B ratio.
;-
function stplot_calc_ebratio, input_time_range, e_var=e_var, b_var=b_var, scale_info=scale_info, output=ebr_var

    foreach var, [e_var,b_var] do begin
        get_data, var, times, vec
        dims = size(vec,dimensions=1)
        ndim = n_elements(dims)
        if ndim eq 2 then begin
            dat = snorm(vec)
        endif else begin
            dat = temporary(vec)
        endelse
        index = where(finite(dat,nan=1),count)
        if count ne 0 then dat[index] = 0

        if n_elements(input_time_range) ne 0 then begin
            time_range = time_double(input_time_range)
        endif else begin
            time_range = minmax(times)
        endelse

        dr0 = sdatarate(times)
        dj = 1d/8
        s0 = 4d*dr0
        s1 = 1800
        if n_elements(scale_info) ne 0 then begin
            if size(scale_info,type=1) eq 8 then scale_info = dictionary(scale_info)
            if scale_info.haskey('s0') then s0 = scale_info['s0']
            if scale_info.haskey('s1') then s1 = scale_info['s1']
            if scale_info.haskey('dj') then dj = scale_info['dj']
        endif

        j1 = floor(alog(s1/s0)/alog(2)/dj)
        s1 = s0*2d^(dj*j1)
        ns = j1+1
        w0 = 6d
        cdelta = 0.776d
        psi0 = !dpi^(-0.25)

        mor = wavelet(dat, dr0, pad=1, $
            s0=s0, dj=dj, j=j1, mother='Morlet', param=w0, period=ps, scale=ss)
        mor_var = var+'_mor'
        fs = 1d/ps
        store_data, mor_var, times, mor, fs
        add_setting, mor_var, smart=1, dictionary($
            'display_type', 'spec', $
            'short_name', 'Morlet '+get_setting(var, 'short_name'), $
            'unit', '('+get_setting(var, 'unit')+')!U2!N' )
        
        index = where(times ge time_range[0] and times le time_range[1], tnrec)
        psd = abs(mor)^2
        psd = psd[index,*]
        gws = total(psd,1)/tnrec^2
        ngws = (gws/ss)*(dr0*dj/cdelta)*tnrec

        options, mor_var, 'ps', ps      ; period in sec.
        options, mor_var, 'fs', fs      ; frequency in Hz.
        options, mor_var, 'gws', gws
        options, mor_var, 'ngws', ngws
    endforeach

    e_mor_var = e_var+'_mor'
    b_mor_var = b_var+'_mor'
    e_gws = get_setting(e_mor_var, 'gws')
    b_gws = get_setting(b_mor_var, 'gws')
    ebr = sqrt(e_gws/b_gws)*1e3

    if n_elements(ebr_var) eq 0 then begin
        prefix = get_prefix(e_var)
        ebr_var = prefix+'ebratio'
    endif
    store_data, ebr_var, ps, ebr
    return, ebr_var

end