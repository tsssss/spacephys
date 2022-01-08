;+
; Check the HOPE L2 data to see if we can figure out the direction of the azimuthal flow.
;-

;---Load data and settings.
    test = 1

    time_range = time_string(['2015-02-18/02:05','2015-02-18/02:15'])
    time_range = time_double(['2015-02-18','2015-02-19'])
    probe = 'a'
    prefix = 'rbsp'+probe+'_'



;---Check angle between each pixel and west.
    rbsp_read_orbit, time_range, probe=probe
    r_gse = get_var_data(prefix+'r_gse', times=times)
    r_sm = cotran(r_gse, times, 'gse2sm')
    r_hat = sunitvec(r_sm)      ; radially outward.
    ntime = n_elements(times)
    ndim = 3
    z_hat = fltarr(ntime,ndim)  ; z-axis.
    z_hat[*,2] = 1
    w_hat = vec_cross(r_hat, z_hat) ; westward, w = r x z.
    w_hat = sunitvec(w_hat)
    west_sm = w_hat
    west_uvw = cotran(west_sm, times, 'sm2uvw', probe=probe)

    rbsp_read_sc_vel, time_range, probe=probe
    v_gse = get_var_data(prefix+'v_gse', times=times)
    v_uvw = cotran(v_gse, times, 'gse2uvw', probe=probe)

    npixel = 5
    pixels = findgen(npixel)
    pixel_polar_angles = [72d,36,0,-36,-72] ; deg.
    rad = constant('rad')
    foreach pixel_id, pixels do begin
        pixel_uvw = fltarr(ntime,3)
        pixel_polar_angle = pixel_polar_angles[pixel_id]*rad
        pixel_uvw[*,0] = cos(pixel_polar_angle)
        pixel_uvw[*,2] = sin(pixel_polar_angle)

        angle_west = sang(pixel_uvw, west_uvw, degree=1)
        pixel_str = string(pixel_id,format='(I0)')
        var = 'angle_west_pixel'+pixel_str
        store_data, var, times, angle_west, limits={$
            ytitle:'deg', labels:'Pixel '+pixel_str, yrange:[0,180] }

        angle_v = sang(pixel_uvw, v_uvw, degree=1)
        var = 'angle_v_pixel'+pixel_str
        store_data, var, times, angle_v, limits={$
            ytitle:'deg', labels:'Pixel '+pixel_str, yrange:[0,180] }
    endforeach

    w_uvw = fltarr(ntime,3)
    w_uvw[*,2] = 1
    w_sm = cotran(w_uvw, times, 'uvw2sm', probe=probe)
    store_data, 'w_sm', times, w_sm, limits={$
        labels:constant('xyz'), colors:constant('rgb'), ytitle:'W in SM'}


;---Load HOPE L2.
    file = '/Volumes/Research/data/rbsp/rbspa/hope/level2/sectors_rel04/2015/rbspa_rel04_ect-hope-sci-l2_20150218_v6.1.0.cdf'

    ; [ntime,nen_bin,nsector,npixel].
    fpdu = cdf_read_var('FPDU', filename=file)
    times = convert_time(cdf_read_var('Epoch_Ion', filename=file), from='epoch', to='unix')
    en_bins = cdf_read_var('HOPE_ENERGY_Ion', filename=file)

    ntime = n_elements(times)
    nen_bin = n_elements(en_bins[0,*])


    ; Time of "HOPE" spin.
    time_step = 11.362  ; sec.
    times = times-time_step*0.5

    ; Expand energy-sum flux per sector to 1D flux.
    nsector = 16
    sec_time_step = time_step/nsector
    fillval = !values.f_nan
    foreach pixel_id, findgen(npixel) do begin
        uts = []
        flux = []
        foreach time, times, time_id do begin
            the_flux = reform(fpdu[time_id,*,*,pixel_id])
            the_uts = time+sec_time_step*(make_bins([-1,nsector],1)+0.5)
            the_flux = [fillval,total(the_flux,1),fillval]
            uts = [uts,the_uts]
            flux = [flux,the_flux]
        endforeach
        store_data, 'fpdu_flux_pixel'+string(pixel_id,format='(I0)'), $
            uts, flux, limits={ylog:1}
    endforeach


    ; Spec of energy-sum flux as a function of time and sector.
    foreach pixel_id, findgen(npixel) do begin
        spec = reform(total(fpdu[*,*,*,pixel_id],2))
        pixel_str = string(pixel_id,format='(I0)')
        store_data, 'fpdu_pixel'+pixel_str, $
            times, spec, findgen(nsector), $
            limits={spec:1, no_interp:1, ylog:0, zlog:1, zrange:[1e6,1e10], $
            ytitle:'Sector #', ztitle:'H+ flux', yticklen:-0.01, xticklen:-0.05}
    endforeach


    ; Plot the spectrogram of flux vs sector and time, for each pixel.
    vars = 'fpdu_pixel'+string(findgen(npixel),format='(I0)')
    nvar = n_elements(vars)
    sgopen, 0, xsize=6, ysize=6
    margins = [12,4,10,2]
    test_time_range = time_double(['2015-02-18/02:00','2015-02-18/02:20'])
    poss = sgcalcpos(nvar,margins=margins, xchsz=xchsz, ychsz=ychsz)
    tplot, vars, position=poss, trange=test_time_range
    fig_labels = letters(nvar)
    for ii=0,nvar-1 do begin
        tpos = poss[*,ii]
        msg = fig_labels[ii]+'. Pixel '+string(ii+1,format='(I0)')
        tx = xchsz*2
        ty = tpos[3]-ychsz*0.7
        xyouts, tx,ty,/normal, msg
    endfor
    stop





    test_time_range = time_double(['2015-02-18/02:00','2015-02-18/02:20'])
;    test_time_range = time_double(['2015-02-18/02:14','2015-02-18/02:18'])
    test_pixels = string(findgen(npixel),format='(I0)')
;    test_pixels = '1'
    foreach pixel_str, test_pixels do begin
        flux_var = 'fpdu_flux_pixel'+pixel_str
        angle_var = 'angle_west_pixel'+pixel_str
        flux = get_var_data(flux_var, times=times, in=test_time_range)
        angle = get_var_data(angle_var, at=times)
        plot, angle, flux, psym=1, $
            ;ylog=1, yrange=[1e3,1e13], ystyle=1, $
            xlog=0, xrange=[0,180], xstyle=1
        ;tplot, ['fpdu_flux_pixel','angle_west_pixel','fpdu_pixel']+pixel_str, trange=test_time_range
        stop
    endforeach

    stop

    ; Select the center pixel.
    energy_index = 40
    pixel_index = 3
    ; [ntime,nen_bin,nsector,npixel].
    fpdu = cdf_read_var('FPDU', filename=file)
    times = convert_time(cdf_read_var('Epoch_Ion', filename=file), from='epoch', to='unix')
    en_bins = cdf_read_var('HOPE_ENERGY_Ion', filename=file)
    spec = reform(fpdu[*,*,*,pixel_index])
    spec = total(fpdu[*,*,*,1:3],4)
;    for ii=0,15 do begin
;        store_data, 'fpdu_'+string(ii,format='(I02)'), times, spec[*,*,ii], en_bins, $
;            limits={spec:1, no_interp:1, ylog:1, zlog:1}
;    endfor

    ntime = n_elements(times)
    nen_bin = n_elements(en_bins[0,*])
    nsector = 16
    data = fltarr(ntime*(nsector+2), nen_bin)
    uts = dblarr(ntime*(nsector+2))
    time_step = sdatarate(times)
    time_step = 11.*1
    for ii=0ull,ntime-1 do begin
        i0 = ii*(nsector+2)
        data[i0+1:i0+nsector,*] = transpose(reform(spec[ii,*,*]))
        data[i0,*] = fillval
        data[i0+nsector+1,*] = fillval
        uts[i0:i0+nsector-1] = times[ii]+time_step*findgen(nsector)/nsector
        uts[i0] = times[ii]-time_step/nsector
        uts[i0+nsector+1] = times[ii]+time_step+time_step/nsector
    endfor
    var = 'fpdu_sec'
    index = where(uts ne 0)
    uts = uts[index]
    data = data[index,*]
    val = reform(en_bins[0,*])
    store_data, var, uts, data, val, $
        limits={spec:1, no_interp:1, ylog:1, zlog:1}

    zlim, var, 1e4, 1e9
    ylim, var, minmax(val)
    tplot, var







end
