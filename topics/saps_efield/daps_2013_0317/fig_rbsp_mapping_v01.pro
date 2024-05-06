
function fig_rbsp_mapping_v01, test=test, event_info=event_info

    version = 'v01'
    id = '2013_0317'

    if n_elements(event_info) eq 0 then event_info = saps_efield_load_data(id)

    mlt_image_var = 'thg_asf_mlt_image_rect'
    mlt_images = get_var_data(mlt_image_var, times=times, settings=settings)
    ntime = n_elements(times)
    
    
    ; KEO
    keo_var = 'thg_asf_keo'
    ytitle = 'MLat!C(deg)'
    
    mlat_range = [55,75]
    mlt_range = [1,5]
    mlt_bins = settings['mlt_bins']
    mlat_bins = settings['mlat_bins']
    index = where_pro(mlt_bins, '[]', mlt_range, count=nbin)
    ;bins = mlt_bins[index]
    ;keo = fltarr(ntime,nbin)
    keo = total(mlt_images[*,index,*],2)/nbin

    
    fmlt_var = 'rbspa_fmlt_igrf_t01_north'
    fmlt = get_var_data(fmlt_var, at=times)
    fmlt_del = [-2,0.5]
    fmlt_del = [-0.1,0.1]
    nbin = n_elements(mlat_bins)
    keo = fltarr(ntime,nbin)
    foreach time, times, time_id do begin
        index = where_pro(mlt_bins,'[]',fmlt[time_id]+fmlt_del,count=count)
        if count eq 0 then continue
        keo[time_id,*] = total(mlt_images[time_id,index,*],2)/count
    endforeach
    
    store_data, keo_var, times, keo, mlat_bins, limits={ytitle:ytitle,spec:1,zlog:0,color_table:49}
    
    
    ; EWO
    ewo_var = 'thg_asf_ewo_omega_band'
    ytitle = 'MLT!C(h)'
    
    mlt_range = [0,6]
    mlat_range = [60,66]
    index = where_pro(mlat_bins, '[]', mlat_range, count=nbin)
    ewo = total(mlt_images[*,*,index],3)/nbin
    
    store_data, ewo_var, times, ewo, mlt_bins, limits={ytitle:ytitle,yrange:mlt_range,spec:1,zlog:0,color_table:49}
    
    ewo_var = 'thg_asf_ewo_streamer'
    ytitle = 'MLT!C(h)'

    mlt_range = [0,6]
    mlat_range = [68,72]
    index = where_pro(mlat_bins, '[]', mlat_range, count=nbin)
    ewo = total(mlt_images[*,*,index],3)/nbin

    store_data, ewo_var, times, ewo, mlt_bins, limits={ytitle:ytitle,yrange:mlt_range,spec:1,zlog:0,color_table:49}
    
    ewo_vars = 'thg_asf_ewo_'+['omega_band','streamer']
    options, ewo_vars, 'zrange', [0,5000]
    
    keo_var = 'thg_asf_keo'
    options, keo_var, 'zrange', [0,5000]
    
    fig_size = [5,6]
    sgopen, 0, size=fig_size, xchsz=xchsz, ychsz=ychsz
    plot_vars = ['thg_asf_ewo_'+['omega_band','streamer'],'thg_asf_keo','rbspa_'+['edot0_fac','e_spec']]
    nplot_var = n_elements(plot_vars)
    margins = [12,4,10,1]
    poss = sgcalcpos(nplot_var, margins=margins)
    plot_tr = time_double(['2013-03-17/07:30','2013-03-17/09:30'])
    plot_tr = time_double(['2013-03-17/05:30','2013-03-17/09:30'])
    tplot, plot_vars, trange=plot_tr, position=poss
    
    yrange = mlt_range
    foreach var, 'thg_asf_ewo_'+['omega_band','streamer'] do begin
        pid = where(plot_vars eq var, count)
        if count ne 0 then begin
            set_axis, var, position=poss[*,pid], xrange=plot_tr, yrange=yrange
            fmlt = get_var_data(fmlt_var, times=times, in=plot_tr)
            plots, times, fmlt, color=sgcolor('red')
        endif
    endforeach
    
    var = 'thg_asf_keo'
    yrange = [55,75]
    pid = where(plot_vars eq var, count)
    if count ne 0 then begin
        set_axis, var, position=poss[*,pid], xrange=plot_tr, yrange=yrange
        
        fmlat_vars = list()
        external_models = ['t89','t01']
        ;external_models = ['t01']
        foreach external_model, external_models do begin
            fmlat_vars.add, 'rbspa_fmlat_'+['dipole']+'_'+external_model+'_north'
        endforeach
        fmlat_vars = fmlat_vars.toarray()
        nfmlat_var = n_elements(fmlat_vars)
        colors = get_color(nfmlat_var)
        
        foreach fmlat_var, fmlat_vars, tid do begin
            fmlat = get_var_data(fmlat_var, times=times, in=plot_tr)+5
            plots, times, fmlat, color=colors[tid]
            tx = times[0]
            ty = fmlat[0]
            tmp = convert_coord(times[0],fmlat[0], data=1, to_normal=1)
            tx = tmp[0]
            ty = tmp[1]+ychsz*0.2
            msg = strupcase(external_models[tid])
            xyouts, tx,ty,msg, normal=1, color=colors[tid]
        endforeach
        plots, minmax(times), 61+[0,0], linestyle=2, sgcolor('red')
    endif

    bar_times = time_double(['2013-03-17/08:55:21','2013-03-17/09:05:39'])
    timebar, bar_times, linestyle=1, color=sgcolor('red')
stop


end

print, fig_rbsp_mapping_v01(test=1, event_info=event_info)
end