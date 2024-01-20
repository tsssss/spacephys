;+
; Get and calculate the spatially averaged energy fluxes.
;-

function fig_2015_0416_0800_eflux_v01, event_info=event_info

test = 1

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
    base_name = 'fig_'+id+'_outflow_'+version+'.pdf'
    
    
    fig_size = [6,6]
    
    sgopen, 0, size=fig_size, xchsz=xchsz, ychsz=ychsz
    uniform_ticklen = -ychsz*0.15*fig_size[0]    
    
    

;---RBSP settings.
    rbsp_info = event_info.rbsp.rbspa
    time_range = event_info.time_range
    time_range = time_double(['2015-04-16/08:05','2015-04-16/08:40'])
    prefix = rbsp_info['prefix']
    probe = rbsp_info['probe']
    
    
    rbspa_o_eflux_var = prefix+'o_eflux_map'

    species = 'o'
    files = rbsp_load_hope(time_range, id='l3%pa', probe=probe, errmsg=errmsg)

    var_list = list()
    suffix = (species eq 'e')? '_Ele': '_Ion'
    time_var = 'Epoch'+suffix
    energy_var = 'HOPE_ENERGY'+suffix
    flux_var = strupcase('f'+species+'du')
    var_list.add, dictionary($
        'in_vars', [energy_var,flux_var], $
        'time_var_name', time_var, $
        'time_var_type', 'Epoch' )

    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    pitch_angles = cdf_read_var('PITCH_ANGLE', filename=files[0])
    npitch_angle = n_elements(pitch_angles)
    get_data, flux_var, times, fluxs
    fillval = 0
    index = where(abs(fluxs) ge 1e30 or fluxs eq 0, count)
    if count ne 0 then fluxs[index] = fillval
    ; flux unit is '#/cm!E2!N-s-sr-keV'
    energys = get_var_data(energy_var)
    ; energy unit is eV.
    energys = reform(energys[0,*])
    ;pa_range = [135,180]
    pa_range = [0d,30]
    en_range = [1e2,3e4]
    en_index = where_pro(energys,'[]', en_range, count=nen)
    pa_index = where_pro(pitch_angles,'[]', pa_range, count=npa)
    o_eflux = fltarr(n_elements(times))
    for ii=0,nen-1 do begin
        for jj=0,npa-1 do begin
            dtheta = 18*constant('rad')
            sr = 2*!dpi*dtheta*sin(pitch_angles[jj]*constant('rad'))
            o_eflux += fluxs[*,en_index[ii],pa_index[jj]]*(energys[en_index[ii]])^2*1e-3*sr
        endfor
    endfor
    o_eflux *= 1.6e-19*1e4*1e3  ; convert from eV/cm^2-s to mW/m^2.
    store_data, rbspa_o_eflux_var, times, o_eflux
    
    external_model = 't04s'
    internal_model = 'dipole'
    bf_var = prefix+'bf_gsm_'+external_model+'_'+internal_model+'_south'
    b0_var = prefix+'b0_gsm'
    b0_gsm = get_var_data(b0_var, at=times)
    bf_gsm = get_var_data(bf_var, at=times)
    cmap = snorm(bf_gsm)/snorm(b0_gsm)
    o_eflux_map = o_eflux*cmap
    store_data, rbspa_o_eflux_var, times, o_eflux_map
    
    
    
    

;---DMSP settings.
    dmsp_info = event_info.dmsp.dmspf19
    dmsp_color = dmsp_info.sc_color
    prefix = dmsp_info['prefix']
    probe = dmsp_info['probe']
    dmsp_arc_time_range = time_double(['2015-04-16/08:03:32','2015-04-16/08:03:42'])
    
    ssusi_id = 'energy'
    if ssusi_id eq '1356' then begin
        ssusi_wavelength = '135.6 nm'
    endif else if ssusi_id eq '1216' then begin
        ssusi_wavelength = '121.6 nm'
    endif else if ssusi_id eq 'energy' then begin
        ssusi_wavelength = 'energy'
    endif else ssusi_wavelength = strupcase(ssusi_id)

    dmsp_plot_time_range = time_double(['2015-04-16/08:00','2015-04-16/08:07'])
    dmsp_mlt_image_var = dmsp_read_mlt_image(dmsp_plot_time_range+[-1800,0], probe=probe, id=ssusi_id)
;    mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id='1216')
    mlat_vars = dmsp_read_mlat_vars(dmsp_plot_time_range, probe=probe, errmsg=errmsg)
    ele_spec_var = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(dmsp_plot_time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(dmsp_plot_time_range, probe=probe, species='e', errmsg=errmsg)
    db_xyz_var = dmsp_read_bfield_madrigal(dmsp_plot_time_range, probe=probe)
    r_var = dmsp_read_orbit(dmsp_plot_time_range, probe=probe)

    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]


    ; Get ssusi eflux using pixel location.
    min_mlat = 50.
    mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
    mlons = get_var_data(mlon_var, at=the_times)
    ; use aacgm mlon.
    mlts = aacgm_mlon2mlt(mlons, the_times)
    sc_tt = (mlts*15-90)*constant('rad')
    sc_rr = (90-abs(mlats))/(90-min_mlat)
    sc_xx = sc_rr*cos(sc_tt)
    sc_yy = sc_rr*sin(sc_tt)
    ntime = n_elements(the_times)
    ssusi_eflux = fltarr(ntime)

    get_data, dmsp_mlt_image_var, times, mlt_images, limits=lim
    tmp = min(times-mean(dmsp_plot_time_range), abs=1, time_id)
    mlt_image = reform(mlt_images[time_id,*,*])
    pixel_mlat = lim.pixel_mlat
    pixel_mlt = lim.pixel_mlt

    min_mlat = 50
    pixel_tt = (pixel_mlt*15-90)*constant('rad')
    pixel_rr = (90-pixel_mlat)/(90-min_mlat)
    pixel_xx = pixel_rr*cos(pixel_tt)
    pixel_yy = pixel_rr*sin(pixel_tt)

    mlt_image = mlt_image[*]
    for ii=0,ntime-1 do begin
        dis = sqrt((pixel_xx-sc_xx[ii])^2+(pixel_yy-sc_yy[ii])^2)
        tmp = min(dis[*], abs=1, index)
        ssusi_eflux[ii] = mlt_image[index]
    endfor
    ssusi_eflux_var = prefix+'ssusi_eflux'
    index = where(ssusi_eflux[1:-1]-ssusi_eflux[0:-2] ne 0)
    store_data, ssusi_eflux_var, the_times[index+1], ssusi_eflux[index+1]
    dmsp_eflux_yrange = [1e-1,1e2]
    add_setting, ssusi_eflux_var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', dmsp_eflux_yrange, $
        'display_type', 'scalar', $
        'short_name', 'Aurora eflux', $
        'unit', 'mW/m!U2!N' )


    ; Get electron eflux, map eflux and convert unit.
    external_model = 't89'
    internal_model = 'dipole'
    bmod_var = geopack_read_bfield(r_var=r_var, models=external_model, igrf=0, suffix='_'+internal_model, t89_par=2, coord='gsm')
    bfmod_var = prefix+'bf_gsm_t89_dipole'
    if check_if_update(bfmod_var) then begin
        vinfo = geopack_trace_to_ionosphere(r_var, models=external_model, $
            igrf=0, south=1, north=0, refine=1, suffix='_'+internal_model)
    endif
    cmap = snorm(get_var_data(bfmod_var, times=times))/snorm(get_var_data(bmod_var, at=times))
    cmap_var = prefix+'cmap'
    store_data, cmap_var, times, cmap

    get_data, prefix+'e_eflux', times, eflux, limits=lim
    cmap = get_var_data(cmap_var, at=times)
    var = prefix+'e_eflux_map'
    theta_loss_cones = [45,60]
    ndim = n_elements(theta_loss_cones)
    ntime = n_elements(times)
    eflux_map = fltarr(ntime,ndim+1)
    foreach theta_loss_cone, theta_loss_cones, lc_id do begin
        sr = !dpi*sin(theta_loss_cone*constant('rad'))^2
        eflux_map[*,lc_id+1] = eflux*cmap*sr*1.6e-19*1e4*1e3  ; convert from eV/cm^2-s-sr to mW/m^2.
    endforeach
    store_data, var, times, eflux_map
    yrange = [2e-1,9e2]
    constant = 10d^[0,1,2]
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', yrange, $
        'constant', constant, $
        'display_type', 'stack', $
        'labels', ['Aurora',tex2str('theta')+'!DLC!N='+string(theta_loss_cones,format='(I0)')+'!U'+tex2str('circ')], $
        'colors', sgcolor(['red','green','blue']), $
        'ytitle', '(mW/m!U2!N)' )
    
    
    var = prefix+'e_en_spec'
    get_data, var, times, data, vals
    index = where(data eq 0 or finite(data,nan=1), count)
    if count ne 0 then begin
        data[index] = 0.01
        store_data, var, times, data, vals
    endif
    
    ytickv = 10d^[2,3,4]
    ytickn = '10!U'+string([2,3,4],format='(I0)')
    zrange = [1e5,5e8]
    ztickv = 10d^[5,6,7,8]
    zticks = n_elements(ztickv)-1
    zminor = 9
    ztickn = '10!U'+string([5,6,7,8],format='(I0)')
    ztickn[0:*:2] = ' '
    zticklen = uniform_ticklen/xchsz*1/fig_size[0]
    options, var, 'constant', ytickv
    options, var, 'ytickv', ytickv
    options, var, 'yticks', n_elements(ytickv)-1
    options, var, 'yminor', 10
    options, var, 'ytickname', ytickn
    options, var, 'ytitle', 'Energy!C(eV)'
    options, var, 'zcharsize', label_size
    options, var, 'zrange', zrange
    options, var, 'ztickname', ztickn
    options, var, 'ztickv', ztickv
    options, var, 'zticks', zticks
    options, var, 'zminor', zminor
    options, var, 'zticklen', zticklen
    options, var, 'color_table', 65



;---electron.
    e_eflux_var = prefix+'e_eflux_map'
    e_efluxs = get_var_data(e_eflux_var, in=dmsp_arc_time_range, times=e_times)
    avg_e_efluxs = total(e_efluxs[*,1:2],1)/n_elements(e_efluxs[*,0])
    foreach theta_loss_cone, theta_loss_cones, lc_id do begin
        msg = string(theta_loss_cone,format='(I0)')+' deg, '+string(avg_e_efluxs[lc_id],format='(F5.1)')+' mW/m!U2!N'
        print, msg
    endforeach
    
;---arc.
    interp_time, ssusi_eflux_var, to=e_eflux_var
    get_data, ssusi_eflux_var, arc_times, arc_efluxs
    index = where_pro(arc_times-1, '[]', dmsp_arc_time_range)
    avg_arc_eflux = mean(arc_efluxs[index])
    msg = 'Arc eflux: '+string(avg_arc_eflux,format='(F5.1)')
    print, msg
    

;---O+ outflow.
    
    
    
    
;---Pflux.
    rbspa_pflux_var = 'rbspa_pfdot0_fac_spinfit_phasef_map'
    rbspa_eflux_time_range = time_double(['2015-04-16/08:07:00','2015-04-16/08:08:30'])
    avg_rbspa_pflux_eflux = mean((get_var_data(rbspa_pflux_var, in=rbspa_eflux_time_range))[*,0])
    msg = 'RBSP-A pflux: '+string(avg_rbspa_pflux_eflux,format='(F5.1)')
    print, msg
    
    rbspb_pflux_var = 'rbspb_pfdot0_fac_survey_map'
    rbspb_eflux_time_range = time_double(['2015-04-16/08:10:00','2015-04-16/08:11:30'])
    rbspb_eflux_time_range = time_double(['2015-04-16/08:10:04','2015-04-16/08:11:00'])
    avg_rbspb_pflux_eflux = mean((get_var_data(rbspb_pflux_var, in=rbspb_eflux_time_range))[*,0])
    msg = 'RBSP-B pflux: '+string(avg_rbspb_pflux_eflux,format='(F5.1)')
    print, msg
    ;       |  S   | e-        | arc   | O+
    ; avg   | 57.3 | 55.2      | 62.5  |
    ; peak  | 350  | 198,132   | 87    |
    
    
    eflux_vars = [$
        'dmspf19_'+['e_eflux_map','ssusi_eflux'], $
        'rbspa_'+['o_eflux_map','pfdot0_fac_spinfit_phasef_map'],'rbspb_pfdot0_fac_survey_map']
    foreach eflux_var, eflux_vars do begin
        eflux_map = (get_var_data(eflux_var, times=times))[*,0]
        if eflux_var eq 'dmspf19_e_eflux_map' then begin
            eflux_map = (get_var_data(eflux_var, times=times))[*,-1]
            stop
        endif
        index = where(finite(eflux_map,nan=1),count)
        if count ne 0 then eflux_map[index] = 0
        eflux_map_mid = (eflux_map[1:-1]+eflux_map[0:-2])*0.5
        dtimes = times[1:-1]-times[0:-2]

        eflux_int = [0,eflux_map_mid*dtimes]
        for ii=0,n_elements(times)-2 do eflux_int[ii+1] = eflux_int[ii]+eflux_map_mid[ii]*dtimes[ii]
        store_data, eflux_var+'_int', times, abs(eflux_int)
    endforeach

    stop



end


print, fig_2015_0416_0800_eflux_v01(event_info=event_info)
end