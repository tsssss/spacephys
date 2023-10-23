;+
; Eflux at themis, dmsp, and aurora.
;-

function fig_2017_0309_0700_eflux_v01, plot_file, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2017_0309_0700'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
test = 0

    themis_info = event_info.themis.thd
    prefix = themis_info['prefix']
    probe = themis_info['probe']
    sc_name = themis_info['sc_name']
    sc_color = themis_info['sc_color']

    time_range = time_double(['2017-03-09/07:20','2017-03-09/07:35'])
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_eflux_'+sc_name+probe+'_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0


;---Load dmsp data.
    dmsp_info = event_info.dmsp.dmspf18
    dmsp_color = dmsp_info.sc_color

    probe = dmsp_info['probe']
    prefix = dmsp_info['prefix']
    ssusi_id = 'energy'
    dmsp_mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id=ssusi_id)
;    mlt_image_var = dmsp_read_mlt_image(time_range, probe=probe, id='1216')
    mlat_vars = dmsp_read_mlat_vars(time_range, probe=probe, errmsg=errmsg)
    ele_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='e', errmsg=errmsg)
    ion_spec_var = dmsp_read_en_spec(time_range, probe=probe, species='p', errmsg=errmsg)
    ele_eflux_var = dmsp_read_eflux(time_range, probe=probe, species='e', errmsg=errmsg)


    ; Calculate the ssusi auroral eflux at dmsp footpoint.
    time_id = 1
    get_data, dmsp_mlt_image_var, times, mlt_images, limits=lim
    mlt_image = reform(mlt_images[time_id,*,*])
    pixel_mlat = lim.pixel_mlat
    pixel_mlt = lim.pixel_mlt

    min_mlat = 50
    pixel_tt = (pixel_mlt*15-90)*constant('rad')
    pixel_rr = (90-pixel_mlat)/(90-min_mlat)
    pixel_xx = pixel_rr*cos(pixel_tt)
    pixel_yy = pixel_rr*sin(pixel_tt)

    dmsp_orbit_time_range = time_double('2017-03-09/'+['07:30','07:35'])

    mlat_var = mlat_vars[0]
    mlt_var = mlat_vars[1]
    mlon_var = mlat_vars[2]
    mlats = get_var_data(mlat_var, in=dmsp_orbit_time_range, times=the_times)
    mlons = get_var_data(mlon_var, at=the_times)
    ; use aacgm mlon.
    mlts = aacgm_mlon2mlt(mlons, the_times)
    sc_tt = (mlts*15-90)*constant('rad')
    sc_rr = (90-abs(mlats))/(90-min_mlat)
    sc_xx = sc_rr*cos(sc_tt)
    sc_yy = sc_rr*sin(sc_tt)
    ntime = n_elements(the_times)
    sussi_eflux = fltarr(ntime)
    mlt_image = mlt_image[*]
    for ii=0,ntime-1 do begin
        dis = sqrt((pixel_xx-sc_xx[ii])^2+(pixel_yy-sc_yy[ii])^2)
        tmp = min(dis[*], abs=1, index)
        sussi_eflux[ii] = mlt_image[index]
    endfor
    ssusi_eflux_var = prefix+'ssusi_eflux'
    store_data, ssusi_eflux_var, the_times, sussi_eflux
    dmsp_eflux_yrange = [1e-1,1e2]
    add_setting, ssusi_eflux_var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', dmsp_eflux_yrange, $
        'display_type', 'scalar', $
        'short_name', 'Aurora eflux', $
        'unit', 'mW/m!U2!N' )
    
    
    ; map eflux and convert unit.
    get_data, prefix+'e_eflux', times, eflux, limits=lim
    var = prefix+'e_eflux_map'
    cmap = 1.4  ; this is for 800 km.
    theta_loss_cone = 45d
    sr = !dpi*sin(theta_loss_cone*constant('rad'))^2
    eflux_map = eflux*cmap*sr*1.6e-19*1e4*1e3  ; convert from eV/cm^2-s-sr to mW/m^2.
    store_data, var, times, eflux_map
    add_setting, var, smart=1, dictionary($
        'ylog', 1, $
        'yrange', dmsp_eflux_yrange, $
        'display_type', 'scalar', $
        'short_name', tex2str('Gamma')+'!De!N', $
        'unit', 'mW/m!U2!N' )

    dmsp_eflux_var = prefix+'eflux_combo'
    ntime = n_elements(times)
    loss_cones = [30d,50]
    nloss_cone = n_elements(loss_cones)
    labels = tex2str('theta')+'!DLC!N = '+string(loss_cones,format='(I0)')+' deg'
    efluxs = fltarr(ntime,nloss_cone)
    foreach loss_cone,loss_cones, lc_id do begin
        sr = !dpi*sin(loss_cone*constant('rad'))^2
        eflux_map = eflux*cmap*sr*1.6e-19*1e4*1e3  ; convert from eV/cm^2-s-sr to mW/m^2.
        efluxs[*,lc_id] = eflux_map
    endforeach
    interp_time, ssusi_eflux_var, times
    ssusi_eflux = get_var_data(ssusi_eflux_var)
    efluxs = [[ssusi_eflux],[efluxs]]
    index = where(times lt time_double('2017-03-09/07:30'))
    efluxs[index,*] = !values.f_nan
    labels = ['Aurora',labels]
    colors = sgcolor(['red','green','blue'])
    store_data, dmsp_eflux_var, times, efluxs
    add_setting, dmsp_eflux_var, smart=1, dictionary($
        'display_type', 'stack', $
        'ylog', 1, $
        'yrange', dmsp_eflux_yrange, $
        'ytitle', '(mW/m!U2!N)', $
        'labels', labels, $
        'colors', colors )

    themis_info = event_info.themis.thd

    prefix = themis_info['prefix']
    probe = themis_info['probe']
    sc_name = themis_info['sc_name']
    sc_color = themis_info['sc_color']


    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_eflux_'+sc_name+probe+'_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0


    ; Settings.
    fac_labels = ['||',tex2str('perp')+','+['west','out']]
    model_settings = themis_info['model_setting']
    external_model = 't89'
    internal_model = 'dipole'
    models = model_settings['models']

    ; the times when E_dsl is more reliable.
    ; Calculate Edot0 angle.
    b0_dsl_var = prefix+'b0_themis_dsl'
    b0_dsl = get_var_data(b0_dsl_var, times=times)
    edot0_angle_var = prefix+'edot0_angle'
    edot0_angle = asin(b0_dsl[*,2]/snorm(b0_dsl))*constant('deg')
    store_data, edot0_angle_var, times, edot0_angle
    add_setting, edot0_angle_var, smart=1, dictionary($
        'display_type', 'scalar', $
        'short_name', 'edot0 angle', $
        'unit', 'deg', $
        'constant', [20,25] )
    edot0_angle = get_var_data(prefix+'edot0_angle', times=times, in=time_range)
    index = where(edot0_angle ge 25, count)
    trs = times[time_to_range(index,time_step=1)]
    durs = trs[*,1]-trs[*,0]
    index = where(durs ge 180, ntr)
    snapshot_trs = trs[index,*]
    nsnapshot = n_elements(snapshot_trs[*,0])



    keo_var = 'thg_asf_keo'
    mlt_image_rect_var = themis_asf_read_mlt_image_rect(time_range, get_name=1)
    fmlt_var = prefix+'fmlt_'+internal_model+'_north'
    dmlt = 0.1
    mlat_range = [65,71]

    get_data, mlt_image_rect_var, times, data, limits=lim
    ntime = n_elements(times)
    model_index = where(models eq external_model)
    fmlts = (get_var_data(fmlt_var, at=times))[*,model_index]

    mlt_bins = lim.mlt_bins
    mlat_bins = lim.mlat_bins
    mlat_index = where_pro(mlat_bins, '[]', mlat_range, count=nybin)
    if nybin eq 0 then begin
        errmsg = 'Invalid MLat range ...'
        return, retval
    endif
    ybins = mlat_bins[mlat_index]

    hehe = fltarr(ntime,nybin)
    foreach time, times, time_id do begin
        the_mlt_range = fmlts[time_id]+[-1,1]*dmlt
        mlt_index = where_pro(mlt_bins, '[]', the_mlt_range, count=nxbin)
        hehe[time_id,*] = total(data[time_id,mlt_index,mlat_index],2)/nxbin
    endforeach
    yrange = mlat_range
    ytitle = 'MLat (deg)'
    
    if n_elements(ystep) eq 0 then ystep = 5
    ytickv = make_bins(yrange, ystep, inner=1)
    yticks = n_elements(ytickv)-1
    yminor = ystep
    store_data, keo_var, times, hehe, ybins
    add_setting, keo_var, smart=1, dictionary($
        'display_type', 'spec', $
        'short_name', 'ASI Count', $
        'unit', '#', $
        'ytitle', 'MLat!C(deg)', $
        'zrange', [5e2,2e4], $
        'ztickv', [1e3,1e4], $
        'zticks', 1, $
        'ztickname', '10!U'+['3','4'], $
        'zminor', 9, $
        'zlog', 1, $
        'zticklen', zticklen, $
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor, $
        'color_table', 49 )


    ; themis pflux at 100 km.
    pflux_setting = themis_info['pflux_setting']
    scale_info = pflux_setting['scale_info']
    pf_fac_var = prefix+'pfdot0_fac_map'
    get_data, pf_fac_var, times, pf_fac
    pf_combo = [[-pf_fac[*,0]],[snorm(pf_fac)]]
    ntime = n_elements(times)
    flags = fltarr(ntime)
    for ii=0,nsnapshot-1 do begin
        index = where_pro(times, '[]', snapshot_trs[ii,*], count=count)
        if count ne 0 then flags[index] = 1
    endfor
    index = where(flags eq 0, count)
    if count ne 0 then pf_combo[index,*] = !values.f_nan
    th_pf_var = prefix+'pf_combo'
    store_data, th_pf_var, times, pf_combo
    add_setting, th_pf_var, smart=1, dictionary($
        'display_type', 'stack', $
        'ytitle', '(mW/m!U2!N)', $
        'labels', ['Earthward S', '|S|'], $
        'ylog', 1, $
        'yrange', [0.1,300], $
        'colors', sgcolor(['red','black']) )
    
    b_var = prefix+'z_gsm'
    ntime = n_elements(times)
    ndim = 3
    z_gsm = fltarr(ntime,ndim)
    z_gsm[*,2] = 1
    store_data, b_var, times, z_gsm
    add_setting, b_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', 'Z', $
        'coord', 'gsm', $
        'unit', '#' )
    b_var = prefix+'b_gsm'
    define_fac, b_var, prefix+'r_gsm', time_var=b_var
    vars = prefix+'p_'+['vbulk','eflux','enthalpy','keflux','hflux']
    foreach var, vars do begin
        to_fac, var+'_gsm', to=var+'_fac'
    endforeach
    
    
    ; plot.
    fig_size = [6,4]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    margins = [12,4,10,1]
    poss = sgcalcpos(3, margins=margins)
    vars = [th_pf_var,dmsp_eflux_var]
    options, vars, 'yrange', [0.2,300]
    options, vars, 'constant', [1,10,100]
    options, vars, 'ylog', 1
    
    
    plot_vars = [keo_var,th_pf_var,dmsp_eflux_var]
    nplot_var = n_elements(plot_vars)
    labels = letters(nplot_var)+') '+['Aurora!CN-hem','S abv AAR','e- eflux!Cblw AAR']
    tplot, plot_vars, trange=time_range, position=poss
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*10
        ty = tpos[3]-ychsz*0.6
        msg = labels[pid]
        xyouts, tx,ty,normal=1, msg
    endfor

    ; Add FMLat.
    fmlat_var = prefix+'fmlat_'+internal_model+'_north'
    the_fmlat_var = prefix+'fmlat_north_smooth'
    dfmlat = [0d,1]
    model_index = (where(models eq external_model))[0]
    get_data, fmlat_var, uts, fmlats
    fmlats = fmlats[*,model_index]
    store_data, the_fmlat_var, uts, fmlats
    pid = where(plot_vars eq keo_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]

        fmlats = get_var_data(the_fmlat_var, times=times)
        get_data, keo_var, limits=lim
        xrange = time_range
        yrange = lim.yrange
        plot, xrange, yrange, $
            xstyle=5, xrange=xrange, $
            ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
;            alpha = 0.5
;            white = sgcolor('white')
;            the_color = alpha*sc_color+(1-alpha)*white
        the_color = sc_color
        oplot, times, fmlats+dfmlat[0], color=the_color, linestyle=2
        oplot, times, fmlats+dfmlat[1], color=the_color, linestyle=2
        tx = time_range[0]
        ty = interpol(fmlats,times,tx)+mean(dfmlat)
        tmp = convert_coord(tx,ty, data=1, to_normal=1)
        tx = tmp[0]+xchsz*0.5
        ty = tmp[1]-ychsz*0.3
        xyouts, tx,ty,normal=1, strupcase(sc_name+'-'+probe), charsize=label_size, alignment=0, color=sc_color
    endif
    
    ; Add ssusi
    pid = where(plot_vars eq dmsp_eflux_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        
        get_data, dmsp_eflux_var, times, efluxs, limits=lim
        xrange = time_range
        yrange = lim.yrange
        xlog = 0
        ylog = lim.ylog
        plot, xrange, yrange, $
            xlog=xlog, xstyle=5, xrange=xrange, $
            ylog=ylog, ystyle=5, yrange=yrange, $
            position=tpos, nodata=1, noerase=1
        
        ssusi_eflux = efluxs[*,0]
        the_color = lim.colors[0]
        index = where(ssusi_eflux ne 0, count)

        oplot, times[index], ssusi_eflux[index], psym=1, symsize=0.25, color=the_color
    endif
    
    
    if keyword_set(test) then stop
    sgclose
    return, plot_file

end

print, fig_2017_0309_0700_eflux_v01(event_info=event_info)
end
