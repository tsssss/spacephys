; Read PAD from cdaweb.

function tracers_read_pad_ele, input_time_range, probe=probe, $
    update=update, get_name=get_name, errmsg=errmsg, _extra=extra
    compile_opt idl2

    errmsg = ''
    retval = !null
    prefix = 'ts'+probe+'_'
    suffix = '_cdaweb'

    var_info = prefix+'pad_ele'+suffix
    if keyword_set(get_name) then return, var_info
    if keyword_set(update) then del_data, var_info
    time_range = time_double(input_time_range)
    if ~check_if_update(var_info, time_range) then return, var_info

    files = tracers_load_ace(time_range, probe=probe, errmsg=errmsg, id='cdaweb%l2')
    if errmsg ne '' then return, retval

    ; Read B field in ts_mag.
    b_var = tracers_read_bfield(time_range, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval

    ; Read PAD variables.
    var_list = list()
    time_var = 'Epoch'
    anode_var = prefix+'l2_ace_TSCS_anode_angle'
    anode_bins = cdf_read_var(anode_var, filename=files[0])
    en_bin_var = prefix+'l2_ace_energy'
    energy_bins = cdf_read_var(en_bin_var, filename=files[0])
    nenergy_bin = n_elements(energy_bins)
    de = abs(energy_bins[1:nenergy_bin-1]-energy_bins[0:nenergy_bin-2])
    de_e = mean(de/energy_bins[0:nenergy_bin-2])



    flux_var = prefix+'l2_ace_def'
    var_list.add, dictionary($
        'in_vars', flux_var, $
        'time_var_name', time_var, $
        'time_var_type', 'tt2000' )
    read_vars, time_range, files=files, var_list=var_list, errmsg=errmsg
    if errmsg ne '' then return, retval

    fluxs = var_get_data(flux_var, times=times)
    fluxs *= de_e   ; Convert from 1/cm^2-s-sr to eV/cm^2-s-sr-eV.
    unit = 'eV/cm!U2!N-s-sr-eV'
    b_ts_mag = var_get_data(b_var, at=times, settings=settings)
    b_tscs = cotran_pro(b_ts_mag, probe=probe, coord_msg='ts_'+['mag','tscs'])
    settings['coord'] = 'ts_tscs'
    b_var = prefix+'b_tscs'
    b_var = var_store(b_var, b_tscs, times, settings=settings)

    tmp_var = prefix+'pad_tmp'
    store_data, tmp_var, times, fluxs[*,40,0]
    options, tmp_var, ylog=1, yrange=[1e3,1e10]

    ; Calc pitch angles.
    ntime = n_elements(times)
    nanode = n_elements(anode_bins)
    ndim = 3
    rad = constant('rad')
    anode_dirs = fltarr(nanode,ndim)
    anode_dirs[*,0] = 0
    anode_dirs[*,1] = sin(anode_bins*rad)
    anode_dirs[*,2] = cos(anode_bins*rad)
;   PA sign is physically wrong if do this.
;    anode_dirs = -anode_dirs ; The flux direction is opposite to the anode direction.
    pa_angles = fltarr(ntime, nanode)
    for ii=0,ntime-1 do begin
        for jj=0,nanode-1 do begin
            pa_angles[ii,jj] = sang(anode_dirs[jj,*], b_tscs[ii,*], degree=1)
        endfor
    endfor

    spin_period = 3d
    sa_times = smkarthm(times[0],times[ntime-1],spin_period,'dx')
    sa_times = sa_times+spin_period*0.5
    nsa_time = n_elements(sa_times)
    pa_range = [0d,180]
    pa_bin_size = 15d
    pa_bins = smkarthm(pa_range[0], pa_range[1], pa_bin_size, 'dx')
    pa_bin_centers = pa_bins[0:-2]+pa_bin_size*0.5
    npa_bin = n_elements(pa_bin_centers)
    pad_fluxs = fltarr(nsa_time, npa_bin, nenergy_bin)
    for ii=0,nsa_time-1 do begin
        sa_index = where_pro(times, '[]', sa_times[ii]+[-1,1]*spin_period*0.5, count=count)
        if count eq 0 then continue
        npa = count*nanode
        the_fluxs = fltarr(npa, nenergy_bin)
        for jj=0,nenergy_bin-1 do begin
            the_fluxs[*,jj] = fluxs[sa_index,jj,*]
        endfor
        the_pas = reform(pa_angles[sa_index,*], count*nanode)
        for jj=0,npa_bin-1 do begin
            in_pa = where_pro(the_pas, '[]', pa_bins[jj]+[-1,1]*pa_bin_size*0.5, count=count)
            if count eq 0 then continue
            the_fluxs_in_bin = the_fluxs[in_pa,*]
            pad_fluxs[ii,jj,*] = total(the_fluxs_in_bin, 1)/count
        endfor
    endfor

    energy_range = [1e3, 1e5]
    ;energy_range = [50,1e3]
    energy_index = where_pro(energy_bins, '[]', energy_range, count=count)
    pa_specs = total(pad_fluxs[*,*,energy_index], 3)
    pa_spec_var = prefix+'ele_pa_spec'+suffix
    var_info = var_store(pa_spec_var, pa_specs, sa_times, pa_bin_centers)
    ytickv = make_bins(pa_range, 45)
    yticks = n_elements(ytickv)-1
    yminor = 5
    ztitle = '('+unit+')'
    add_setting, var_info, smart=1, dictionary($
        'display_type', 'spec', $
        'ylog', 0, $
        'yminor', yminor, $
        'yrange', pa_range, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'ytitle', 'PA!C(deg)', $
        'ztitle', ztitle, $
        'zlog', 1)
    
    en_spec_var = prefix+'ele_en_spec'+suffix
    var_info = var_store(en_spec_var, total(pad_fluxs,2,nan=1), sa_times, energy_bins)
    add_setting, var_info, smart=1, dictionary($
        'display_type', 'spec', $
        'ylog', 1, $
        'ytitle', 'Energy!C(eV)', $
        'ztitle', ztitle, $
        'zlog', 1)

    plot_vars = [pa_spec_var, en_spec_var]
    options, plot_vars, zrange=[1e4,1e10]
;    tplot, plot_vars
;    stop

    return, var_info

end


compile_opt idl2
time_range = ['2026-02-16/04:00','2026-02-16/04:05']
probe = '2'
pad_var = tracers_read_pad_ele(time_range, probe=probe)
print, pad_var
end