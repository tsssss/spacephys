
time_range = ['2016-03-05/12:50','2016-03-06/02:40']
date = '2016-03-07'
mission_probe = ['mms','4']

secofday = constant('secofday')
tr_list = list()
tr_list.add, time_double(date)+[-1,0]*secofday
tr_list.add, time_double(date)+[0,1]*secofday

mission = mission_probe[0]
probe = mission_probe[1]
prefix = mission+probe+'_'
project_info = micro_injection_stat_load_project()
sample_energys = project_info['sample_energys']
scale_info = project_info['scale_info']
routine = mission+'_read_kev_electron_cdaweb'

plot_var_list = list()
foreach tr, tr_list, tid do begin
    time_range = time_double(tr)
    data_time_range = time_range+[-1,1]*scale_info['s1']*4
    en_spec_var = call_function(routine, data_time_range, probe=probe)

    var_list = list()
    var_list.add, en_spec_var

    foreach sample_energy, sample_energys do begin
        energy_str = string(sample_energy, format='(I0)')

        ; Obtain the flux for the given energy.
        fluxs = get_var_data(en_spec_var, times=times, energys, settings=settings)
        tmp = min(energys-sample_energy, abs=1, energy_index)
        if energy_index eq 0 or energy_index eq n_elements(energys) then begin
            errmsg = 'Energy '+energy_str+'keV not found in the data.'
            continue
        endif
        the_fluxs = fluxs[*,energy_index]
        flux_var = prefix+'kev_ele_'+energy_str+'kev'
        flux_var = var_store(flux_var, the_fluxs, times)
        unit = settings['unit']
        options, flux_var, ylog=1, ytitle='e- flux!C('+unit+')', labels=energy_str+'keV', yrange=[1e-2,1e6]

        ; Interpolate to common_times.
        time_step = 20d
        common_times = make_bins(data_time_range, time_step)
        foreach var, [flux_var] do begin
            flux = get_var_data(var, at=common_times)
            var = var_store(var, flux, common_times)
        endforeach

        ; Calculate cwt.
        spec_var = stplot_mor_new(flux_var, scale_info=scale_info)
        options, spec_var, zlog=1, zrange=[1e-1,1e8]
        
        ; Log version.
        log_flux_var = prefix+'log_kev_ele_'+energy_str
        
        the_fluxs = get_var_data(flux_var)
        log_fluxs = alog10(the_fluxs)
        index = where(finite(log_fluxs))
        log_fluxs = interpol(log_fluxs[index], common_times[index], common_times)
        log_flux_var = var_store(log_flux_var, log_fluxs, common_times)
        options, log_flux_var, ylog=0, ytitle='Log e- flux!C('+unit+')', labels=energy_str
        
        log_spec_var = stplot_mor_new(log_flux_var, scale_info=scale_info)
        options, log_spec_var, zlog=1, zrange=[1e-4,1e1]

        var_list.add, [flux_var,spec_var,log_flux_var,log_spec_var], extract=1
    endforeach

    vars = var_list.toarray()
    suffix = '_day'+string(tid+1,format='(I0)')
    new_vars = vars+suffix
    foreach var, vars, vid do begin
        var = rename_var(var, output=new_vars[vid])
    endforeach
    plot_var_list.add, new_vars
endforeach


nxpan = n_elements(plot_var_list)
nypan = n_elements(plot_var_list[0])

plot_file = 0
sgopen, plot_file, size=[12,8]
margins = [10,4,12,1]
poss = sgcalcpos(nypan, nxpan, margins=margins, xpad=0)
tplot_options, 'tickinterval', 3600*4

foreach plot_vars, plot_var_list, pid do begin
    tr = tr_list[pid]
    if pid eq 1 then begin
        options, plot_vars, ytitle=' ', ytickformat='(A1)', yminor=0, yticks=1, yticklen=1e-5
    endif else begin
        options, plot_vars, ystyle=9, labels=' '
    endelse
    tplot, plot_vars, trange=tr, position=reform(poss[*,pid,*]), noerase=1, novtitle=1
endforeach


end