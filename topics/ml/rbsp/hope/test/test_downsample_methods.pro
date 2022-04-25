;+
; Test downsampling 
;-

time_range = time_double(['2013-06-07','2013-06-08'])
time_range = time_double(['2015-03-15','2015-03-18'])
probe = 'a'
fillval = 0.01
downsample_methods = ['mean','median','log_mean']
test_time_steps = [1,2,5]*60d
all_species = rbsp_hope_species()
all_species = ['e','p','o','he']
all_species = ['e']

test = 0

foreach time_step, test_time_steps do begin
    common_times = make_bins(time_range+[0,-1]*time_step,time_step)+time_step*0.5
    ntime = n_elements(common_times)
    prefix = 'rbsp'+probe+'_'

    files = rbsp_load_hope(time_range, probe=probe, id='l3%pa')
    foreach species, all_species do begin
        energy_var = 'energy_bins_'+species
        energy_bins = ml_rbsp_hope_read_energy_bins(probe=probe, species=species)

        nenergy = n_elements(energy_bins)
        energy_range = minmax(energy_bins)
        energy_step = mean(energy_bins[1:nenergy-1]/energy_bins[0:nenergy-2])
        energy_edges = sqrt(energy_bins[1:nenergy-1]*energy_bins[0:nenergy-2])
        energy_edges = [energy_edges[0]/energy_step,energy_edges,energy_edges[-1]*energy_step]
        en_spec_var = prefix+'hope_flux_'+species
        flux_var = strupcase('f'+species+'do')

        var_list = list()
        the_time_var = (species eq 'e')? 'Epoch_Ele': 'Epoch_Ion'
        energy_var = (species eq 'e')? 'HOPE_ENERGY_Ele': 'HOPE_ENERGY_Ion'
        var_list.add, dictionary($
            'in_vars', flux_var, $
            'out_vars', en_spec_var, $
            'time_var_name', the_time_var, $
            'time_var_type', 'epoch' )
        var_list.add, dictionary($
            'in_vars', energy_var, $
            'time_var_name', the_time_var, $
            'time_var_type', 'epoch' )
        read_vars, time_range, files=files, var_list=var_list
        get_data, en_spec_var, times, flux, limits=lim, dlimits=dlim
        get_data, energy_var, times, energys
        index = where(flux le 0 or flux ge 1e30, count)
        if count ne 0 then flux[index] = fillval

        ; Calculate the lowres data.
        plot_vars = list(en_spec_var)
        foreach downsample_method, downsample_methods do begin
            flux_lowres = fltarr(ntime,nenergy)+fillval
            foreach time, common_times, time_id do begin
                time_index = where(times ge time-0.5*time_step and $
                    times le time+0.5*time_step, count)
                if count eq 0 then continue
                the_energy = energys[time_index,*]
                the_flux = flux[time_index,*]
                foreach energy, energy_bins, energy_id do begin
                    index = where(the_energy ge energy_edges[energy_id] and $
                        the_energy le energy_edges[energy_id+1], count)
                    if count eq 0 then continue
                    if downsample_method eq 'mean' then begin
                        flux_lowres[time_id,energy_id] = mean(the_flux[index], nan=1)
                    endif else if downsample_method eq 'median' then begin
                        tmp = the_flux[index]
                        flux_lowres[time_id,energy_id] = median(the_flux[index])
                        ;the_index = where(tmp ne fillval, count)
                        ;if count ne 0 then flux_lowres[time_id,energy_id] = median(tmp[the_index])
                    endif else if downsample_method eq 'log_mean' then begin
                        flux_lowres[time_id,energy_id] = 10^mean(alog10(the_flux[index]), nan=1)
                    endif
                endforeach
            endforeach
            lowres_en_spec_var = en_spec_var+'_lowres_'+downsample_method
            store_data, lowres_en_spec_var, common_times, flux_lowres, energy_bins, limits=lim, dlimits=dlim
            plot_vars.add, lowres_en_spec_var
        endforeach
                    
        store_data, en_spec_var, times, flux, energys
        vars = plot_vars.toarray()
        zlim, vars, 1e0, 1e7, 1
        if species eq 'e' then zlim, vars, 1e3, 1e10, 1
        ylim, vars, 4, 4e4, 1
        options, vars, 'spec', 1
        options, vars, 'no_interp', 1
        options, vars, 'xticklen', -0.02
        options, vars, 'yticklen', -0.01
        options, vars, 'ytitle', '(eV)'
        options, vars, 'ztitle', '(#/cm!U2!N-s-sr-keV)'
        options, vars, 'color_table', 61
        
        base = 'fig_test_hope_downsample'+$
            '_'+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_to_')+$
            '_'+'rbsp'+probe+$
            '_'+species+$
            '_'+string(time_step/60,format='(I0)')+'min'+$
            '_v01.pdf'
        fig_file = join_path([srootdir(),'fig',base])
        if keyword_set(test) then fig_file = 0
        sgopen, fig_file, xsize=8, ysize=6, xchsz=xchsz, ychsz=ychsz
        nvar = n_elements(vars)
        margins = [12,4,10,2]
        poss = sgcalcpos(nvar,margins=margins)
        tplot, vars, trange=time_range, position=poss
        fig_labs = letters(nvar)+'. '+['orig',downsample_methods]
        for ii=0,nvar-1 do begin
            tpos = poss[*,ii]
            tx = xchsz*2
            ty = tpos[3]-ychsz*0.7
            xyouts, tx,ty,normal=1, fig_labs[ii]
        endfor
        
        tpos = poss[*,0]
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,normal=1, base
        
        if keyword_set(test) then stop
        sgclose
    endforeach
endforeach

end