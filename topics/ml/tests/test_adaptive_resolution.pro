;+
; Test to generate CDFs at an adaptive temporal resolution.
;-

test_list = list()
test_list.add, dictionary($
    'time_range', ['2015-03-15','2015-03-18'], $
    'probes', 'a', $
    'species', ['o','p'], $
    'energys', [300d] )
    test_list.add, dictionary($
        'time_range', ['2013-06-07','2013-06-08'], $
        'probes', 'a', $
        'species', ['o'], $
        'energys', [300d] )
test = 0

resolutions = [60,300]
foreach info, test_list do begin
    time_range = time_double(info.time_range)
    foreach probe, info.probes do begin
        energy_bins = ml_rbsp_hope_read_energy_bins(species=species, probe=probe)
        prefix = 'rbsp'+probe+'_'

        foreach species, info.species do begin
            foreach energy, info.energys do begin
                tmp = min(energy_bins-energy, energy_index, abs=1)
                energy = energy_bins[energy_index]
                energy_str = string(round(energy),format='(I0)')
                
                species_str = rbsp_hope_species_name(species)
                dis_var = ml_rbsp_read_dis(time_range, probe=probe)
                
                
                spec_vars = []
                foreach resolution, resolutions do begin
                    suffix = '_'+string(resolution/60,format='(I0)')+'min'
                    flux_var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, $
                        species=species, resolution=resolution, errmsg=errmsg)
                    rename_var, flux_var, to=flux_var+suffix
                    spec_vars = [spec_vars,flux_var+suffix]
                endforeach
                rbsp_read_en_spec, time_range, probe=probe
                orig_var = prefix+species+'_en_spec'
                spec_vars = [orig_var,spec_vars]
                zlim, spec_vars, 1e0, 1e7, 1
                ylim, spec_vars, 4, 4e4, 1
                options, spec_vars, 'ytitle', 'Energy (eV)'
                get_data, spec_vars[-1], limits=lim
                options, spec_vars, 'ztitle', lim.ztitle
                ztickv = [1e0,1e3,1e6]
                zticks = n_elements(ztickv)
                options, spec_vars, 'ztickv', ztickv
                options, spec_vars, 'zticks', zticks
                
                options, dis_var, 'ytickv', [1,3,5]
                options, dis_var, 'yticks', 2
                options, dis_var, 'yminor', 2
                options, dis_var, 'yrange', [1,6]
                
                
                get_data, orig_var, times, limits=lim
                ntime = n_elements(times)
                nspec_var = n_elements(spec_vars)
                flux = fltarr(ntime,nspec_var)
                flux_var = prefix+'flux'
                for ii=0,nspec_var-1 do begin
                    data = get_var_data(spec_vars[ii], at=times, energys)
                    if spec_vars[ii] eq orig_var then begin
                        for jj=0,ntime-1 do begin
                            flux[jj,ii] = interpol(data[jj,*],energys[jj,*],energy)
                        endfor
                    endif else begin
                        flux[*,ii] = data[*,energy_index]
                    endelse
                endfor
                index = where(flux le 1e-2, count)
                ztitle = lim.ztitle
                ytitle = ztitle+'!C'+energy_str+' eV'
                if count ne 0 then flux[index] = !values.f_nan
                labels = ['orig','1min','5min']
                colors = sgcolor(['black','silver','green'])
                store_data, flux_var, times, flux, limits={$
                    ylog:1, ytitle:ytitle, labels:labels, colors:colors}
                
                
                plot_file = join_path([googledir(),'works','ml_hope','plots','test_adaptive_resolution',$
                    prefix+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_to_')+'_'+species+'_'+energy_str+'eV_test_adaptive_resolution_v01.pdf'])
                if keyword_set(test) then plot_file = test
                sgopen, plot_file, xsize=7, ysize=8, xchsz=xchsz, ychsz=ychsz
                vars = [spec_vars,dis_var,flux_var]
                nvar = n_elements(vars)
                margins = [12,4,10,1]
                poss = sgcalcpos(nvar, ypans=[1,1,1,0.5,1.5], margins=margins)
                tplot, vars, trange=time_range, position=poss
                labels = letters(nvar)+') '+['Orig','1min','5min','|R|',species_str+' Flux']
                for ii=0,nvar-1 do begin
                    tpos = poss[*,ii]
                    tx = xchsz*2
                    ty = tpos[3]-ychsz*0.8
                    xyouts, tx,ty,normal=1, labels[ii]
                endfor
                if keyword_set(test) then stop
                sgclose
                
            endforeach
        endforeach
    endforeach
endforeach

end