pro rbsp_themis_injection_gen_survey_plot_for_rbsp_injection, filename=txt_file, test=test, plot_dir=plot_dir
    ; This will overwrite existing files and only need to run once.


    test = 1

    pinfo = rbsp_themis_injection_load_project()
    if n_elements(plot_dir) eq 0 then plot_dir = join_path([pinfo['plot_dir'],'rbsp_injection'])

    if n_elements(txt_file) eq 0 then begin
        txt_file = rbsp_themis_injection_identify_rbsp_injection()
    endif
    if file_test(txt_file) eq 0 then begin
        txt_file = rbsp_themis_injection_identify_rbsp_injection(txt_file)
    endif
    log_file = txt_file+'.log'

    lines = read_all_lines(log_file, skip_header=1)
    injections = read_all_lines(txt_file)
    injection_id = 0
    probes = ['a','b']
    time_step = 120
    target_energy = 100.

    foreach line, lines do begin
        infos = strsplit(line,' ',extract=1)
        date = time_double(infos[0],tformat='YYYY_MMDD')
        ninjection = total(float(infos[1:2]))
        info_list = list()
        for ii=0,ninjection-1 do begin
            infos = strsplit(injections[injection_id],' ',extract=1)
            probe = strmid(infos[0],4,1)
            time_range = time_double(infos[1:2])
            info_list.add, list(probe, time_range)
            injection_id += 1
        endfor
        
        time_range = date+[0,constant('secofday')]
        foreach probe, probes do begin
            prefix = 'rbsp'+probe+'_'
            
            ; Load electron flux.
            data_time_range = time_range
            flux_var = rbsp_read_kev_electron(data_time_range, probe=probe, spec=1)
            flux_var_spec = flux_var+'_spec'
            rename_var, flux_var, to=flux_var_spec
            uniform_time, flux_var_spec, time_step
            
            
            ; Prepare the flux and distance.
            fluxs = get_var_data(flux_var_spec, energys, times=times, limits=lim)
            tmp = min(energys-target_energy, energy_index, abs=1)
            flux = fluxs[*,energy_index]
            
            invalid_flux = 1e0
            ntime = n_elements(times)
            pad_rec = 2
            index = where(flux le invalid_flux or finite(flux,nan=1), count)
            if count ne 0 then begin
                ranges = time_to_range(index,time_step=1)
                nrange = n_elements(ranges)*0.5
                for ii=0,nrange-1 do begin
                    i0 = (ranges[ii,0]-pad_rec)>0
                    i1 = (ranges[ii,1]+pad_rec)<(ntime-1)
                    fluxs[i0:i1,*] = !values.f_nan
                endfor
            endif
            store_data, flux_var_spec, times, fluxs, energys
            uniform_time, flux_var_spec, time_step
            add_setting, flux_var_spec, smart=1, dictionary($
                'display_type', 'spec', $
                'zlog', 1, $
                'ylog', 1, $
                'ytitle', 'Energy (keV)', $
                'ztitle', 'e!U-!N flux (#/cm!U2!N-s-sr-keV)', $
                'ytickv', [1e2,1e3], $
                'ytickname', ['10!U2','10!U3'], $
                'yticks', 1, $
                'yminor', 10 )

            ; flux.
            ntime = n_elements(times)
            nenergy = n_elements(energys)
            if energy_index ge 1 and energy_index+2 lt nenergy then begin
                store_data, flux_var, times, fluxs[*,energy_index-1:energy_index+2], energys[energy_index-1:energy_index+2]
            endif else begin
                fluxs = fltarr(ntime,nenergy)
                energys = findgen(nenergy)
                store_data, flux_var, times, fluxs, energys
            endelse
            
            add_setting, flux_var, smart=1, {$
                display_type: 'list', $
                ylog: 1, $
                yrange: [1e1,1e6], $
                ytickv: [1e1,1e2,1e3,1e4,1e5,1e6], $
                ytickname: [' ','10!U2',' ','10!U4',' ','10!U6'], $
                yticks: 5, $
                yminor: 10, $
                color_table: 52, $
                unit: '#/cm!U2!N-s-sr-keV', $
                value_unit: 'keV', $
                short_name: 'e!U-!N flux' }
            options, flux_var, 'colors', sgcolor(['gray','red','green','blue'])

            ; dis.
            r_var = rbsp_read_orbit(time_range, probe=probe)
            get_data, r_var, times, r_vec
            dis_var = prefix+'dis'
            dis = snorm(r_vec)
            store_data, dis_var, times, dis
            add_setting, dis_var, smart=1, dictionary($
                'unit', 'Re', $
                'short_name', '|R|', $
                'display_type', 'scalar', $
                'yticks', 2, $
                'yminor', 2, $
                'ytickv', [2,4,6] )
        endforeach

        base = 'rbsp_themis_injection_survey_plot_for_rbsp_injections_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_v01.pdf'
        year_str = time_string(time_range[0],tformat='YYYY')
        plot_file = join_path([plot_dir,year_str,base])
        if keyword_set(test) then plot_file = 0
        if size(plot_file,type=1) eq 7 then begin
            thick = 20
        endif else begin
            thick = 4
        endelse
        sgopen, plot_file, xsize=8, ysize=8
        plot_vars = []
        ypans = []
        foreach probe, probes do begin
            prefix = 'rbsp'+probe+'_'
            plot_vars = [plot_vars,prefix+['kev_e_flux_spec','kev_e_flux','dis']]
            ypans = [ypans,[1,1,0.5]]
        endforeach
        nvar = n_elements(plot_vars)
        margins = [12,4,10,1]
        poss = sgcalcpos(nvar, margins=margins, ypans=ypans, xchsz=xchsz, ychsz=ychsz)

        tplot, plot_vars, position=poss, noerase=1
        
        foreach probe, probes, probe_id do begin
            tpos = poss[*,probe_id*nvar*0.5+1]
            plot, time_range, [0,1], position=tpos, $
                xstyle=5, xrange=time_range, $
                ystyle=5, yrange=[0,1], $
                noerase=1, nodata=1
                
            foreach info, info_list do begin
                if info[0] ne probe then continue
                plots, info[1], [1,1], thick=thick, color=sgcolor('red')
            endforeach
        endforeach
        
        fig_letters = letters(nvar)+') '
        foreach probe, probes, probe_id do begin
            i0 = probe_id*nvar*0.5
            i1 = i0+nvar*0.5-1
            fig_letters[i0:i1] += ['Spec','Flux','|R|']
        endforeach
        for ii=0,nvar-1 do begin
            tpos = poss[*,ii]
            tx = tpos[0]-xchsz*9
            ty = tpos[3]-ychsz*0.7
            xyouts, tx,ty,normal=1, fig_letters[ii]
            
            msg = (ii lt nvar*0.5)? 'RBSP-A': 'RBSP-B'
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty,normal=1, msg
        endfor

        if keyword_set(test) then stop
        sgclose
    endforeach

    stop



end



rbsp_themis_injection_gen_survey_plot_for_rbsp_injection

end