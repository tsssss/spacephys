;+
; Load HOPE data from Cristian's moments.
; Plot H and O eflux vs Poynting flux.
;-


test = 0

;---Load data.
    id = '2013_0607_event_info'
    if tnames(id) eq '' then load = 1 else load = 0
    if load eq 1 then _2013_0607_load_data3

    hope_data_dir = join_path([srootdir(),'data','hope_moment_cristian'])
    probes = ['a','b']
    species = dictionary($
        'p', 'h', $
        'o', 'o', $
        'he', 'he', $
        'e', 'ele')
    mom_vars = dictionary($
        'V_fac', 'vbulk_fac', $
        'EFlux_fac', 'energy_flux_fac')
    rgb = sgcolor(['red','green','blue'])
    plot_time_range = time_double(['2013-06-07/04:52','2013-06-07/05:02'])

    
;    foreach probe, probes do begin
;        foreach the_species, species.keys() do begin
;            data_file = file_search(join_path([hope_data_dir,'rbsp'+probe]), 'RBSP'+probe+'_*_'+species[the_species]+'.tplot')
;            if file_test(data_file) eq 0 then stop
;            tplot_restore, filename=data_file
;            foreach the_var, mom_vars.keys() do begin
;                new_var = 'rb'+probe+'_'+the_species+'_'+mom_vars[the_var]
;                tplot_rename, the_var, new_var
;                options, new_var, 'colors', rgb
;                get_data, new_var, uts, dat
;                dat[*,2] = -dat[*,2]    ; sign of parallel is wrong, should be anti-parallel.
;                store_data, new_var, uts, dat
;            endforeach
;        endforeach
;    endforeach


;---MHD velocity, use H and O.
    event_info = get_var_data(id)
    foreach probe, probes do begin
        pre0 = 'rb'+probe+'_'
        the_info = (probe eq 'a')? event_info.rbspa: event_info.rbspb
        the_time_range = the_info.time_range
        v_o_fac = get_var_data(pre0+'o_vbulk_fac', in=the_time_range)
        v_o_fac = reform(total(v_o_fac,1)/n_elements(v_o_fac[*,0]))
        v_h_fac = get_var_data(pre0+'p_vbulk_fac', in=the_time_range)
        v_h_fac = reform(total(v_h_fac,1)/n_elements(v_h_fac[*,0]))
        rho_o = the_info.o.n*the_info.o.mass
        rho_h = the_info.h.n*the_info.h.mass
        v_mhd = (v_o_fac*rho_o+v_h_fac*rho_h)/(rho_o+rho_h)
        print, 'V MHD (km/s)' , v_mhd
        print, 'n(O+)/n(O+H+)', the_info.o.n/(the_info.o.n+the_info.h.n)
    endforeach


;---Combine H and O eflux.
    foreach probe, probes do begin
        pre0 = 'rb'+probe+'_'
        
        cmap = get_var_data(pre0+'cmap', in=plot_time_range, limits=lim)
        model_index = where(lim.labels eq default_model)
        cmap0 = mean(cmap[*,model_index])
        
        o_eflux = get_var_data(pre0+'o_energy_flux_fac', times=times)
        p_eflux = get_var_data(pre0+'p_energy_flux_fac', at=times)
        eflux = [[p_eflux[*,2]],[o_eflux[*,2]]]
        eflux = o_eflux[*,2]*cmap0

        the_var = pre0+'ion_eflux'
;        store_data, the_var, times, eflux, limits={$
;            ytitle: 'Eflux!C(mW/m!U2!N)', $
;            colors: sgcolor(['black','red']), $
;            labels: ['H+','O+']+' '+tex2str('Gamma')+'!D||!N'}
        store_data, the_var, times, eflux, limits={$
            ytitle: '@100km!C(mW/m!U2!N)', $
            yrange: [-8,2], $
            ystyle: 1, $
            yticks: 2, $
            yminor: 5, $
            labels: 'O+ '+tex2str('Gamma')+'!D||!C  energy!C  flux'}
            
        pf = get_var_data(pre0+'pf_fac', times=times)
        the_var = pre0+'pf_para'
        yrange = (probe eq 'a')? [-20,60]: [-50,150]
        eflux = pf[*,0]*cmap0
        
        store_data, the_var, times, eflux, limits={$
            ytitle: '@100km!C(mW/m!U2!N)', $
            yrange: yrange, $
            ystyle: 1, $
            yticks: 2, $
            yminor: 4, $
            labels: 'S!D||!N!C  Poynting!C  flux'}
    endforeach

;---Plot the energy fluxes.
    var_suffix = ['o_pa','o_en','ion_eflux','pf_para']
    nypanel = n_elements(var_suffix)
    nxpanel = n_elements(probes)
    fig_xsize = 8
    fig_ysize = 4
    xticklen_chsz = -0.15
    yticklen_chsz = -0.30
    plot_file = join_path([srootdir(),'fig_hope.pdf'])
    if keyword_set(test) then plot_file = test
    margins = [10,4,8,1]
    sgopen, plot_file, xsize=fig_xsize, ysize=fig_ysize, /inch
    poss = sgcalcpos(nypanel,nxpanel, xpad=18, ypad=0.4, xchsz=xchsz, ychsz=ychsz, margins=margins)
    tops = poss[*,0,0]
    xticklen = xticklen_chsz*ychsz/(tpos[3]-tpos[1])
    yticklen = yticklen_chsz*xchsz/(tpos[2]-tpos[0])
    foreach probe, probes do begin
        pre0 = 'rb'+probe+'_'
        options, pre0+var_suffix, 'xticklen', xticklen
        options, pre0+var_suffix, 'yticklen', yticklen
    endforeach

    foreach probe, probes do begin
        pre0 = 'rb'+probe+'_'

        options, pre0+'o_pa', 'ytitle', 'Pitch!C(deg)'
        options, pre0+'o_en', 'ytitle', 'Energy!C(eV)'

        the_poss = (probe eq 'a')? poss[*,0,*]: poss[*,-1,*]
        the_poss = reform(the_poss)
        loadct2, 43
        device, decomposed=0
        tplot, pre0+var_suffix[0:-3], position=the_poss[*,0:-3], /novtitle, /nouttick, /noerase
        device, decomposed=1
        tplot, pre0+var_suffix[-2:*], position=the_poss[*,-2:*], /noerase
        
        for ii=nypanel-3, nypanel-1 do begin
            tpos = the_poss[*,ii]
            tx = tpos[0]+xchsz*0.5
            ty = (ii eq nypanel-1)? tpos[3]-ychsz*1: tpos[1]+ychsz*0.3
            msg = (ii eq nypanel-1)? 'Earthward': 'Outflowing'
            xyouts, tx,ty,msg, /normal
        endfor

        for ii=0, nypanel-1 do begin
            tpos = the_poss[*,ii]
            tx = tpos[0]-xchsz*8
            ty = tpos[3]-ychsz*constant('full_ychsz')
            xyouts, tx,ty,/normal, probe+'-'+string(ii+1,format='(I0)')+'.'
        endfor
    endforeach

    if keyword_set(test) then stop
    sgclose

end
