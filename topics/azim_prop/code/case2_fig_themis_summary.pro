;+
; Plot the summary plot for THD.
;-


;---Input.
    time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    test_times = time_double(['2014-08-28/10:10','2014-08-28/10:20'])
test = 0


    foreach the_probe, ['a','d','e'] do begin
        probe = 'th'+the_probe
        prefix = probe+'_'


    ;---Load data.
        _2014_0828_10_load_data
        en_spec_var = prefix+'ele_en_spec'
        if tnames(en_spec_var) eq '' then $
            themis_read_ele_en_spec, time_range, probe=the_probe
        en_spec_var = prefix+'ion_en_spec'
        if tnames(en_spec_var) eq '' then $
            themis_read_ion_en_spec, time_range, probe=the_probe
        n_var = prefix+'ele_n'
        if tnames(n_var) eq '' then $
            themis_read_density, time_range, probe=the_probe
        u_var = prefix+'u_gsm'
        if tnames(u_var) eq '' then $
            themis_read_ion_vel, time_range, probe=the_probe

        ; Load tilt.
        project = azim_df_load_project()
        azim_df_load_basic_data, project=project


    ;---Settings.
        options, n_var, 'yrange', [0.03,3]
        foreach var, prefix+'kev_'+['e','h']+'_flux' do begin
            get_data, var, times, fluxes, ebins, limits=lim
            ;nebin = n_elements(ebins)
            step = (var eq prefix+'kev_e_flux')? 3: 3
            fluxes = fluxes[*,0:*:step]
            ebins = ebins[0:*:step]
            store_data, var+'_plot', times, fluxes, ebins
            yrange = 10d^ceil(alog10(minmax(fluxes)))>1
            yrange = (var eq prefix+'kev_e_flux')? [1,1e7]: [1,10e5]
            ytickv_log = make_bins(alog10(yrange),1)
            ytickv = 10d^ytickv_log
            yticks = n_elements(ytickv)-1
            yminor = 10
            ytickn = '10!U'+string(ytickv_log,format='(I0)')
            ytickn[0:*:2] = ' '
            add_setting, var+'_plot', /smart, {$
                display_type: 'list', $
                ylog:1, $
                yrange: yrange, $
                yticks: yticks, $
                ytickv: ytickv, $
                yminor: yminor, $
                ytickname: ytickn, $
                color_table: 52, $
                unit: '#/cm!U2!N-s-sr-keV', $
                value_unit: 'keV', $
                short_name: 'e!U-!N flux'}
            if var eq prefix+'kev_h_flux' then options, var+'_plot', 'ytitle', 'H!U+!N flux!C(#/cm!U2!N-s-sr-keV)'
        endforeach
        options, prefix+['ele','ion']+'_en_spec', 'ytitle', 'Energy!C(eV)'

        the_var = prefix+'theta'
        options, the_var, 'yrange', [-1,1]*25
        options, the_var, 'yticks', 2
        options, the_var, 'yminor', 5

        the_var = prefix+'b_gsm'
        get_data, the_var, times, b_gsm
        yrange = minmax(b_gsm)/10.
        yrange = [floor(yrange[0]),ceil(yrange[1])]*10
        options, the_var, 'yrange', yrange
        options, the_var, 'yticks', 2
        options, the_var, 'yminor', 5

        the_var = prefix+'e_gsm'
        options, the_var, 'yrange', [-1,1]*30
        options, the_var, 'yticks', 2
        options, the_var, 'yminor', 6

        the_var = prefix+'theta'
        options, the_var, 'constant', 0
        options, the_var, 'labels', 'Detrended!C  tile angle'


        ; Pflux.
        get_data, prefix+'pf_fac', times, pffac
        cmap = get_var_data(prefix+'cmap', at=times)
        pfpara = pffac[*,0]*cmap
        tvar = prefix+'pf_para'
        store_data, tvar, times, pfpara, limits={$
            ytitle: '(mW/m!U2!N)', constant: 0, yticks:2, yminor:5, labels:'@100 km'}
        yrange = minmax(pfpara)*0.1
        yrange = [floor(yrange[0]),ceil(yrange[1])]*10
        options, tvar, 'yrange', yrange


        ; Earthward vel.
        get_data, prefix+'u_gsm', times, u_gsm
        r_gsm = get_var_data(prefix+'r_gsm', at=times)
        ur = -vec_dot(u_gsm, sunitvec(r_gsm))
        yrange = minmax(ur)/100
        yrange = [floor(yrange[0]),ceil(yrange[1])]*100
        store_data, prefix+'ur', times, ur, limits={$
            ytitle: '(km/s)', yrange:yrange, $
            yticks:2, yminor:5, constant:0, labels: 'Earthward!C  component'}




    ;---Plot.
        root_dir = join_path([homedir(),'Dropbox','mypapers','df_substorm','plot'])
        ofn = join_path([root_dir,'case2_fig_'+probe+'_summary.pdf'])
        if keyword_set(test) then ofn = 0
        sgopen, 0, xsize=1, ysize=1, xchsz=abs_xchsz, ychsz=abs_ychsz
        sgclose, /wdelete
        pan_xsize = 4
        pan_ysize = 1
        ypad = 0.4
        margins = [10,4,8,2]
        vars = prefix+['ele_en_spec','ion_en_spec','ele_n','b_gsm','e_gsm','pf_para','theta','ur', $
            'kev_e_flux_plot','kev_h_flux_plot']
        nvar = n_elements(vars)
        labels = letters(nvar)+' .'+['Ele','Ion','Density','B GSM','E GSM','S!D||!N','Tilt', 'Ion V','','']
;        if probe ne 'thd' then begin
;            vars = vars[0:-3]
;            nvar -= 2
;            labels = labels[0:-3]
;        endif
        fig_xsize = pan_xsize+total(margins[[0,2]])*abs_xchsz
        fig_ysize = pan_ysize*nvar+(ypad*(nvar-1)+total(margins[[1,3]]))*abs_ychsz

        sgopen, ofn, xsize=fig_xsize, ysize=fig_ysize, xchsz=xchsz, ychsz=ychsz
        tplot_options, 'version', 2

        poss = sgcalcpos(nvar, margins=margins, ypad=ypad)
        tplot, vars, trange=time_range, position=poss
        for ii=0, nvar-1 do begin
            tpos = poss[*,ii]
            tx = xchsz*1
            ty = tpos[3]-ychsz*0.8
            xyouts, tx,ty,/normal, labels[ii]
        endfor
        timebar, test_times, linestyle=1

        tpos = poss[*,0]
        tx = total(tpos[[0,2]])*0.5
        ty = tpos[3]+ychsz*0.5
        xyouts, tx,ty,/normal, alignment=0.5, strupcase(probe)+' summary plot'

        if keyword_set(test) then stop
        sgclose
    endforeach

end
