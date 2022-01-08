
pro fig_esa, ifn, tr, filename = ofn

    if n_elements(ofn) eq 0 then message, 'no output filename ...'
    tplot_restore, filename = ifn
    

    options, 'ele_nflux', 'labels', 'Je'
    options, 'ele_eflux', 'labels', 'KEe'
    
    options, 'ion_en_spec', 'ytitle', 'Energy!C(eV)'
    options, 'ele_en_spec', 'ytitle', 'Energy!C(eV)'
    
    options, 'ion_pa_spec', 'ytitle', 'Pitch!CAngle!C(deg)'
    options, 'ion_pa_spec', 'yrange', [-90,270]
    options, 'ion_pa_spec', 'ystyle', 1
    options, 'ele_pa_spec', 'ytitle', 'Pitch!CAngle!C(deg)'
    options, 'ele_pa_spec', 'yrange', [-90,270]
    options, 'ele_pa_spec', 'ystyle', 1
    options, 'ele_pa_spec', 'yticks', 2
    options, 'ele_pa_spac', 'yminor', 2
    
    options, 'ion_nflux', 'labels', 'Ji'    
    options, 'ion_nflux', 'yrange', [-8e9,0]
    options, 'ion_nflux', 'ystyle', 1
    options, 'ion_nflux', 'yticks', 2
    options, 'ion_nflux', 'ytitle', 'Ji!C(1/cm!U2!N-s)'
    options, 'ele_nflux', 'labels', 'Je'    
    options, 'ele_nflux', 'yticks', 3
    options, 'ele_nflux', 'ytitle', 'Je!C(1/cm!U2!N-s)'
    
    options, 'ion_eflux', 'labels', 'KEi'
    options, 'ion_eflux', 'ytitle', 'KEi!C(mW/m!U2!N)!C '
    options, 'ion_eflux', 'yrange', [-1,5]
    options, 'ion_eflux', 'ystyle', 1
    options, 'ion_eflux', 'yticks', 2
    options, 'ion_eflux', 'yminor', 3
    options, 'ele_eflux', 'labels', 'KEe'
    options, 'ele_eflux', 'ytitle', 'KEe!C(mW/m!U2!N)!C '
    options, 'ele_eflux', 'yrange', [-1,5]
    options, 'ele_eflux', 'ystyle', 1
    options, 'ele_eflux', 'yticks', 2
    options, 'ele_eflux', 'colors', 6
    options, 'ele_eflux', 'yminor', 3
    
    ; default settings.
    device, decomposed = 0 & loadct2, 43
    !p.font = 1 & !p.thick = 2
    !y.charsize = 0.75 & !x.charsize = 0.75 &!z.charsize = 0.65
    time_stamp, /off
    tplot_options, 'ygap', 0.4
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 5
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.7
    
    vars = ['ion_en_spec','ion_pa_spec','ion_nflux','ion_eflux', $
        'ele_en_spec','ele_pa_spec','ele_nflux','ele_eflux']
    
    stplot_merge, ['ion_n','ele_n'], newname = 'density', $
        ytitle = 'n (cm!U-3!N)', colors = [0,6], labels = ['ni','ne']
    options, 'density', 'ylog', 1
    vars = ['ion_en_spec','ion_pa_spec','density','ion_eflux','ele_eflux']
    labs = ['ilat','mlt','dis']
    sgwindow, 0, xsize = 4, ysize = 4.28, /inch
;    sgwindow, 0, xsize = 5, ysize = 5.69, /inch
    wdelete, 0
    sgpsopen, ofn
    sgindexcolor
    tplot, vars, var_label = labs, trange = tr
    xyouts, 0.5, 0.96, '(b) FAST/ESA', /normal, charsize = !p.charsize, alignment = 0.5
    timebar, time_double('1998-09-25/04:27:55'), thick = 2
    timebar, time_double('1998-09-25/04:29:35'), thick = 2
    sgpsclose
    
end

ifn = sdiskdir('Works')+'/works/polarcap/data/fa_sdt_esa_19980925_08278.tplot'
ofn = shomedir()+'/esa.eps'
tr = time_double(['1998-09-25/04:27','1998-09-25/04:32'])
fig_esa, ifn, tr, filename = ofn
end
