
pro test_sg_basic_plot_structure

    fn = sdiskdir('Works')+'/confs/seminar/seminar_2014_1021/data/seminar_dat.tplot'
    tplot_restore, filename = fn
    ofn = shomedir()+'/agu_2014_fig_eflux.eps'
    eventid = '1998_1001_02'
    logfile = sdiskdir('Works')+'/works/cusp/cusp_list_of_conjun.log'
    info = cusp_read_conjun_list(logfile, event = eventid)
    potr = info.polar.plot_time
    fatr = info.fast.plot_time
    potrcusp = info.polar.cusp_time
    fatrcusp = info.fast.cusp_time
    
    ; **** plot settings.
    posl = [0.1,0.15,0.45,0.95]
    posr = [0.6,0.15,0.95,0.95]
    faclabs = ['v','p','b']
    rgb = [6,4,2]
    !p.font = 1 & !y.charsize = 0.7 & !x.charsize = 0.7 &!z.charsize = 0.6
    thick = 2 & !p.thick = 2 & !x.thick = 2 & !y.thick = 2 & !z.thick = 2
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.7
    time_stamp, /off
    
    options, ['po_','fa_']+'ele_keflux_map', 'labels', 'KEe'
    options, ['po_','fa_']+'ion_keflux_map', 'labels', 'KEi'
    options, ['po_','fa_']+'pf_fac_mat_para_map', 'labels', 'S!D//!N'
    povars = 'po_'+['ele_keflux','ion_keflux','pf_fac_mat_para']+'_map'
    favars = 'fa_'+['ele_keflux','ion_keflux','pf_fac_mat_para']+'_map'
    options, [povars,favars], 'ytitle', '(mW/m!U2!N)'
    ylim, povars, -10, 15
    ylim, favars, -10, 15
    
    sgpsopen, ofn, xsize = 7, ysize = 5, /inch
    sgindexcolor, ct = 43
    ; polar.
    vars = povars
    labs = 'po_'+['ilat','mlt','dis']
    nvar = n_elements(vars)
    pos = transpose(sgcalcpos(nvar, position=posl))
    tplot, vars, trange = potr, position = pos, /noerase, var_label = labs, title = 'Polar (P)'
    timebar, potrcusp, thick = thick
    ; fast.
    vars = favars
    labs = 'fa_'+['ilat','mlt','dis']
    nvar = n_elements(vars)
    pos = transpose(sgcalcpos(nvar, position=posr))
    tplot, vars, trange = fatr, position = pos, /noerase, var_label = labs, title = 'FAST (F1)'
    timebar, fatrcusp, thick = thick
    sgpsclose, /pdf

end