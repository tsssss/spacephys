
pro fig_poynt, ifn, tr, filename = ofn, ps = ps, png = png

    labfac = ['v','p','b']
    rgb = [6,4,2]
    !p.font = 0 & !y.charsize = 0.7 & !x.charsize = 0.7 &!z.charsize = 0.6
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 5
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.7
    time_stamp, /off
    
    ifnfa = sdiskdir('Works')+'/works/polarcap/data/fa_sdt_fld_19980925_08278.tplot'
    trfa = time_double(['1998-09-25/04:27','1998-09-25/04:32'])
    posfa = [0.55,0.15,0.95,0.95]
    ifnpo = sdiskdir('Works')+'/works/polarcap/data/po_sdt_fld_1998092415.sdt'
    trpo = time_double(['1998-09-25/05:00','1998-09-25/06:00'])
    pospo = [0.1,0.15,0.5,0.95]
    ofn = shomedir()+'/field.'
    if keyword_set(ps) then ofn += 'eps'
    if keyword_set(png) then ofn += 'png'

    sgwindow, 0, xsize = 7, ysize = 4, /inch
    if keyword_set(ps) then sgpsopen, ofn
    if keyword_set(png) then sgpsopen, ofn
    sgindexcolor, ct = 43
    
    ; plot fast.
    fast_sdt_prep_poynting_flux, ifnfa
    stplot_calc_poynting_flux, 'de_facv', 'db_facv', 'pf_facv', $
        method = 'mat', filter = [3,10,30,100], scaleinfo = [0.5,250,60]
    
    vars1 = 'de_facv_mat_comp'+labfac
    vars2 = 'db_facv_mat_comp'+labfac
    stplot_split, 'de_facv_mat', newnames = vars1, $
        labels = 'dE'+labfac, ytitles = 'dE'+labfac+' (mV/m)'
    stplot_split, 'db_facv_mat', newnames = vars2, $
        labels = 'dB'+labfac, ytitles = 'dB'+labfac+' (nT)'
    vars = [vars1,vars2]
    nvar = n_elements(vars)
    
    pos = transpose(sgcalcpos(nvar, position=posfa))
    
    options, vars1, 'ytitle', '(mV/m)'
    options, vars2, 'ytitle', '(nT)'
    options, 'mlt', 'ytitle', 'MLT (hr)'
    options, 'ilat', 'ytitle', 'ILat (deg)'
    options, 'dis', 'ytitle', 'Dist (Re)'
    ylim, vars1, -200, 200, 0
    ylim, vars2, -400, 400, 0
    options, vars, 'yticks', 2
    options, vars, 'yminor', 4
    tplot, vars, trange = trfa, var_label = ['ilat','mlt','dis'], position = pos
    timebar, time_double('1998-09-25/04:27:55'), thick = 2
    timebar, time_double('1998-09-25/04:29:35'), thick = 2
    
    ; plot polar.
    
    if keyword_set(ps) then begin sgpsclose & wdelete, 0 & endif
    if keyword_set(png) then begin sgzclose & wdelete, 0 & endif

end

fig_poynt, /ps
fig_poynt
end