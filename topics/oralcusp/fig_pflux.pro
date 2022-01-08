
pro fig_pflux, fn0, ps = ps, png = png, reload = reload

    if n_elements(fn0) eq 0 then fn0 = dialog_pickfile()
    if n_elements(tnames()) eq 1 then tplot_restore, filename = fn0
    if keyword_set(reload) then tplot_restore, filename = fn0
    ofn = shomedir()+'/pflux.'
    if keyword_set(ps) then ofn += 'eps'
    if keyword_set(png) then ofn += 'png'
    
    ; default settings.
    labfac = ['v','p','b']
    rgb = [6,4,2]
    !p.font = 1 & !y.charsize = 0.7 & !x.charsize = 0.7 &!z.charsize = 0.6
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.7
    time_stamp, /off
    
    posr = [0.6,0.15,0.95,0.95]
    posl = [0.1,0.15,0.45,0.95]
    prepo = 'po_' & prefa = 'fa_'
    trfa = time_double(['1998-09-25/04:27','1998-09-25/04:32'])
    trpo = time_double(['1998-09-25/05:00','1998-09-25/06:00'])
    
    sgwindow, 0, xsize = 7, ysize = 5, /inch
    if keyword_set(ps) then sgpsopen, ofn
    if keyword_set(png) then sgpsopen, ofn
    sgindexcolor, ct = 43
    xyouts, sgfont('','times')
    
    ; **** polar.
    pre = prepo
    ftrs = ['!C'+['10-1500','10-40','40-180','180-1500']+' sec']
    lbls = ['total','high','mid','low']
    vars = [pre+'pf_fac_mat',tnames(pre+'pf_fac_matf?')]
    nvar = n_elements(vars)
    for i = 0, nvar-1 do $
        options, vars[i], 'ytitle', 'S fac!C(mW/m!U2!N)'+ftrs[i]
    pos = transpose(sgcalcpos(6, position = posl))
    
    options, vars, 'labels', ['v','p','b']
    options, vars, 'yticks', 3
    options, vars, 'yminor', 5
    ylim, vars, -0.5, 1.0, 0
    
    labs = pre+['ilat','mlt','dis']
    options, pre+'mlt', 'ytitle', 'MLT (hr)'
    options, pre+'ilat', 'ytitle', 'ILat (deg)'
    options, pre+'dis', 'ytitle', 'Dist (Re)'
    
    tplot, vars, trange = trpo, var_label = labs, position = pos
    timebar, time_double('1998-09-25/05:27'), thick = 2
    timebar, time_double('1998-09-25/05:40'), thick = 2
    
    ; **** fast.
    pre = prefa
    ftrs = ['!C'+['3-100','3-10','10-30','30-100']+' sec']
    lbls = ['total','high','mid','low']
    vars = [pre+'pf_fac_mat',tnames(pre+'pf_fac_matf?')]
    nvar = n_elements(vars)
    for i = 0, nvar-1 do $
        options, vars[i], 'ytitle', 'S fac!C(mW/m!U2!N)'+ftrs[i]
    pos = transpose(sgcalcpos(6, position = posr))
    
    options, vars, 'labels', ['v','p','b']
    options, vars, 'yticks', 3
    options, vars, 'yminor', 5
    ylim, vars, -5, 10, 0
    
    labs = pre+['ilat','mlt','dis']
    options, pre+'mlt', 'ytitle', 'MLT (hr)'
    options, pre+'ilat', 'ytitle', 'ILat (deg)'
    options, pre+'dis', 'ytitle', 'Dist (Re)'
    
    tplot, vars, trange = trfa, var_label = labs, position = pos, /noerase
    timebar, time_double('1998-09-25/04:27:55'), thick = 2
    timebar, time_double('1998-09-25/04:29:35'), thick = 2
    
    ; add title.
    xyouts, 0.75, 0.97, '(b) FAST', /normal, charsize = !p.charsize, alignment = 0.5
    xyouts, 0.25, 0.97, '(a) Polar', /normal, charsize = !p.charsize, alignment = 0.5

    if keyword_set(ps) then begin sgpsclose & wdelete, 0 & endif
    if keyword_set(png) then begin sgzclose & wdelete, 0 & endif
end

case susrhost() of
    'Sheng@Xps': fn0 = shomedir()+'/Downloads/oral_cusp_data.tplot'
endcase
fig_pflux, fn0, /ps
fig_pflux, fn0
end