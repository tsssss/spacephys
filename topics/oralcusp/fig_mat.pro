pro fig_mat, fn0, ps = ps, png = png, reload = reload

    if n_elements(fn0) eq 0 then fn0 = dialog_pickfile()
    if n_elements(tnames()) eq 1 then tplot_restore, filename = fn0
    if keyword_set(reload) then tplot_restore, filename = fn0
    ofn = shomedir()+'/mat.'
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
    pre = 'po_'
    
    sgwindow, 0, xsize = 7, ysize = 5, /inch
    if keyword_set(ps) then sgpsopen, ofn
    if keyword_set(png) then sgpsopen, ofn
    sgindexcolor, ct = 43
    xyouts, sgfont('','times')
    
    labs = pre+['ilat','mlt','dis']
    options, pre+'mlt', 'ytitle', 'MLT (hr)'
    options, pre+'ilat', 'ytitle', 'ILat (deg)'
    options, pre+'dis', 'ytitle', 'Dist (Re)'
    
    ; **** dEv.
    vars = [pre+'de_fac',pre+'de_fac_comp1_mat',pre+'de_fac_mat',$
        tnames(pre+'de_fac_matf?')]
    tmp = vars[[0,2,3,4,5]] & idx = 0
    ftrs = ['','!C'+['10-1500','10-40','40-180','180-1500']+' sec']
    lbls = ['origin','total','high','mid','low']
    for i = 0, n_elements(tmp)-1 do $
        stplot_index, tmp[i], idx, newname = tmp[i]+'_v', $
        ytitle = 'dEv!C(mV/m)'+ftrs[i], labels = lbls[i]
    vars[[0,2,3,4,5]] += '_v'
    nvar = n_elements(vars)
    pos = transpose(sgcalcpos(nvar, position=posl))
    
    vname = pre+'de_fac_comp1_mat'
    options, vname, 'ytitle', 'dEv MAT!C(sec)'
    options, vname, 'ztitle', '(mV/m)'
    options, vname, 'no_interp'
    zlim, vname, -2, 2, 0
    ylim, pre+'de_fac_v', -130, 50, 0
    
    vname = tnames(pre+'de_fac_matf?')
    ftrs = ['10-40','40-180','180-1500']+'!Csec'
    options, vname, 'ytitle', 'dEv fac!C(mV/m)'
    for i = 0, n_elements(vname)-1 do options, vname[i], 'labels', ftrs[i]

    tplot, vars, var_label = labs, position = pos
    timebar, time_double('1998-09-25/05:27'), thick = 2
    timebar, time_double('1998-09-25/05:40'), thick = 2
    ftrs = [10,40,180,1500]
    for i = 0, n_elements(ftrs)-1 do $
        stplot_bar, pre+'db_fac_comp1_mat', ftrs[i], position = pos[*,1]
    
    ; **** dBp.
    vars = [pre+'db_fac',pre+'db_fac_comp2_mat',pre+'db_fac_mat',$
        tnames(pre+'db_fac_matf?')]
    tmp = vars[[0,2,3,4,5]] & idx = 1
    ftrs = ['','!C'+['10-1500','10-40','40-180','180-1500']+' sec']
    lbls = ['origin','total','high','mid','low']
    
;    vars = [pre+'db_fac_mat',pre+'db_fac_comp2_mat',$
;        tnames(pre+'db_fac_matf?')]
;    tmp = vars[[0,2,3,4]] & idx = 1
;    ftrs = ['!C'+['10-1500','10-40','40-180','180-1500']+'sec']
;    lbls = ['total','high','mid','low']
    for i = 0, n_elements(tmp)-1 do $
        stplot_index, tmp[i], idx, newname = tmp[i]+'_p', $
        ytitle = 'dBp!C(nT)'+ftrs[i], labels = lbls[i]
    vars[[0,2,3,4,5]] += '_p'
    nvar = n_elements(vars)
    pos = transpose(sgcalcpos(nvar, position=posr))
    
    vname = pre+'db_fac_comp2_mat'
    options, vname, 'ytitle', 'dBp MAT!C(nT)'
    options, vname, 'ztitle', '(nT)'
    options, vname, 'no_interp'
    zlim, vname, -6, 6, 0
    ylim, pre+'db_fac_p', -100, 80, 0
   
    tplot, vars, var_label = labs, position = pos, /noerase
    timebar, time_double('1998-09-25/05:27'), thick = 2
    timebar, time_double('1998-09-25/05:40'), thick = 2
    ftrs = [10,40,180,1500]
    for i = 0, n_elements(ftrs)-1 do $
        stplot_bar, pre+'db_fac_comp1_mat', ftrs[i], position = pos[*,1]
    
    ; add title.
    xyouts, 0.25, 0.97, '(b) Polar dEv', /normal, charsize = !p.charsize, alignment = 0.5
    xyouts, 0.75, 0.97, '(a) Polar dBp', /normal, charsize = !p.charsize, alignment = 0.5
    
    if keyword_set(ps) then begin sgpsclose & wdelete, 0 & endif
    if keyword_set(png) then begin sgzclose & wdelete, 0 & endif
end

case susrhost() of
    'Sheng@Xps': fn0 = shomedir()+'/Downloads/oral_cusp_data.tplot'
endcase
fig_mat, fn0, /ps
end